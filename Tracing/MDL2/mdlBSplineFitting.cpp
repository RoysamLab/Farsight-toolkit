/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include "mdlBSplineFitting.h"

namespace mdl
{

//Constructor
BSplineFitting::BSplineFitting(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	region = m_inputImage->GetBufferedRegion();
	sizeX = region.GetSize(0);
	sizeY = region.GetSize(1);
	sizeZ = region.GetSize(2);
	numPix = sizeX*sizeY*sizeZ;

	debug = false;
	splineOrder = 3;
	splineLevels = 7;

	//input
	nodes = NULL;
	bbpairs = NULL;
}

BSplineFitting::~BSplineFitting()
{
	m_inputImage=NULL;
}

bool BSplineFitting::Update()
{
	if(!m_inputImage || !nodes || !bbpairs)
		return false;

	this->findBranches();
	this->smoothBranches();
	this->detectExtraSpines();
	
	//Write out smooth backbone and extra spines files:
	if(debug)
	{
		std::cerr << "DONE\n";
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes_out);
		fhdl->SetLines(&bbpairs_out);
		fhdl->Write("SmoothBackbone.vtk");
		delete fhdl;
	}//end if debug

	return true;
}

void BSplineFitting::findBranches()
{
	if(!nodes || !bbpairs)
		return;

	int num_nodes = (int)nodes->size();

	//Graph properties for each node: degree and connectivity
	std::vector<Graphprop> graphPointInfo;
	//Initialize the graph containers (for each node):
	for(int i=0; i<num_nodes; i++)
    {
		Graphprop p;
		p.deg = 0;
		p.outVert.clear();
		graphPointInfo.push_back(p);
    }

    //Iterate through bbpairs and populate the graph info
	//we need the degree of each node and a list of up to 8 connected nodes
	for(int i=0; i<(int)bbpairs->size(); ++i)
	{
		int nodeIndex1 = bbpairs->at(i).first;
		int nodeIndex2 = bbpairs->at(i).second;

		int tmpdeg = graphPointInfo.at(nodeIndex1).deg;
		graphPointInfo.at(nodeIndex1).deg = tmpdeg+1;
		graphPointInfo.at(nodeIndex1).outVert.insert(nodeIndex2);

		tmpdeg = graphPointInfo.at(nodeIndex2).deg;
		graphPointInfo.at(nodeIndex2).deg = tmpdeg+1;
		graphPointInfo.at(nodeIndex2).outVert.insert(nodeIndex1);
	}
   
	//Create a set of branch points
	branchPts.clear();
   
	//We need to get rid of junction points 
	//so that we can break up backbone into segments
	//This is important for segments between junction points.
	//We'll save the removed nodes for later use!
	for (int i=0;i<num_nodes;i++)
	{
		if (graphPointInfo.at(i).deg >=3)
		{
			branchPts.insert( i );	//Save the node
			//iterate through each neighbor
			SetType::iterator it;
			SetType bset = graphPointInfo.at(i).outVert;
			for (it = bset.begin(); it != bset.end(); it++)
			{
				int nbrPt = (*it);
				(graphPointInfo.at(nbrPt).deg)-- ; //decrement the degree
				//Romove me as a neighbor:
				graphPointInfo.at(nbrPt).outVert.erase(i);
			} // end for
		}// end if 
	}// end for


	branches.clear();
	int NumBranches = 0;

	for (int i=0; i<num_nodes; i++)		//Iterate through all nodes
    {
		if (graphPointInfo.at(i).deg == 1)	//I am endpoint
		{
			NumBranches++;		//So must be a new branch (increment counter)
			branches[NumBranches].push_back(i);
			int lastpoint = i;
			int nextpoint = *graphPointInfo.at(i).outVert.begin();
			// if not end of branch, keep going.
			while (graphPointInfo.at(nextpoint).deg == 2) 
			{
				branches[NumBranches].push_back(nextpoint);

				if (*graphPointInfo.at(nextpoint).outVert.begin() == lastpoint)
				{	//next point was were I just was, so use other neighbor
					lastpoint = nextpoint;
					nextpoint = *graphPointInfo.at(nextpoint).outVert.rbegin();
				}
				else 
				{
					lastpoint = nextpoint;
					nextpoint = *graphPointInfo.at(nextpoint).outVert.begin(); 
				}
			} // end while 

			//Moved inside if statement so I loose branch points
			//branches[NumBranches].push_back(nextpoint);
			if(graphPointInfo.at(nextpoint).deg == 1)
			{
				//Set degree = 0 so I don't start this trace again!
				graphPointInfo.at(nextpoint).deg = 0;
				branches[NumBranches].push_back(nextpoint);
			}
		}// end if
	}// end for 

	if(debug)
	{
		std::cerr << "There are " <<  NumBranches  << " branches" << std::endl;
		for(int i=1; i<=NumBranches; ++i)
		{
			//std::cerr << i << ":" << branches[i].size() << " points\n";
		}
	}
}

void BSplineFitting::smoothBranches()
{
	nodes_out.clear();
	bbpairs_out.clear();

	//Add the branch points to the output nodes list:
	std::set<int>::iterator it3;
	for(it3=branchPts.begin(); it3!=branchPts.end(); it3++)
	{
		nodes_out.push_back( nodes->at(*it3) );
	}
	
	float discretePt = 1;				//Sampling rate

	//**********************************************
	// MAIN LOOP THROUGH EACH BRANCH
	int NumBranches = (int)branches.size();
	for(int b=1; b<=NumBranches; b++)
	{
		ListType list = branches[b];	//current points list in this branch
		std::vector<fPoint3D> points;	//hold the points in each branch
		ListType::iterator it;
		for(it=list.begin(); it!=list.end(); it++)
		{
			fPoint3D point = nodes->at(*it); //Get the point
			//point.y *= -1;	//invert y coordinate
			points.push_back(point);
		}

		int inNumPts = (int)points.size();
		int outNumPts = (float)inNumPts / discretePt;
		if(inNumPts == 0 || outNumPts == 0)
			continue;

		//Do the BSpline Algorithm:
		std::vector<fPoint3D> newPts;
		newPts = bbBSplineFitting(points, outNumPts, splineOrder, splineLevels);

		//**************
		//Update output
		int startIdx = (int)nodes_out.size();
		nodes_out.insert(nodes_out.end(),newPts.begin(),newPts.end());
		int endIdx = (int)nodes_out.size()-1;

		for(int i=startIdx; i<endIdx; i++)
		{
			bbpairs_out.push_back( pairE(i,i+1) );
		}
	}	//MAIN LOOP END
}

std::vector<fPoint3D> BSplineFitting::bbBSplineFitting(
		std::vector<fPoint3D> inPts, int numOut, unsigned int order, unsigned int levels)
{
	std::vector<fPoint3D> retVect;

	const unsigned int ParametricDimension = 1;
	const unsigned int DataDimension = 3;
	typedef double RealType;
	typedef itk::Vector<RealType, DataDimension> VectorType;
	typedef itk::Image<VectorType, ParametricDimension> bImageType;  
	typedef itk::PointSet<VectorType, ParametricDimension> PointSetType;

	PointSetType::Pointer pointSet = PointSetType::New();  
	
	int numIn = (int)inPts.size();
	for(int i=0; i<numIn; ++i)
    {
		VectorType V;
		V[0] = inPts.at(i).x;
		V[1] = inPts.at(i).y;
		V[2] = inPts.at(i).z;

		RealType t = (RealType)i / (RealType)(numIn-1);
		PointSetType::PointType point;
		point[0] = t;
		pointSet->SetPoint( i, point );
		pointSet->SetPointData( i, V );
    }

	// Instantiate the filter and set the parameters
	typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, bImageType>  FilterType;
	FilterType::Pointer filter = FilterType::New();
  
	// Define the parametric domain
	bImageType::SpacingType spacing;  
	spacing.Fill( 0.002 ); // 0.001
	bImageType::SizeType size;  
	size.Fill( static_cast<unsigned int>( 1.0/spacing[0] ) + 1 );
	bImageType::PointType origin;  
	origin.Fill( 0.002 );// 0.0

	filter->SetSize( size );
	filter->SetOrigin( origin );
	filter->SetSpacing( spacing );
	filter->SetInput( pointSet );

	filter->SetSplineOrder( order );  
	FilterType::ArrayType ncps;
	ncps.Fill( 4 );  
	filter->SetNumberOfControlPoints( ncps );
	filter->SetNumberOfLevels( levels );
	filter->SetGenerateOutputImage( false );
	try 
	{
		filter->Update();
    }
	catch (...) 
    {
		std::cerr << "Test 2: itkBSplineScatteredDataImageFilter exception thrown" << std::endl;
		return retVect;
    }
  
	for(int i=0; i<numOut; ++i)
    {
		RealType t = (RealType)i / (RealType)(numOut-1);
		PointSetType::PointType point;
		point[0] = t;

		VectorType V; 
		filter->Evaluate( point, V );
		if(abs(V[0])>0.1 && abs(V[1]) > 0.1 && abs(V[2]) > 0.1)
		{
			fPoint3D point;
			point.x = V[0];
			point.y = V[1];
			point.z = V[2];
			retVect.push_back(point);
		}
	}
	
	return retVect;

}

void BSplineFitting::detectExtraSpines()
{
	/*
	for(LOOP THROUGH EACH BRANCH????)
	{
	//- Compute missing spines from spline function and original 3D points--//
		
		//find minimum x
		int mintemp = points.at(0).x;
		for (int i=0; i<inNumPts; i++)
		{	 
		   if (points.at(i).x < mintemp)
			   mintemp = points.at(i).x;
		}

		float * Min_dist2spline = new float[inNumPts];
		int * Index_Min_dist2spline = new int[inNumPts];
		if(!Min_dist2spline || !Index_Min_dist2spline)
			return false;
		
		for (int i=0; i<inNumPts; i++)
		{
			int x_val_index = selffloor( (points[i].x - mintemp) / discretePt );
			if(x_val_index <= 0)       
				x_val_index=0;
			if(x_val_index >= inNumPts-1)
				x_val_index=inNumPts-1; 

			//Initialize the min dist and corresponding index with a parallel point on the spline function
			float minDist = (float)dist2pts(points[i], newPts[x_val_index]);
			int minIndex = x_val_index;
			for( int j=0; j<outNumPts; j++)
			{  
				//If the point on spline falls within the local range
				if(abs(points[i].x - newPts[j].x) < minDist)
				{
					//evaluate the current distance of two points
					float dis = (float) dist2pts(points[i], newPts[j]);
					if (dis < minDist)
					{
						minDist = dis;
						minIndex = j;
					}
				}
			}
			Min_dist2spline[i] = minDist;
			Index_Min_dist2spline[i] = minIndex;
		}// end for

		//--- Find all the spines with local max distance-----//  
		for (int i=0; i<inNumPts; i++)
		{
			int IsLocalMax = 1;
			for(int j=0; j<inNumPts; j++)
			{
				if (dist2pts(points[j], points[i]) < 5)
				{
					//Consider only neighbor points
					if (Min_dist2spline[j] > Min_dist2spline[i])  
					{
						IsLocalMax=0;  
						break; 
					}
				}
			}// end subfor

			//Get intensity values from image:
			ImageType::IndexType index1;
			index1[0] = points[i].x;
			index1[1] = points[i].y;
			index1[2] = points[i].z;
			ImageType::IndexType index2;
			long idx=Index_Min_dist2spline[i];
			index2[0] = newPts[idx].x;
			index2[1] = newPts[idx].y;
			index2[2] = newPts[idx].z;
			PixelType pix1 = m_inputImage->GetPixel(index1);
			PixelType pix2 = m_inputImage->GetPixel(index2); //I got negative here

			//Average intensity:
			int aveInt = (pix1 + pix2)/2; 
      
			// 2, 1.5 - threshold for the detect of missing spines  
			if ( IsLocalMax == 1 && Min_dist2spline[i]>= 1.3 && aveInt > 20)
			{
				newSPNpts.push_back(points.at(i));
				int idx = Index_Min_dist2spline[i];
				newSPNpts.push_back(newPts.at(idx));
			} // end if IsLocalMax ...
        
		} //end for i ... NumPoints

		delete[] Min_dist2spline;
		delete[] Index_Min_dist2spline;
	)
	*/
}


//---------------------Round to nearest integer function----------------------//
int BSplineFitting::round(float number)
{  
	int temp = int(number);
	if ((number-temp)>0.5) 
    {
		return (int(number+float(0.5)));
    }
	else 
    {
		return (temp);
    }
}

// y = floor(a) rounds fi object a to the nearest integer in the direction of
// negative infinity and returns the result in fi object y.
int BSplineFitting::selffloor(float number)
{
	if (number >0 || number ==0) 
    {
		return (int)number;
    }
	else
    {
		return (int)(number-0.5); 
    }
}

double BSplineFitting::dist2pts(fPoint3D p1, fPoint3D p2)
{
	double h;
	double dx = (double)(p1.x - p2.x);
	double dy = (double)(p1.y - p2.y);
	double dz = (double)(p1.z - p2.z);
	h=sqrt((dx*dx)+(dy*dy)+(dz*dz));
	return h;
}


}

