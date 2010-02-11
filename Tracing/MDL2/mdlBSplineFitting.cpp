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
	splineLevels = 8;

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

	int num_nodes = (int)nodes->size();

	//Note: indexOnBackbone in format xx.ooo where xx is index of backbone, ooo is order of pts.
	//Note: branch numbering starts at 1
	//DUDE: what if more than 1000 points are on a branch????
	std::vector<float> indexOnBackbone;
	std::vector<Graphprop> graphPointInfo;

	//Initialize the graph containers (for each node):
	for(int i=0; i<num_nodes; i++)
    {
		Graphprop p;
		p.deg = 0;
		for(int j=0; j<8; ++j)
			p.outVert[j] = 0;
		graphPointInfo.push_back(p);
		indexOnBackbone.push_back( 0.0 );
    }
	//Iterate through bbpairs and populate the graph info
	//we need the degree of each node and a list of up to 8 connected nodes
	for(int i=0; i<(int)bbpairs->size(); ++i)
	{
		int nodeIndex1 = bbpairs->at(i).first;
		int nodeIndex2 = bbpairs->at(i).second;

		int tmpdeg = graphPointInfo.at(nodeIndex1).deg;
		graphPointInfo.at(nodeIndex1).deg = tmpdeg+1;
		int mod = tmpdeg % 8;
		graphPointInfo.at(nodeIndex1).outVert[mod] = nodeIndex2;

		tmpdeg = graphPointInfo.at(nodeIndex2).deg;
		graphPointInfo.at(nodeIndex2).deg = tmpdeg+1;
		mod = tmpdeg % 8;
		graphPointInfo.at(nodeIndex2).outVert[mod] = nodeIndex1;
	}

	//-----Get rid of junction point so that breaking up backbone segments------//
	for (int i=0; i<num_nodes; i++)
    {
		if(graphPointInfo.at(i).deg >= 3)	//I am a branch point
		{
			//Iterate through each branch
			for (int j=0; j<graphPointInfo.at(i).deg; j++)
			{
				int outpt = graphPointInfo.at(i).outVert[j]; //The current neighbor
				graphPointInfo.at(outpt).deg--;	//decr its count
				if ( graphPointInfo.at(outpt).deg == 1) //If neighbor is down to 1
				{
					if (graphPointInfo.at(outpt).outVert[0] == i) //And I am that neighbor
					{
						//Remove junction point from neighbor points
						graphPointInfo.at(outpt).outVert[0] = graphPointInfo.at(outpt).outVert[1]; 
					} //end if
				} //end if
			} // end for
		} // end if
    }// end for

	//This is the number of points on the branch
	std::map<int,int> NumPtsOnBranch;
	int NumBranches = 0;

	for (int i=0; i<num_nodes; i++)		//Iterate through all nodes
    {
		if (graphPointInfo.at(i).deg == 1)	//If deg = 1, I am endpoint
		{
			NumBranches++;		//So must be a new branch (increment counter)
			float indexBBpts = 1.0; //counts points on this branch
			indexOnBackbone.at(i) = NumBranches + indexBBpts/1000;
			NumPtsOnBranch[NumBranches]++;

			int lastpoint = i;
			int nextpoint = graphPointInfo.at(i).outVert[0]; // the first is [0]
			// if not end of branch, keep going.
			while (graphPointInfo.at(nextpoint).deg == 2) 
			{
				indexBBpts++;	//new point on branch
				indexOnBackbone.at(nextpoint) = NumBranches + indexBBpts/1000;
				NumPtsOnBranch[NumBranches]++;

				if (graphPointInfo.at(nextpoint).outVert[0] == lastpoint)
				{
					lastpoint = nextpoint;
					nextpoint = graphPointInfo.at(nextpoint).outVert[1];
				}
				else 
				{
					lastpoint = nextpoint;
					nextpoint = graphPointInfo.at(nextpoint).outVert[0]; 
				}
			} // end while 

			//Look for end of branch:
			if (graphPointInfo.at(nextpoint).deg == 1)
			{
				graphPointInfo.at(nextpoint).deg = 0;	//disconnect ??
				indexBBpts = indexBBpts + 1;
				indexOnBackbone.at(nextpoint) = NumBranches + indexBBpts/1000;
				NumPtsOnBranch[NumBranches]++;
			} // end if
		}// end if
	}// end for 
	//---Done assigning index to backones--------------//

	if(debug)
		std::cerr << "There are " <<  NumBranches  << " branches" << std::endl; 
	
	//**********************************************************************
	//**********************************************************************
	//BEGIN CURVE FITTING:

	//The number points we will downsample to:
	float discretePt = 1; // sampling rate
	int NewNumAllPoints = (int)( (float)num_nodes/discretePt ); 

	//This records the new points:
	//Add 1 for the end of array flag
	Point3D * pointsVTK = new Point3D[NewNumAllPoints+1];		//New bb points
	Point3D * posExtraSpines = new Point3D[NewNumAllPoints+1];	//New spine points
	if(!pointsVTK || !posExtraSpines)
		return false;

	//initialization 
	for (int i=0; i<NewNumAllPoints; i++)
    {
		pointsVTK[i].x = 0;
		pointsVTK[i].y = 0;
		pointsVTK[i].z = 0;
		posExtraSpines[i].x = -9999;	//End of array
		posExtraSpines[i].y = -9999;
		posExtraSpines[i].z = -9999;
    }
	pointsVTK[NewNumAllPoints].x = -9999;//end flags
	pointsVTK[NewNumAllPoints].y = -9999;
	pointsVTK[NewNumAllPoints].z = -9999; 

	int numVTKpts = 0;
	int numExtraSpines = 0;

	//**********************************************
	// MAIN LOOP THROUGH EACH BRANCH
	for(int b=1; b<=NumBranches; b++)
	{
		//This will be used to hold the points in each branch
		int ptsCount = NumPtsOnBranch[b];
		Point3D * points = new Point3D[ptsCount]; //num_nodes is max possible
		if(!points)
			return false;

		for(int i=0; i<num_nodes; ++i)
		{
			//If this point is on current branch
			float curIdx = indexOnBackbone.at(i);
			//If this point is on current branch
			if( selffloor(curIdx) == b ) 
			{
				int tmpindex = round((curIdx - selffloor(curIdx))*1000)-1;
				// In order to exchange x-coord and y-coord
				points[tmpindex].x = nodes->at(i).x;
				points[tmpindex].y = -1*nodes->at(i).y; //FLIP
				points[tmpindex].z = nodes->at(i).z;
			} // end if on backbone
		} //end for num_nodes

		//NumBranchPts[b] = ptsCount;			//points in this branch;
		int NumPoints = ptsCount / discretePt;	//# points I want to end up with
		
		Point3D * xyz_vals = new Point3D[NumPoints];
		float * Min_dist2spline = new float[NumPoints];
		int * Index_Min_dist2spline = new int[NumPoints];
		if(!xyz_vals || !Min_dist2spline || !Index_Min_dist2spline)
			return false;

		//Do the BSpline Algorithm:
		bool fastVersion = false;
		if(fastVersion)
		{
		}
		else
		{
			BackboneBSplineFitting(points, ptsCount, xyz_vals, NumPoints, splineOrder, splineLevels);
		}

		//Update pointsVTK
		for (int i=numVTKpts; i < numVTKpts+NumPoints; i++)
		{
			pointsVTK[i].x = xyz_vals[i-numVTKpts].x;
			pointsVTK[i].y = xyz_vals[i-numVTKpts].y;
			pointsVTK[i].z = xyz_vals[i-numVTKpts].z;
			if(i < numVTKpts+3 || i > numVTKpts+NumPoints-3)
			{
				pointsVTK[i].x = points[i-numVTKpts].x;
				pointsVTK[i].y = points[i-numVTKpts].y;
				pointsVTK[i].z = points[i-numVTKpts].z;
			}
		}

		//- Compute missing spines from spline function and original 3D points--//
		int mintemp = points[0].x;
		for (int i = 0; i<NumPoints-1; i++)
		{	 
		   if (points[i].x < mintemp)
			   mintemp = points[i].x;
		}
		
		for (int i = 0;i<NumPoints-1;i++)
		{
			int x_val_index = selffloor((points[i].x-mintemp)/discretePt);
			if(x_val_index <= 0)       
				x_val_index=0;
			if(x_val_index >= NumPoints-1)
				x_val_index=NumPoints-1; 

			//Initialize the min dist and corresponding index with a parallel point on the spline function
			float minDist = (float) dist2pts(points[i], xyz_vals[x_val_index]);
			int minIndex = x_val_index;
			for( int j = 0;j<NumPoints-1;j++)
			{  
				//If the point on spline falls within the local range
				if(abs(points[i].x - xyz_vals[j].x)< minDist)
				{
					//evaluate the current distance of two points
					float dis = (float) dist2pts(points[i], xyz_vals[j]);
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
		for (int i = 0;i<NumPoints-1;i++)
		{
			int IsLocalMax = 1;
			for(int j = 0;j<NumPoints-1;j++)
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
			index2[0] = xyz_vals[idx].x;
			index2[1] = xyz_vals[idx].y;
			index2[2] = xyz_vals[idx].z;
			PixelType pix1 = m_inputImage->GetPixel(index1);
			PixelType pix2 = m_inputImage->GetPixel(index2); //I got negative here

			//Average intensity:
			int aveInt = (pix1 + pix2)/2; 
      
			// 2, 1.5 - threshold for the detect of missing spines  
			if ( IsLocalMax == 1 && Min_dist2spline[i]>= 1.3 && aveInt > 20)
			{
				numExtraSpines = numExtraSpines + 1;
			
				posExtraSpines[numExtraSpines*2-2].x = points[i].x;
				posExtraSpines[numExtraSpines*2-2].y = points[i].y;
				posExtraSpines[numExtraSpines*2-2].z = points[i].z;
        
				int idx = Index_Min_dist2spline[i];
				posExtraSpines[numExtraSpines*2-1].x =  xyz_vals[idx].x;
				posExtraSpines[numExtraSpines*2-1].y =  xyz_vals[idx].y;
				posExtraSpines[numExtraSpines*2-1].z =  xyz_vals[idx].z;
			} // end if IsLocalMax ...
        
		} //end for i ... NumPoints
		
		numVTKpts += NumPoints;
		if(debug)
			std::cerr << "there are " << numVTKpts << "points on branch " << b << std::endl;

		delete[] xyz_vals;
		delete[] Min_dist2spline;
		delete[] Index_Min_dist2spline;
		delete[] points;

	}//end for loop through branches

	if(debug)
		std::cerr << "THE LOOP IS OVER" << std::endl;

	
	//Write out smooth backbone and extra spines files:
	if(debug)
	{
		FILE * outExspine = fopen("ExtraSpine.vtk", "w");
		if(outExspine != NULL)
		{
			//-Output the extra spines to a vtk file-//
			fprintf(outExspine, "# vtk DataFile Version 3.0\n");
			fprintf(outExspine,"MST of skel\n");
			fprintf(outExspine,"ASCII\n");
			fprintf(outExspine,"DATASET POLYDATA\n");

			fprintf(outExspine,"POINTS %d float\n", 2 * numExtraSpines);
			for (int i = 0; i<numExtraSpines*2;i++)
				fprintf(outExspine,"%d %d %d\n",posExtraSpines[i].x,-posExtraSpines[i].y,posExtraSpines[i].z);
   
			fprintf(outExspine, "LINES %d %d\n", numExtraSpines, numExtraSpines*3);
			for (int i = 0;i<numExtraSpines;i++)
				fprintf(outExspine, "2 %d %d\n", 2*i, 2*i+1); // 
			fclose(outExspine);
		}

		FILE * outbackbone = fopen("SmoothBackbone.vtk", "w");
		if (outbackbone != NULL)
		{
			//-Output the smooth Backbone to a vtk file-//
			fprintf(outbackbone, "# vtk DataFile Version 3.0\n");
			fprintf(outbackbone,"MST of skel\n");
			fprintf(outbackbone,"ASCII\n");
			fprintf(outbackbone,"DATASET POLYDATA\n");

			fprintf(outbackbone,"POINTS %d float\n", numVTKpts);
			for (int i = 0; i<numVTKpts;i++)
			{
				// - pointsVTK[i].y  is just for view
				//if (abs(pointsVTK[i].x) > 0.1 && abs(pointsVTK[i].y) > 0.1  && abs(pointsVTK[i].z) > 0.1 ) 
				fprintf(outbackbone,"%d %d %d\n",pointsVTK[i].x, -pointsVTK[i].y, pointsVTK[i].z);
			}
   
			fprintf(outbackbone, "LINES %d %d\n", numVTKpts-NumBranches, (numVTKpts-NumBranches)*3);
			int indexPts = 0;
			for (int i = 1; i<= NumBranches; i++)
			{
				for(int j =0; j < NumPtsOnBranch[i]-1; j++)
				{
					fprintf(outbackbone, "2 %d %d\n", j+indexPts, j+indexPts+1);
				}
				indexPts = indexPts + (int)NumPtsOnBranch[i];
			}
			fclose(outbackbone);
		}
	}//end if debug
	

	//MEMORY CLEAN UP

	//delete[] flagOnBackbone; //for backbone index
	//delete[] graphInfPoints;

	//delete[] numBranchpts;
	delete[] pointsVTK;
	delete[] posExtraSpines;

	return true;
}

//std::vector<Point3D> BSplineFitting::bbBSplineFitting(std::vector<Point3D

int BSplineFitting::BackboneBSplineFitting(Point3D * Points, int Num,
                           Point3D *SamplePoints, int Num2, 
						   int order, int levels)
{
	const unsigned int ParametricDimension = 1;
	const unsigned int DataDimension = 3;
	typedef double RealType;
	typedef itk::Vector<RealType, DataDimension> VectorType;
	typedef itk::Image<VectorType, ParametricDimension> bImageType;  
	typedef itk::PointSet<VectorType, ParametricDimension> PointSetType;
	PointSetType::Pointer pointSet = PointSetType::New();  
	int k=0;
	double delt = (double) 1.0/(double)Num +0.0000001;

	for ( RealType t = 0.0; t <= 1.0+1e-10; t += delt ) // <= 
 // for ( RealType t = 0.0; t < 1.0; t += delt )
    {
		unsigned long i = pointSet->GetNumberOfPoints();
		//i=Num;
		PointSetType::PointType point;
		point[0] = t;
		pointSet->SetPoint( i, point );   
		

		VectorType V;
		V[0] = Points[k].x;
		V[1] = Points[k].y;
		V[2] = Points[k].z;
		k++;
		
		pointSet->SetPointData( i, V );
    }

  // Instantiate the filter and set the parameters
	typedef itk::BSplineScatteredDataPointSetToImageFilter
						<PointSetType, bImageType>  FilterType;
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
  
	delt = (double)1.0/(double)Num2+0.000001;
	k=0;
	try 
	{
		filter->Update();
    }
	catch (...) 
    {
		std::cerr << "Test 2: itkBSplineScatteredDataImageFilter exception thrown" << std::endl;
		return EXIT_FAILURE;
    }
  
	for ( RealType t = 0.0; t <= 1.0+1e-10; t += delt) // < =
	// for ( RealType t = 0.0; t < 1.0; t += delt )
    {
		PointSetType::PointType point;
		point[0] = t;

		VectorType V; 
		filter->Evaluate( point, V );
		if(abs(V[0])>0.1 && abs(V[1]) > 0.1 && abs(V[2]) > 0.1)
		{
			SamplePoints[k].x = V[0];
			SamplePoints[k].y = V[1];
			SamplePoints[k].z = V[2];
			k++;
		}
	}
	return EXIT_SUCCESS;
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

double BSplineFitting::dist2pts(Point3D p1, Point3D p2)
{
	double h;
	double dx = (double)(p1.x - p2.x);
	double dy = (double)(p1.y - p2.y);
	double dz = (double)(p1.z - p2.z);
	h=sqrt((dx*dx)+(dy*dy)+(dz*dz));
	return h;
}

}