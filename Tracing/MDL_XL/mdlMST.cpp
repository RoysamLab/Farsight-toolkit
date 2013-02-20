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

#include "mdlMST.h"
#include "mdlVolumeProcess.h"

namespace mdl
{

//Constructor
MST::MST(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	region = m_inputImage->GetBufferedRegion();
	sizeX = region.GetSize(0);
	sizeY = region.GetSize(1);
	sizeZ = region.GetSize(2);
	numPix = sizeX*sizeY*sizeZ;

	debug = false;
	useVoxelRounding = true;
	edgeRange = 10;		//Some default values:
	power = 1;
	PruneThreshold = 4.0;
	Alpha = 0.5;

	//input
	skeletonPoints = NULL;
	nodeGraph = NULL;
	mstGraph = NULL;

	nodes.clear();
	edgeArray.clear();
	edgeWeight.clear();
	spanningTree.clear();
	nodeDegree.clear();
}



MST::~MST()
{
	m_inputImage=NULL;
	skeletonPoints=NULL;

	if(nodeGraph)
	{
		delete nodeGraph;
		nodeGraph = NULL;
	}
	if(mstGraph)
	{
		delete mstGraph;
		mstGraph = NULL;
	}
}


void MST::SetSkeletonPoints(std::vector<fPoint3D> * sp)
{
	skeletonPoints = sp;
}

// for spine

void MST::SetVesselMap(ImageType::Pointer VesselMap)
{
    m_VesselMap =  VesselMap;
	ImageType::RegionType Vesselregion;
	Vesselregion =  VesselMap->GetBufferedRegion();
    if (sizeX != Vesselregion.GetSize(0) || sizeY != Vesselregion.GetSize(1) || sizeZ != Vesselregion.GetSize(2))
	{
		m_VesselMap = m_inputImage;
	}
}

bool MST::CreateGraphAndMST(int type )
{
	if(!m_inputImage || !skeletonPoints)
		return false;

	this->skeletonPointsToNodes(useVoxelRounding);

	this->nodesToEdges( type );

	this->minimumSpanningTree();

	return true;
}

int MST::roundToInt(double v)
{
	double intpart;

	if( modf(v, &intpart) < 0.5 )
		return (int)floor(v);
	else
		return (int)ceil(v);
}

//Because the skeleton Points may be float values and may have the same
// voxel that they are closest to.
//So I make sure that each voxel is in my list only once!!!
bool MST::skeletonPointsToNodes(bool roundToNearestVoxel)
{
	if(!skeletonPoints)
		return false;

	nodes.clear();
	if(!roundToNearestVoxel)	//Just use the skeleton points
	{
		nodes.insert(nodes.begin(), skeletonPoints->begin(), skeletonPoints->end());
		if(debug)
			std::cerr << "Number of Nodes = " << nodes.size() << std::endl;
		return true;
	}

	bool * nodeAdded = new bool[numPix];
	if(!nodeAdded)	//Couldn't allocate memory
		return false;

	//- Initialize to zero
	for(long idx=0; idx<numPix; idx++)   
	{ 
		nodeAdded[idx] = false;
	}

	for(int i=0; i<(int)skeletonPoints->size(); ++i)
	{
		fPoint3D fnode = skeletonPoints->at(i);
		int x = roundToInt(fnode.x);
		int y = roundToInt(fnode.y);
		int z = roundToInt(fnode.z);
		long idx = (z)*sizeX*sizeY + (y)*sizeX + (x);
		if(idx >= numPix)
			continue;

		if(!nodeAdded[idx])
		{
			nodeAdded[idx] = true;
			fPoint3D p = {x,y,z};
			nodes.push_back( p );
		}
	}

	delete[] nodeAdded;

	if(debug)
		std::cerr << "Number of Nodes = " << nodes.size() << std::endl;

	return true;
}

//I'm going to convert the nodes into edges and edge weights.
//I only make edges between nodes that are within the edgeRange
bool MST::nodesToEdges(int type)
{
	if( (int)nodes.size() == 0)
		return false;

	int num_nodes = (int)nodes.size();

	ImageType::Pointer dtImage = NULL;
	if(type == 2)
	{
		mdl::VolumeProcess *volProc = new mdl::VolumeProcess();
		volProc->SetInput(this->m_inputImage);
		volProc->RunDistanceTransform();
		dtImage = volProc->GetOutput();
		delete volProc;
	}

    //float rate = (float) (sizeZ*sizeZ)/(sizeX*sizeY); 
	float rate = 0.1;
	//Iterate through the nodes and create edges within the specified range.
	for(int i=0; i<num_nodes; ++i)
	{
		//for(int j=1; j<i; ++j)
		for(int j=0; j<i; ++j)
		{
			fPoint3D n1 = nodes.at(i);
			fPoint3D n2 = nodes.at(j);

			if( abs(n1.x - n2.x) > (float)edgeRange )
				continue;
			if( abs(n1.y - n2.y) > (float)edgeRange )
				continue;
			if( abs(n1.z - n2.z) > (float)(edgeRange))
				continue;

			//If I'm here, then I've found a close enough node
			edgeArray.push_back( pairE(i+1,j+1) ); //add an edge (count starts at 1 for nodes)

			float weight = 10000;
			if(type == 1)
			{
				weight = getXiaosongEdgeWeight(n1, n2, power);
			}
			else if(type == 2)
			{
				weight = getEdgeWeight(n1, n2, dtImage);
			}
			else if(type == 3)
			{
				weight = getEdgeWeight(n1, n2, m_inputImage); //(doesn't use power parameter)
			}
			else if(type == 4)
			{
				weight = getGeodesicEdgeWeight(n1, n2, m_inputImage);
			}
			edgeWeight.push_back(weight);
		}//end for j
	}//end for i

	if(debug)
		std::cerr << "Finished making edges = " << edgeArray.size() << std::endl;

	return true;
}

float MST::getXiaosongEdgeWeight(fPoint3D n1, fPoint3D n2, double pwr)
{
	if(!m_inputImage)
		return 10000;

	//Now compute edge weight:
	float densityFactor = 0;
	if(pwr >= 0.1)
	{
		//Get the two intensity values (of closest voxel - no interpolation)
		ImageType::IndexType index1;
		index1[0] = roundToInt((double)n1.x);
		index1[1] = roundToInt((double)n1.y);
		index1[2] = roundToInt((double)n1.z);
		ImageType::IndexType index2;
		index2[0] = roundToInt((double)n2.x);
		index2[1] = roundToInt((double)n2.y);
		index2[2] = roundToInt((double)n2.z);
		PixelType pix1 = m_inputImage->GetPixel(index1);
		PixelType pix2 = m_inputImage->GetPixel(index2);

		//Adjust the intensity values
		float i1adj = (float)pow(double(pix1)+0.001, pwr);
		float i2adj = (float)pow(double(pix2)+0.001, pwr);

		//Compute the densityFactor:
		float term1 = (float)pow(double(i1adj+i2adj),double(1.05));
		// AAAHH THIS IS GROSS AND DOES NOT MAKE SENSE
		//float term2 = pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5));
		float term2 = 0;
		densityFactor = fabs(term1 + term2);
	}

	float dx = n1.x - n2.x;
	float dy = n1.y - n2.y;
	float dz = n1.z - n2.z;
	float num = (float)sqrt(float(dx*dx + dy*dy + dz*dz));
	float den = densityFactor*.02 + 1;
	float weight = num / den;
	return weight;
}

//By Xiao Liang (modified by Isaac Abbott):
float MST::getEdgeWeight(fPoint3D n1, fPoint3D n2, ImageType::Pointer img)
{
	if(!img)
		return 10000;

	ImageType::IndexType index1;
	index1[0] = roundToInt((double)n1.x);
	index1[1] = roundToInt((double)n1.y);
	index1[2] = roundToInt((double)n1.z);
	ImageType::IndexType index2;
	index2[0] = roundToInt((double)n2.x);
	index2[1] = roundToInt((double)n2.y);
	index2[2] = roundToInt((double)n2.z);
	
	PixelType dti = img->GetPixel(index1);
	PixelType dtj = img->GetPixel(index2);
			
	float dx = n1.x - n2.x;
	float dy = n1.y - n2.y;
	float dz = n1.z - n2.z;
    //float term1 = (float)pow(double(dti+dtj),double(1.05));
	float term1 = dti + dtj; 
	float densityFactor = fabs(term1);
	float num = (float)sqrt(float(dx*dx + dy*dy + dz*dz));
	float weight = 2*num / densityFactor;

	return weight;
}

float MST::getGeodesicEdgeWeight(fPoint3D n1, fPoint3D n2, ImageType::Pointer img)
{
	if(!img)
		return 10000;

	float dx = n2.x - n1.x;
	float dy = n2.y - n1.y;
	float dz = n2.z - n1.z;

	//Find 5 points along the edge
	std::vector<fPoint3D> testPoints;
	for(double t=0.0; t<=1; t+=1.0/4.0)
	{
		fPoint3D tPoint;
		tPoint.x = (dx)*t + n1.x;
		tPoint.y = (dy)*t + n1.y;
		tPoint.z = (dz)*t + n1.z;
		testPoints.push_back(tPoint);
	}

	double startI = 0;	//intensity of start point
	double avg=0;		//average intensity of other points
	double div = testPoints.size()-1;
	for(int i=0; i<(int)testPoints.size(); ++i)
	{
		fPoint3D tP = testPoints.at(i);
		ImageType::IndexType index1;
		index1[0] = roundToInt((double)tP.x);
		index1[1] = roundToInt((double)tP.y);
		index1[2] = roundToInt((double)tP.z);
		PixelType val = img->GetPixel(index1);
		if(i==0)
			startI = (double)val;
		else
			avg += (double)val/div;
	}

	//Now compute the weight. 
	//Short distance means small weight
	//But big difference in intensities means large weight
	float term1 = 2 * (float)sqrt(float(dx*dx + dy*dy + dz*dz));
	float term2 = (float)sqrt(float((startI-avg)*(startI-avg)));
	float weight = term1 + term2;
	return weight;
}

bool MST::minimumSpanningTree()
{
	int num_edges = (int)edgeArray.size();
	int num_nodes = (int)nodes.size();

	if( num_edges <= 0 )
		return false;

	//Convert std c++ vector to c array:
	pairE * edge_array = new pairE[num_edges];
	float * edge_wght = new float[num_edges];
	if(!edge_array || !edge_wght)
		return false;

	for(int i = 0; i<num_edges; ++i)
	{
		edge_array[i] = edgeArray.at(i);
		edge_wght[i] = edgeWeight.at(i);
	}

	edgeArray.clear();
	edgeWeight.clear();

	//Create graph of all nodes:
	if(nodeGraph) 
		delete nodeGraph;
	nodeGraph = new Graph(edge_array, edge_array + num_edges, edge_wght, num_nodes);
	if(!nodeGraph) 
		return false;

	delete[] edge_array;
	delete[] edge_wght;

	//typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightType;
	//WeightType weight = get(boost::edge_weight, *nodeGraph);
	
	// MST algorithm
	spanningTree.clear();
	boost::kruskal_minimum_spanning_tree(*nodeGraph, back_inserter(spanningTree));

	if(debug)
		std::cerr << "kruskal_minimum_spanning_tree(MST) is finished!" << std::endl;

	// Create a graph for the initial MST
	if(mstGraph)
		delete mstGraph;
	mstGraph = new Graph(1);
	if(!mstGraph)
		return false;

	int num_edge_MST = 0;
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		add_edge(source(*ei, *nodeGraph), target(*ei, *nodeGraph), *mstGraph);
		num_edge_MST++;
	}

	if(debug)
		std::cerr << "MST edges = " << num_edge_MST << std::endl;

	if(debug)
	{
		std::vector<pairE> mstLines;

		//Get the lines from the mst:
		typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
		Edge_iter ei, ei_end; 
		for (tie(ei, ei_end) = edges(*mstGraph); ei != ei_end; ++ei)
		{
			pairE ne( (int)source(*ei, *mstGraph)-1, (int)target(*ei, *mstGraph)-1 );
			mstLines.push_back( ne );
		}

		std::cerr << "Number of mst lines = " << mstLines.size() << std::endl;

		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes);
		fhdl->SetLines(&mstLines);
		fhdl->Write("InitialMST.vtk");
		delete fhdl;
	}

	return true;
}

bool MST::ErodeAndDialateNodeDegree(int mophStrength)
{
	if((int)spanningTree.size() == 0 || !nodeGraph)
	{
		std::cerr << "missing tree or g\n";
		return false;
	}

	int num_nodes = (int)nodes.size();

	int * degree_nodes = new int[num_nodes+1];
	int * degree_nodes_buffer = new int[num_nodes+1];
	if(!degree_nodes || !degree_nodes_buffer)
	{
		std::cerr << "Could not allocate degree buffers\n";
		return false;
	}

	//Inititialize for all the vertices.
	//Actual vertex index starts from 1
	for(int i=0; i<num_nodes+1; ++i)
	{
		degree_nodes[i] = 0;
		degree_nodes_buffer[i] = 0;
	}

	// create initial degree_nodes array
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		degree_nodes[source(*ei, *nodeGraph)] ++;
		degree_nodes[target(*ei, *nodeGraph)] ++;
	}

	// -- Erosion and Dilation of MST
	int times_erosion = mophStrength;

	int * edge_eroded = new int[num_nodes*2];
	if(!edge_eroded)
	{
		std::cerr << "Could not allocate edge_eroded\n";
		return false;
	}

	int num_edge_eroded = 0;
	while (times_erosion != 0) 
	{
		times_erosion--;
		for (int i=1; i<=num_nodes; i++)   
			degree_nodes_buffer[i] = degree_nodes[i];

		for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
		{
			if (degree_nodes_buffer[source(*ei, *nodeGraph)]>0 && degree_nodes_buffer[target(*ei, *nodeGraph)]>0)  
			{
				if (degree_nodes_buffer[source(*ei, *nodeGraph)]==1 || degree_nodes_buffer[target(*ei, *nodeGraph)]==1)  
				{
					degree_nodes[source(*ei, *nodeGraph)] --;
					degree_nodes[target(*ei, *nodeGraph)] --;
          
					// Save the edges eroded in a stack-like array. Each edge takes two elements of the array
					edge_eroded[num_edge_eroded*2]  = (int)source(*ei, *nodeGraph); // Saving of eroded edges
					edge_eroded[num_edge_eroded*2+1] = (int)target(*ei, *nodeGraph);
					num_edge_eroded ++;
				}//end if	
			}// end if
		}// end for (vector...
	}// end while (times_erosion != 0)

	if(debug)
		std::cerr << "Erosion of MST is finished!" << std::endl;

	// Dilation the MST by counting up the degree of nodes
	while (num_edge_eroded !=0) 
	{
		num_edge_eroded --;
		int edge_source = edge_eroded[num_edge_eroded*2];  // Read the stored eroded edges
		int edge_target = edge_eroded[num_edge_eroded*2+1];
		if ((degree_nodes[edge_source]+ degree_nodes[edge_target]) == 1)  
		{  // if a branch tip edge
			degree_nodes[edge_source] ++;
			degree_nodes[edge_target] ++;
		}
	}

	if(debug)
		std::cerr << "Dilation of MST is finished!" << std::endl;

	//copy int output vectors:
	nodeDegree.clear();
	for(int i=0; i<num_nodes; ++i)
	{
		nodeDegree.push_back( degree_nodes[i+1] );
	}

	delete[] edge_eroded;
	delete[] degree_nodes_buffer;
	delete[] degree_nodes;
		
	if(debug)
		std::cerr << " DIALATED NODES = " << nodeDegree.size() << std::endl;

	if(debug)
	{
		FILE * fout = fopen("degrees.txt", "w");
		if (fout != NULL)
		{
			for(int i=0; i<(int)nodeDegree.size(); ++i)
			{
				fprintf(fout,"%d %f\n", i, (float)nodeDegree.at(i));
			}
			fclose(fout); 
		}
	}
	return true;
}

std::vector<pairE> MST::BackboneExtract()
{
	std::vector<pairE> retLines;

	if(!mstGraph || (int)nodeDegree.size() == 0)
	{
		std::cerr << "missing tree or mstGraph or node degrees\n";
		return retLines;
	}

	//Get the backbone lines from the mst
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(*mstGraph); ei != ei_end; ++ei)
	{
		int n1 = (int)source(*ei, *mstGraph)-1; //node 1
		int n2 = (int)target(*ei, *mstGraph)-1;  //node 2
		if(nodeDegree.at(n1)!=0 && nodeDegree.at(n2)!=0)
		{
			pairE ne( n1, n2 );
			retLines.push_back( ne );
		}
   }

	if(debug)
		std::cerr << "Number of backbone lines = " << retLines.size() << std::endl;

	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes);
		fhdl->SetLines(&retLines);
		fhdl->Write("BackboneCandidate.vtk");
		delete fhdl;
	}

	return retLines;
}

//SearchFirst
std::vector<pairE> MST::SearchFirstandSecondLevelBranch(void)
{  
	std::vector<pairE> retLines;
    int num_nodes = (int)nodes.size();

	if(num_nodes == 0 || (int)nodeDegree.size() == 0 || !mstGraph)
	{return retLines;
	}

	// Create a Backbone vertice flag array
	bool * vertBackbone = new bool[num_nodes+1];
	for (int i=0; i<num_nodes; i++)   
	{
		if (nodeDegree.at(i) >= 1)  
			vertBackbone[i] = true;
		else
			vertBackbone[i] = false;
	}

    Graph prunedGraph (num_nodes+1);
	prunedGraph = morphGraphPrune(mstGraph, &nodes, PruneThreshold);
    
    typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	typedef boost::graph_traits < Graph >::vertex_iterator Vertex_iter;
	typedef boost::graph_traits<Graph>::out_edge_iterator Outedge_iter;
    Outedge_iter  outei, outedge_end;
	Outedge_iter  outei2, outedge_end2;
    Outedge_iter  outei3, outedge_end3;
	Vertex_iter vi, vend;

	// Spine Candidate graph created, is to save the possible spine
	Graph msTreeSpineCandidate(num_nodes+1);    

    const int MAXNumBranch = 10;
	//int index_vert;
    int vertsCurBranch2[MAXNumBranch][2000];  // For the 2nd level branch at BB (at most MAXNumBranch 2nd level branches, at most 2000 vertices)
    
	int vertsCurBr_Index2[MAXNumBranch];      // when array at [0], it is the 1st level branch at BB
    for (int i=0; i<MAXNumBranch; i++)  vertsCurBr_Index2[i]=0;  // Initialize to zeros

	typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, prunedGraph);  // get index map of vertices

  // for all vertex in the graph
    for(boost::tie(vi, vend) = vertices(prunedGraph); vi != vend; ++vi)
     { 
       vertsCurBr_Index2[0] = 0;
	   // Get index from the graph IndexMap
       int index_vert = int(index[*vi]); 
       // if the vertex is not on Backbone, continue
       if (vertBackbone[index_vert]==0) 
		   continue;  
       int outdegree = (int)out_degree(*vi, prunedGraph);
       int numBranch_on_Backbone = outdegree - nodeDegree.at(index_vert);
       // if it has at least one branch out of BackBone  
       if (numBranch_on_Backbone <= 0)  continue;  

       // For each out branch (edge) on the BackBone
       for (boost::tie(outei, outedge_end) = out_edges(*vi, prunedGraph); outei != outedge_end; outei++)
        {
        int targ = (int)target(*outei, prunedGraph);
	    // the edge is on BackBone, continue
        if (vertBackbone[targ])  continue;  
      
	    vertsCurBranch2[0][vertsCurBr_Index2[0]] = (int)source(*outei, prunedGraph);
        vertsCurBr_Index2[0]++;
        vertsCurBranch2[0][vertsCurBr_Index2[0]] = (int)target(*outei, prunedGraph);

        while (out_degree(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph), prunedGraph) == 2)
         {
          for (boost::tie(outei2, outedge_end2) =
              out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph),prunedGraph);
              outei2 != outedge_end2; ++outei2)
            {
                if (target(*outei2, prunedGraph) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
                 continue;
                vertsCurBranch2[0][vertsCurBr_Index2[0]+1] = (int)target(*outei2, prunedGraph);
            } // end for
		  vertsCurBr_Index2[0]++; 
	     } // end while 

        // find the first level branch and add them to msTreeSpineCandidate
        for (int j = 1; j <= vertsCurBr_Index2[0]; j++) 
        {
         add_edge(vertsCurBranch2[0][j-1], vertsCurBranch2[0][j], msTreeSpineCandidate);  
        } // end for

       // Begin to check the 2nd level of branches located at BB
       int ind2Brch = 0;
        // For each 2nd level branch starting from the end of 1st level branch
       for (boost::tie(outei2, outedge_end2) =out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph),prunedGraph);
           outei2 != outedge_end2; ++outei2)
         {
          if (target(*outei2, prunedGraph) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
             continue;  // continue if the out edge belongs to the old branch
          ind2Brch++;
          vertsCurBranch2[ind2Brch][0] = vertsCurBranch2[0][vertsCurBr_Index2[0]];
          vertsCurBr_Index2[ind2Brch] = 1;
          vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]] = (int)target(*outei2, prunedGraph);
          // Search for the end of 2nd level branch
          while (out_degree(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], prunedGraph), prunedGraph) == 2)
          {
          //curVert = vertex(vertsCurBranch[vertsCurBr_Index], msTree);
            for (boost::tie(outei3, outedge_end3) = out_edges(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], prunedGraph), prunedGraph);
                                                outei3 != outedge_end3; ++outei3) 
		     {
               if (target(*outei3, prunedGraph) == (unsigned int)vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]-1])
                continue;
               vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]+1] = (int)target(*outei3, prunedGraph);
             } // end for
              vertsCurBr_Index2[ind2Brch]++;
           } // end while 
     
        for (int j = 1; j <= vertsCurBr_Index2[ind2Brch]; j++) 
         {
          add_edge(vertsCurBranch2[ind2Brch][j-1], vertsCurBranch2[ind2Brch][j], msTreeSpineCandidate);   // add branch for the 2nd level
         }
      } // End of 2nd level branch
	} // // End of each out edge
  } // End of all vertice
    
    //msTreeSpineCandidate = morphGraphPrune(msTreeSpineCandidate, &nodes, PruneThreshold);

    //Get the spine lines from the Graph DetectedSpine
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(msTreeSpineCandidate); ei != ei_end; ++ei)
	{
		int n1 = (int)source(*ei, msTreeSpineCandidate)-1; //node 1
		int n2 = (int)target(*ei, msTreeSpineCandidate)-1;  //node 2
		pairE ne( n1, n2 );
	    retLines.push_back( ne );
    }

	if(debug)
		std::cerr << "Number of  lines = " << retLines.size() << std::endl;

	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes);
		fhdl->SetLines(&retLines);
		fhdl->Write("FirstandSecondLevelBranch.vtk");
		delete fhdl;
	}

    delete []vertBackbone;
	return retLines;
}

std::vector<pairE> MST::SpineExtract()
{
	std::vector<pairE> retLines;

	const int MAXNumBranch = 10;
    double mahalanobis_dist[MAXNumBranch];
    double mahalanobis_dist_nonSpine[MAXNumBranch];
	double meanDensityBranch[MAXNumBranch];   // Suppose at most MAXNumBranch branches at the 2nd level branch from BB
    double meanVesselBranch[MAXNumBranch];
    double length_leaf[MAXNumBranch];
    double aveDensityBranch[MAXNumBranch];   // Suppose at most MAXNumBranch branches at the 2nd level branch from BB
    double aveVesselBranch[MAXNumBranch];
    double length_2leaf[MAXNumBranch];  // length of two level branches
	int indVert,indVert_last;
    int slsz = sizeX*sizeY;   // slice size
    int sz = slsz*sizeZ;
	
	int num_nodes = (int)nodes.size();

	if(num_nodes == 0 || (int)nodeDegree.size() == 0 || !mstGraph)
	{return retLines;
	}

	// Create a Backbone vertice flag array
	bool * vertBackbone = new bool[num_nodes+1];
	for (int i=0; i<num_nodes; i++)   
	{
		if (nodeDegree.at(i) >= 1)  
			vertBackbone[i] = true;
		else
			vertBackbone[i] = false;
	}

	
	Graph prunedGraph (num_nodes+1);
	prunedGraph = morphGraphPrune(mstGraph, &nodes, PruneThreshold);
    //Graph prunedGraph = *mstGraph;
    typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	typedef boost::graph_traits < Graph >::vertex_iterator Vertex_iter;
	typedef boost::graph_traits<Graph>::out_edge_iterator Outedge_iter;
    Outedge_iter  outei, outedge_end;
	Outedge_iter  outei2, outedge_end2;
    Outedge_iter  outei3, outedge_end3;
	Vertex_iter vi, vend;
	// Spine Candidate graph created, is to save the possible spine
	Graph msTreeSpineCandidate(num_nodes+1);     
    Graph DetectedSpine(num_nodes+1); 
	int NumberNodesofSpineCandidate=0;
    int NumberNodesofRealSpine=0;
    
    int vertsCurBranch2[MAXNumBranch][2000];  // For the 2nd level branch at BB (at most MAXNumBranch 2nd level branches, at most 2000 vertices)
    int vertsCurBr_Index2[MAXNumBranch];      // when array at [0], it is the 1st level branch at BB

    for (int i=0; i<MAXNumBranch; i++)  vertsCurBr_Index2[i]=0;  // Initialize to zeros
	
    double sample[3];
	int num_leaves = 0;
	typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, prunedGraph);  // get index map of vertices
	MDLClassifier LDA_RealSpine(3);
    MDLClassifier LDA_NonSpine(3);
	int LDA_t1= LDA_RealSpine.MeanVectorandVarianceMatrix((char *)RealSpineFeatureFilename);
	int LDA_t2= LDA_NonSpine.MeanVectorandVarianceMatrix((char *)NonSpineFeatureFilename);
	if (debug)
	{
      if (LDA_t1 >0)
      {
	   std::cout << " There is Real-Spine Feature Sample,We do machine learning based classification" << std::endl;
      }
      else 
      {
	  std::cout << " There is not Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
      }
      if (LDA_t2 >0)
      {
	  std::cout << " There is Non-Spine Feature Sample,We do machine learning based classification" << std::endl;
      }
      else 
      {
	   std::cout << " There is Non-Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
      }
	}
 

  // for all vertex in the graph
   for(boost::tie(vi, vend) = vertices(prunedGraph); vi != vend; ++vi)
    { 
    vertsCurBr_Index2[0] = 0;

	// Get index from the graph IndexMap
    int index_vert = int(index[*vi]); 

    // if the vertex is not on Backbone, continue
    if (vertBackbone[index_vert]==0) 
		continue;  

    int outdegree = (int)out_degree(*vi, prunedGraph);

    int numBranch_on_Backbone = outdegree - nodeDegree.at(index_vert);

    // if it has at least one branch out of BackBone  
    if (numBranch_on_Backbone <= 0)  continue;  
    // For each out branch (edge) on the BackBone

    for (boost::tie(outei, outedge_end) = out_edges(*vi, prunedGraph); outei != outedge_end; outei++)
      {
      int targ = (int)target(*outei, prunedGraph);

	  // the edge is on BackBone, continue
      if (vertBackbone[targ])  continue;  
      meanDensityBranch[0] = 0;
      meanVesselBranch[0] = 0;
      vertsCurBranch2[0][vertsCurBr_Index2[0]] = (int)source(*outei, prunedGraph);
      vertsCurBr_Index2[0]++;
      vertsCurBranch2[0][vertsCurBr_Index2[0]] = (int)target(*outei, prunedGraph);
      num_leaves++;

      while (out_degree(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph), prunedGraph) == 2)
       {
        for (boost::tie(outei2, outedge_end2) =
             out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph),prunedGraph);
             outei2 != outedge_end2; ++outei2)
          {
            if (target(*outei2, prunedGraph) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
              continue;
            vertsCurBranch2[0][vertsCurBr_Index2[0]+1] = (int)target(*outei2, prunedGraph);
          } // end for
        
		vertsCurBr_Index2[0]++;
        
	   } // end while 

      // Evaluate with MDL if the branch is chosen
      bool branchChosen = 1;
      length_leaf[0] = 0;

      for (int j = 0; j <= vertsCurBr_Index2[0]; j++)
        {
         indVert = vertsCurBranch2[0][j];
     
	     // output the locations 
         //fprintf(fclass_identify, "%d  %6.2f %6.2f %6.2f\n", num_leaves, vertexPos[indVert].x, vertexPos[indVert].y, vertexPos[indVert].z);
        
          if (j==0)indVert_last = indVert;
          if(indVert >= (int)this->nodes.size())indVert = indVert_last;

		  fPoint3D vertexPosStart,vertexPosEnd;
		  vertexPosStart.x = this->nodes.at(indVert).x;
		  vertexPosStart.y = this->nodes.at(indVert).y;
          vertexPosStart.z = this->nodes.at(indVert).z;
		  vertexPosEnd.x = this->nodes.at(indVert_last).x;
          vertexPosEnd.y = this->nodes.at(indVert_last).y;
		  vertexPosEnd.z = this->nodes.at(indVert_last).z;

          double length_edge = (vertexPosStart.x-vertexPosEnd.x)*(vertexPosStart.x-vertexPosEnd.x);
		  length_edge+= (vertexPosStart.y-vertexPosEnd.y)*(vertexPosStart.y-vertexPosEnd.y);
		  length_edge+= (vertexPosStart.z-vertexPosEnd.z)*(vertexPosStart.z-vertexPosEnd.z);

		  length_edge = sqrt(length_edge);
		  length_leaf[0] += length_edge;
		  indVert_last = indVert;
          
          ImageType::IndexType index;
		  index[0] = roundToInt((double) vertexPosStart.x);
		  index[1] = roundToInt((double) vertexPosStart.y);
		  index[2] = roundToInt((double) vertexPosStart.z);
	      PixelType pix1 = m_inputImage->GetPixel(index);
		  PixelType pix2 = m_VesselMap->GetPixel(index);
		  meanDensityBranch[0] += pix1;
          meanVesselBranch[0] += pix2;
        } // end for 
      
	    meanDensityBranch[0] = meanDensityBranch[0]/ (vertsCurBr_Index2[0]+1);
        meanVesselBranch[0] = meanVesselBranch[0] / (vertsCurBr_Index2[0]+1);

        sample[0] =  meanDensityBranch[0];
        sample[1] =  length_leaf[0];
        sample[2] =  meanVesselBranch[0];

    
	   if (LDA_t1>0)
          mahalanobis_dist[0] = LDA_RealSpine.MahalanobisDist(sample); 
       else 
          mahalanobis_dist[0] = LDA_RealSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 1);
       if (LDA_t2>0)
          mahalanobis_dist_nonSpine[0] = LDA_NonSpine.MahalanobisDist(sample);
       else 
          mahalanobis_dist_nonSpine[0] = LDA_NonSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 0);

       int ind2Brch = 0;

        // For each 2nd level branch starting from the end of 1st level branch
       for (boost::tie(outei2, outedge_end2) =out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], prunedGraph),prunedGraph);
           outei2 != outedge_end2; ++outei2)
         {
          if (target(*outei2, prunedGraph) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
             continue;  // continue if the out edge belongs to the old branch
          ind2Brch++;
          vertsCurBranch2[ind2Brch][0] = vertsCurBranch2[0][vertsCurBr_Index2[0]];
          vertsCurBr_Index2[ind2Brch] = 1;
          vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]] = (int)target(*outei2, prunedGraph);
          // Search for the end of 2nd level branch
          while (out_degree(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], prunedGraph), prunedGraph) == 2)
          {
          //curVert = vertex(vertsCurBranch[vertsCurBr_Index], msTree);
            for (boost::tie(outei3, outedge_end3) = out_edges(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], prunedGraph), prunedGraph);
                                                outei3 != outedge_end3; ++outei3) 
		    {
              if (target(*outei3, prunedGraph) == (unsigned int)vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]-1])
                continue;
               vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]+1] = (int)target(*outei3, prunedGraph);
            } // end for
              vertsCurBr_Index2[ind2Brch]++;
           } // end while 

        // Compute the feature-based description length for the 2nd level branch
          length_leaf[ind2Brch] = 0;
          meanDensityBranch[ind2Brch] = 0;
          meanVesselBranch[ind2Brch] = 0;
          num_leaves++;  // for training purpose
          
		  for (int j = 0; j <= vertsCurBr_Index2[ind2Brch]; j++)
          {
            indVert = vertsCurBranch2[ind2Brch][j];
          
		    // second level feature 
            //fprintf(fclass_identify, "%d  %6.2f %6.2f %6.2f\n", num_leaves, vertexPos[indVert].x, vertexPos[indVert].y, vertexPos[indVert].z);
         
            if (j==0) indVert_last = indVert;
			if(indVert >= (int)this->nodes.size())indVert = indVert_last;
				
            fPoint3D vertexPosStart,vertexPosEnd;
			vertexPosStart.x = this->nodes.at(indVert).x;
			vertexPosStart.y = this->nodes.at(indVert).y;
            vertexPosStart.z = this->nodes.at(indVert).z;
			vertexPosEnd.x = this->nodes.at(indVert_last).x;
            vertexPosEnd.y = this->nodes.at(indVert_last).y;
			vertexPosEnd.z = this->nodes.at(indVert_last).z;

            double length_edge = (vertexPosStart.x-vertexPosEnd.x)*(vertexPosStart.x-vertexPosEnd.x);
		    length_edge+= (vertexPosStart.y-vertexPosEnd.y)*(vertexPosStart.y-vertexPosEnd.y);
			length_edge+= (vertexPosStart.z-vertexPosEnd.z)*(vertexPosStart.z-vertexPosEnd.z);

			length_edge = sqrt(length_edge);
			length_leaf[ind2Brch] += length_edge;
			indVert_last = indVert;
			ImageType::IndexType index;

		    index[0] = roundToInt((double) vertexPosStart.x);
		    index[1] = roundToInt((double) vertexPosStart.y);
		    index[2] = roundToInt((double) vertexPosStart.z);
	        PixelType pix1 = m_inputImage->GetPixel(index);
		    PixelType pix2 = m_VesselMap->GetPixel(index);
		    meanDensityBranch[ind2Brch] += pix1;
            meanVesselBranch[ind2Brch] += pix2;   
          } // end for

         meanDensityBranch[ind2Brch] = meanDensityBranch[ind2Brch]/ (vertsCurBr_Index2[ind2Brch]+1);
         meanVesselBranch[ind2Brch] = meanVesselBranch[ind2Brch] / (vertsCurBr_Index2[ind2Brch]+1);

         // Compute the average features of two level branches
         length_2leaf[ind2Brch] = length_leaf[ind2Brch] + length_leaf[0];
         aveDensityBranch[ind2Brch] = (meanDensityBranch[ind2Brch]*length_leaf[ind2Brch] +  meanDensityBranch[0]*length_leaf[0]);
         aveDensityBranch[ind2Brch] = aveDensityBranch[ind2Brch] / length_2leaf[ind2Brch];
         aveVesselBranch[ind2Brch] = (meanVesselBranch[ind2Brch]*length_leaf[ind2Brch] + meanVesselBranch[0]*length_leaf[0]);
         aveVesselBranch[ind2Brch] = aveVesselBranch[ind2Brch] / length_2leaf[ind2Brch];
	     sample[0] =  aveDensityBranch[ind2Brch];
         sample[1] =  length_2leaf[ind2Brch];
         sample[2] =  aveVesselBranch[ind2Brch];
        
	
		 if (LDA_t1>0)
          mahalanobis_dist[ind2Brch] = LDA_RealSpine.MahalanobisDist(sample); 
         else 
          mahalanobis_dist[ind2Brch] = LDA_RealSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 1);
         if (LDA_t2>0)
          mahalanobis_dist_nonSpine[ind2Brch] = LDA_NonSpine.MahalanobisDist(sample);
         else 
          mahalanobis_dist_nonSpine[ind2Brch] = LDA_NonSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 0);

      } // End of 2nd level branch

    
     //----------------------------MDL fitness for spine --------------------------------------//
 
      int MDL_minIndex = -1; // Indicate that empty model set is chosen
      double sum_mahalanobis_nonSpine = 0;
      for (int i = 0; i<= ind2Brch; i++)
        sum_mahalanobis_nonSpine += mahalanobis_dist_nonSpine[i];
      
	  double MDL_min = sum_mahalanobis_nonSpine;

      // 2. Only 1st level branch model set
      double MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[0] + mahalanobis_dist[0];
      MDL = MDL + (1-Alpha)*(1/Alpha)*(-14.0/3.0);      // alpha represents the model description length of one branch
      if (MDL < MDL_min)
        {
        MDL_min = MDL;
        MDL_minIndex = 0;  
        }
      // 3. Two level branch model set, including 1st level and 2nd level branches
      for (int i = 1; i <= ind2Brch; i++)
        {
        MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[i] + mahalanobis_dist[i];
        MDL = MDL + (1-Alpha)*(1/Alpha)*(-14.0/3.0);
        if (MDL <  MDL_min)
          {
          MDL_min = MDL;
          MDL_minIndex =  i;
          }
        } //end for
    
      // Adding Spine model - Only add the branches with the minimal MDL
       if (MDL_minIndex >= 0)
        {
         for (int j = 1; j <= vertsCurBr_Index2[0]; j++)
           {
            add_edge(vertsCurBranch2[0][j-1], vertsCurBranch2[0][j], DetectedSpine); 
           }
         if (MDL_minIndex >= 1)
           {
           for (int j = 1; j <= vertsCurBr_Index2[MDL_minIndex]; j++)
            {
             add_edge(vertsCurBranch2[MDL_minIndex][j-1],
                     vertsCurBranch2[MDL_minIndex][j],  DetectedSpine);
            } // end for
          }// end if
        }// end if(MDL_minIndex >= 0)

      // 1. Empty model set is chosen (no branch)
       sum_mahalanobis_nonSpine = 0;
       for (int i = 0; i<= ind2Brch; i++) 
       {
        sum_mahalanobis_nonSpine += mahalanobis_dist_nonSpine[i];
       }

      // 2. Only 1st level branch model set
       MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[0] + mahalanobis_dist[0];
       MDL_minIndex = 0; // Indicate that first level branch model set is chosen
       MDL_min = MDL;

      // 3. Two level branch model set, including 1st level and 2nd level branches
       for (int i = 1; i <= ind2Brch; i++) 
       {
         MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[i] + mahalanobis_dist[i];
         if (MDL <  MDL_min) 
         {
          MDL_min = MDL;
          MDL_minIndex =  i;
         }
       }// end for 
     
	   //fprintf(fout_MDL, "%10d %20f %20f\n", ind2Brch, sum_mahalanobis_nonSpine, MDL_min);     
   }    // End of each out edge 
  } // End of all vertice


    //Get the spine lines from the Graph DetectedSpine
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(DetectedSpine); ei != ei_end; ++ei)
	{
		int n1 = (int)source(*ei, DetectedSpine)-1; //node 1
		int n2 = (int)target(*ei, DetectedSpine)-1;  //node 2
		pairE ne( n1, n2 );
		retLines.push_back( ne );
    }

	if(debug)
		std::cerr << "Number of spine lines = " << retLines.size() << std::endl;
	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes);
		fhdl->SetLines(&retLines);
		fhdl->Write("MDLSpine.vtk");
		delete fhdl;
	}
	return retLines;
}

// Pruning on graph for trivia branches
// Removing branches with length below certain threshold
Graph MST::morphGraphPrune(Graph *graph, std::vector<fPoint3D> *nodes, float lengthThreshold)
{ 
    typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	typedef boost::graph_traits < Graph >::vertex_iterator Vertex_iter;
	typedef boost::graph_traits<Graph>::out_edge_iterator Outedge_iter;

	Edge_iter   ei, ei_end;
	Outedge_iter  outei, outedge_end;
	Vertex_iter vi, vend;

	//I probably want to create a copy of the graph somehow??
	//Graph ng = *graph;
    Graph ng = *graph;
    /*
    for (tie(ei, ei_end) = edges(*graph); ei != ei_end; ++ei) 
	{
		add_edge(source(*ei, *graph), target(*ei, *graph), ng);
	}
    */
	int curBranchVerts[2000];
	int curBrVerts_Index = 0;
	float length_edge = 0; 
	int indVert, indVert_last;
    int branchChosen;
	float length_leaf = 0;

	// Consider all leaves
	int num_leaves = 0;

    // for all vertex in the graph
	for(boost::tie(vi, vend) = vertices(ng); vi != vend; ++vi) 
	{ 
		curBrVerts_Index = 0;
		if (out_degree(*vi, ng) == 1)  
		{ // if it is a leaf tip
			for (boost::tie(outei, outedge_end) = out_edges(*vi, ng); outei != outedge_end; ++outei)
			{
				curBranchVerts[curBrVerts_Index] = (int)source(*outei, ng);
				curBrVerts_Index++;
				curBranchVerts[curBrVerts_Index] = (int)target(*outei, ng);
			} // end for

			while (out_degree(vertex(curBranchVerts[curBrVerts_Index], ng), ng) == 2) 
			{
				
				for (boost::tie(outei, outedge_end) = out_edges(vertex(curBranchVerts[curBrVerts_Index], ng), ng) ; outei != outedge_end; ++outei)
				{
					if (target(*outei, ng) == (unsigned int)curBranchVerts[curBrVerts_Index-1])
							continue;
					curBranchVerts[curBrVerts_Index+1] = (int)target(*outei, ng);
				} //end for
				curBrVerts_Index++;
			} // end while
				
			branchChosen = 1;// Evaluate with MDL if the branch is chosen
				
			length_leaf = 0;  
			for (int j = 0; j <= curBrVerts_Index; j++)
			{
				indVert = curBranchVerts[j];
				if(indVert >= (int)nodes->size())
				{
					std::cerr << "not a node\n";
					indVert = indVert_last;
				}

				if (j==0) 
				{
					indVert_last = indVert;
				} // end if
                    
				fPoint3D vertexPosStart,vertexPosEnd;
				vertexPosStart.x = nodes->at(indVert).x;
				vertexPosStart.y = nodes->at(indVert).y;
                vertexPosStart.z = nodes->at(indVert).z;
				vertexPosEnd.x = nodes->at(indVert_last).x;
                vertexPosEnd.y = nodes->at(indVert_last).y;
				vertexPosEnd.z = nodes->at(indVert_last).z;

                length_edge = (vertexPosStart.x-vertexPosEnd.x)*(vertexPosStart.x-vertexPosEnd.x);
				length_edge+= (vertexPosStart.y-vertexPosEnd.y)*(vertexPosStart.y-vertexPosEnd.y);
				length_edge+= (vertexPosStart.z-vertexPosEnd.z)*(vertexPosStart.z-vertexPosEnd.z);

				length_edge = sqrt(length_edge);
				length_leaf += length_edge;
				indVert_last = indVert;
				//if(debug)
					//std::cout << "length_leaf = " << length_leaf << std::endl;
			} //end for

			//if(debug)
			//	std::cout << "I am here" << std::endl;

			if (length_leaf < lengthThreshold)    
			{  
				branchChosen = 0;     
			}
			// Remove the branch based on leaf-length critierion
			if (branchChosen == 0)  
			{
				for (int j = 1; j <= curBrVerts_Index; j++) 
				{
					remove_edge(curBranchVerts[j-1], curBranchVerts[j], ng);
				}
			}
		}
	}

	if (debug)
		std::cerr << "The Graph Pruning is done with respect to the threshold = " << lengthThreshold << std::endl;

    if (debug)
	{
		std::vector<pairE> retLines;
		Edge_iter ei, ei_end; 
		for (tie(ei, ei_end) = edges(ng); ei != ei_end; ++ei)
		{
			pairE ne( (int)source(*ei, ng)-1, (int)target(*ei, ng)-1 );
			retLines.push_back( ne );
		}

		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(nodes);
		fhdl->SetLines(&retLines);
		fhdl->Write("PrunedGraph.vtk");
		delete fhdl;
	}

	return ng;
}


}