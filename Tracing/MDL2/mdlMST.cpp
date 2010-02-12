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

bool MST::CreateGraphAndMST()
{
	if(!m_inputImage || !skeletonPoints)
		return false;

	this->skeletonPointsToNodes(useVoxelRounding);
	this->nodesToEdges();
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
bool MST::nodesToEdges()
{
	if( (int)nodes.size() == 0)
		return false;

	int num_nodes = (int)nodes.size();

	//Iterate through the nodes and create edges within the specified range.
	for(int i=0; i<num_nodes; ++i)
	{
		//for(int j=1; j<i; ++j)
		for(int j=0; j<i; ++j)
		{
			fPoint3D n1 = nodes.at(i);
			fPoint3D n2 = nodes.at(j);
			float dx = n1.x - n2.x;
			float dy = n1.y - n2.y;
			float dz = n1.z - n2.z;

			if( abs(dx) > (float)edgeRange )
				continue;
			if( abs(dy) > (float)edgeRange )
				continue;
			if( abs(dz) > (float)edgeRange )
				continue;

			//If I'm here, then I've found a close enough node
			edgeArray.push_back( pairE(i+1,j+1) ); //add an edge (count starts at 1 for nodes)

			//Now compute edge weight:
			float densityFactor = 0;
			if(power >= 0.1)
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
				float i1adj = (float)pow(double(pix1)+0.001, power);
				float i2adj = (float)pow(double(pix2)+0.001, power);

				//Compute the densityFactor:
				float term1 = (float)pow(double(i1adj+i2adj),double(1.05));
				// AAAHH THIS IS GROSS AND DOES NOT MAKE SENSE
				//float term2 = pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5));
				float term2 = 0;
				densityFactor = fabs(term1 + term2);
			}

			float num = (float)sqrt(float(dx*dx + dy*dy + dz*dz));
			float den = densityFactor*.02 + 1;
			float weight = num / den;
			edgeWeight.push_back(weight);
		}//end for j
	}//end for i

	if(debug)
		std::cerr << "Finished making edges = " << edgeArray.size() << std::endl;

	return true;
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

//MDL based Spine Extraction
std::vector<pairE> MST::SpineExtract()
{
	std::vector<pairE> retLines;

	int num_nodes = (int)nodes.size();
	if(num_nodes == 0 || (int)nodeDegree.size() == 0 || !mstGraph)
		return retLines;

	// Create a Backbone vertice flag array
	bool * vertBackbone = new bool[num_nodes+1];
	for (int i=0; i<num_nodes; i++)   
	{
		if (nodeDegree.at(i) >= 1)  
			vertBackbone[i] = true;
		else
			vertBackbone[i] = false;
	}

	/////
	//NOT FINISHED YET::

	Graph prunedGraph = morphGraphPrune(mstGraph, &nodes, 50.0);

	//Xiao Liang: you can work on finishing this function by
	//re-writing the code in MDABasedSplineExtraction.
	//
	//Your output is a vector of line-pairs corresponding to node indexes.
	//

	return retLines;
}

// Pruning on graph for trivia branches
// Removing branches with length below certain threshold
Graph MST::morphGraphPrune(Graph *graph, std::vector<fPoint3D> *nodes, float lengthThreshold)
{ 
    typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	typedef boost::graph_traits < Graph >::vertex_iterator Vertex_iter;
	typedef boost::graph_traits<Graph>::out_edge_iterator Outedge_iter;

	//I probably want to create a copy of the graph somehow??
	Graph ng = *graph;

	Edge_iter   ei, ei_end;
	Outedge_iter  outei, outedge_end;

	int curBranchVerts[2000];
	int curBrVerts_Index = 0;

	float length_edge = 0; 
	int indVert, indVert_last;
    int branchChosen;
	float length_leaf = 0;

	// Consider all leaves
	int num_leaves = 0;

    // for all vertex in the graph
	Vertex_iter vi, vend;
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
				boost::tie(outei, outedge_end) = out_edges(vertex(curBranchVerts[curBrVerts_Index], ng), ng);
				for ( ; outei != outedge_end; ++outei)
				{
					if (target(*outei, ng) == (unsigned int)curBranchVerts[curBrVerts_Index-1])
							continue;
					curBranchVerts[curBrVerts_Index+1] = (int)target(*outei, ng);
				} //end for
				curBrVerts_Index++;
			} // end while
				
			branchChosen = 1;// Evaluate with MDL if the branch is chosen
				
			length_leaf = 0;  // must reset the value, by xiaoliang
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
				vertexPosStart.x = this->nodes.at(indVert).x;
				vertexPosStart.y = this->nodes.at(indVert).y;
                vertexPosStart.z = this->nodes.at(indVert).z;
				vertexPosEnd.x = this->nodes.at(indVert_last).x;
                vertexPosEnd.y = this->nodes.at(indVert_last).y;
				vertexPosEnd.z = this->nodes.at(indVert_last).z;

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