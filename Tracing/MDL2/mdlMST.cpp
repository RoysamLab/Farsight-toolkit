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
	edgeRange = 10;		//Some default values:
	power = 1;

	//input
	skeletonPoints = NULL;

	edge_array = NULL;
	edge_wght = NULL;
	num_edges = 0;
	g = NULL;

	//output
	nodes.clear();
	spanningTree.clear();
	nodeDegree.clear();
}

MST::~MST()
{
	m_inputImage=NULL;
	skeletonPoints=NULL;

	if(edge_array)
	{
		delete[] edge_array;
		edge_array = NULL;
	}
	if(edge_wght)
	{
		delete[] edge_wght;
		edge_wght = NULL;
	}
	if(g)
	{
		delete g;
		g = NULL;
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

	this->skeletonPointsToNodes();
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
bool MST::skeletonPointsToNodes()
{
	if(!skeletonPoints)
		return false;

	bool * nodeAdded = new bool[numPix];
	if(!nodeAdded)	//Couldn't allocate memory
		return false;

	//- Initialize to zero
	for(long idx=0; idx<numPix; idx++)   
	{ 
		nodeAdded[idx] = false;
	}

	nodes.clear();

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

	int num_nodes = (int)nodes.size();

	if(debug)
		std::cerr << "Number of Nodes = " << num_nodes << std::endl;

	return true;
}

//I'm going to convert the nodes into edges and edge weights.
//I only make edges between nodes that are within the edgeRange
bool MST::nodesToEdges()
{
	if( (int)nodes.size() == 0)
		return false;

	std::vector<pairE> edgeArray;
	std::vector<float> edgeWeight;

	//int num_nodes = (int)nodes.size()-1;
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
				//Get the two intensity values:
				ImageType::IndexType index1;
				index1[0] = (int)n1.x;
				index1[1] = (int)n1.y;
				index1[2] = (int)n1.z;
				ImageType::IndexType index2;
				index2[0] = (int)n2.x;
				index2[1] = (int)n2.y;
				index2[2] = (int)n2.z;
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

	//Convert std c++ vector to c array:
	num_edges = (unsigned int)edgeArray.size();
	
	if(edge_array)
		delete[] edge_array;
	if(edge_wght)
		delete[] edge_wght;
	edge_array = new pairE[num_edges];
	edge_wght = new float[num_edges];

	for(unsigned int i = 0; i<num_edges; ++i)
	{
		edge_array[i] = edgeArray.at(i);
		edge_wght[i] = edgeWeight.at(i);
	}

	edgeWeight.clear();
	edgeArray.clear();

	if(debug)
		std::cerr << "Finished making edges = " << num_edges << std::endl;

	return true;
}

bool MST::minimumSpanningTree()
{
	if( num_edges <= 0 )
		return false;

	int num_nodes = (int)nodes.size();

	spanningTree.clear();

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph * g = new Graph(num_nodes);
	for (size_t j = 0; j < num_edges; ++j) 
	{
		Edge e; 
		bool inserted;
		tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, *g);
		weightmap[e] = edge_w[j];
	}
#else
	if(g)
		delete g;
	g = new Graph(edge_array, edge_array + num_edges, edge_wght, num_nodes);
#endif

	if(!g)
		return false;

	boost::property_map<Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, *g);
	
	// MST algorithm
	boost::kruskal_minimum_spanning_tree(*g, back_inserter(spanningTree));

	if(debug)
		std::cerr << "kruskal_minimum_spanning_tree(MST) is finished!" << std::endl;

	/*
	// Create a graph for the initial MST
	Graph msTree(num_nodes+1);
	int num_edge_MST = 0;
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		add_edge(source(*ei, *g), target(*ei, *g), msTree);
		num_edge_MST++;
	}
	if(debug)
		std::cerr << "MST edges = " << num_edge_MST << std::endl;
	*/

	return true;
}

bool MST::ErodeAndDialateNodeDegree(int mophStrength)
{
	if((int)spanningTree.size() == 0 || !g)
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
		degree_nodes[source(*ei, *g)] ++;
		degree_nodes[target(*ei, *g)] ++;
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
			if (degree_nodes_buffer[source(*ei, *g)]>0 && degree_nodes_buffer[target(*ei, *g)]>0)  
			{
				if (degree_nodes_buffer[source(*ei, *g)]==1 || degree_nodes_buffer[target(*ei, *g)]==1)  
				{
					degree_nodes[source(*ei, *g)] --;
					degree_nodes[target(*ei, *g)] --;
          
					// Save the edges eroded in a stack-like array. Each edge takes two elements of the array
					edge_eroded[num_edge_eroded*2]  = (int)source(*ei, *g); // Saving of eroded edges
					edge_eroded[num_edge_eroded*2+1] = (int)target(*ei, *g);
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

	if((int)spanningTree.size() == 0 || !g || (int)nodeDegree.size() == 0)
	{
		std::cerr << "missing tree or g or node degrees\n";
		return retLines;
	}

	//Create a msTree graph for backbone from the MST generated above
	//BackBone graph created or Graph for the parts of total skeletons 
	Graph msTreeBB(1);
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei)
	{
		if (nodeDegree.at(source(*ei, *g)-1)!=0 && nodeDegree.at(target(*ei, *g)-1)!=0)
		{  
			add_edge(source(*ei, *g), target(*ei, *g), msTreeBB);
		}
    }

	//Get the lines from the backbone mst:
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(msTreeBB); ei != ei_end; ++ei)
	{
		pairE ne( (int)source(*ei, msTreeBB)-1, (int)target(*ei, msTreeBB)-1 );
		retLines.push_back( ne );
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
	if(num_nodes == 0 || (int)spanningTree.size() == 0 || (int)nodeDegree.size() == 0)
		return retLines;

	// Create a graph for the initial MST
	Graph msTree(num_nodes+1);
	int num_edge_MST = 0;
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		add_edge(source(*ei, *g), target(*ei, *g), msTree);
		num_edge_MST++;
	}
	if(debug)
		std::cerr << "MST edges = " << num_edge_MST << std::endl;

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

	return retLines;
}

// this function is used to look the Minimum-Spanning Tree of the Skeletons 
Graph MST::CreateInitialmsTree(void)
{
    
	int num_nodes = (int)nodes.size();
    Graph msTree(num_nodes+1);
	if (num_nodes == 0 || (int)spanningTree.size() == 0 || (int)nodeDegree.size() == 0) 
		return msTree;

   	// Create a graph for the initial MST
	
	int num_edge_MST = 0;
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		add_edge(source(*ei, *g), target(*ei, *g), msTree);
		num_edge_MST++;
	}

    std::vector<pairE> retLines;

	//Get the lines from the backbone mst:
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(msTree); ei != ei_end; ++ei)
	{
		pairE ne( (int)source(*ei, msTree)-1, (int)target(*ei, msTree)-1 );
		retLines.push_back( ne );
   }

	if(debug)
		std::cerr << "Number of InitialSkeletonTree lines = " << retLines.size() << std::endl;

	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&nodes);
		fhdl->SetLines(&retLines);
		fhdl->Write("InitialMST.vtk");
		delete fhdl;
	}

	return msTree;

}

//  this function still has some problem 
Graph MST::morphGraphPrune(Graph msTree, int num_nodes,  float length_Threshold)
{

	
    Graph msTree_buffer(num_nodes+1);
	msTree_buffer = msTree; 
    typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter   ei, ei_end;
	typedef boost::graph_traits < Graph >::vertex_iterator Vertex_iter;
	Vertex_iter vi, vend;
	
	boost::graph_traits<Graph>::out_edge_iterator  outei, outedge_end;

	//std::vector<Point3D> vertexPos;

	int curBranchVerts[2000];

	int curBrVerts_Index = 0;
	float length_edge = 0; 
	int indVert, indVert_last;
	int num_leaves;
    int branchChosen;
	float length_leaf = 0;

	// Consider all leaves
	num_leaves = 0;

	  
    // for all vertex in the graph
	for(boost::tie(vi, vend) = vertices(msTree); vi != vend; ++vi) 
	 { 
			curBrVerts_Index = 0;
			if (out_degree(*vi, msTree) == 1)  
			   { // if it is a leaf tip
				 for (boost::tie(outei, outedge_end) = out_edges(*vi, msTree); outei != outedge_end; ++outei)
				   {
					curBranchVerts[curBrVerts_Index] = (int)source(*outei, msTree);
					curBrVerts_Index++;
					curBranchVerts[curBrVerts_Index] = (int)target(*outei, msTree);
				   } // end for

				  while (out_degree(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree) == 2) 
				   {
					//curVert = vertex(curBranchVerts[curBrVerts_Index], msTree);
					for (boost::tie(outei, outedge_end) = out_edges(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree);
																							outei != outedge_end; ++outei)
					  {
						if (target(*outei, msTree) == (unsigned int)curBranchVerts[curBrVerts_Index-1])
							continue;
						curBranchVerts[curBrVerts_Index+1] = (int)target(*outei, msTree);
					  } //end for
					  curBrVerts_Index++;
				   } // end while
				
				  branchChosen = 1;// Evaluate with MDL if the branch is chosen
				
				//if (i == prunetimes-1) num_leaves++;
                length_leaf = 0;  // must be reset the value,by xiaoliang
				for (int j = 0; j <= curBrVerts_Index; j++)
				{
					indVert = curBranchVerts[j];

					if (j==0) 
					 {
						indVert_last = indVert;
					 } // end if
                    
					Point3D vertexPosStart,vertexPosEnd;
					vertexPosStart.x = this->nodes.at(indVert).x;
					vertexPosStart.y = this->nodes.at(indVert).y;
                    vertexPosStart.z = this->nodes.at(indVert).z;
					vertexPosEnd.x = this->nodes.at(indVert_last).x;
                    vertexPosEnd.y = this->nodes.at(indVert_last).y;
					vertexPosEnd.z = this->nodes.at(indVert_last).z;

                    length_edge = (vertexPosStart.x-vertexPosEnd.x)*(vertexPosStart.x-vertexPosEnd.x);
					length_edge+= (vertexPosStart.y-vertexPosEnd.y)*(vertexPosStart.y-vertexPosEnd.y);
					length_edge+= (vertexPosStart.z-vertexPosEnd.z)*(vertexPosStart.z-vertexPosEnd.z);
					/*
					length_edge = (vertexPos[indVert].x-vertexPos[indVert_last].x)*(vertexPos[indVert].x-vertexPos[indVert_last].x);
					length_edge+= (vertexPos[indVert].y-vertexPos[indVert_last].y)*(vertexPos[indVert].y-vertexPos[indVert_last].y);
					length_edge+= (vertexPos[indVert].z-vertexPos[indVert_last].z)*(vertexPos[indVert].z-vertexPos[indVert_last].z);
                     */
					length_edge = sqrt(length_edge);
					length_leaf += length_edge;
					indVert_last = indVert;
					if(debug)
						std::cout << "length_leaf = " << length_leaf << std::endl;
				} //end for

				if(debug)
					std::cout << "I am here" << std::endl;
				if (length_leaf < length_Threshold)    
				{  
					branchChosen = 0;     
				}
				// Remove the branch based on leaf-length critierion
				if (branchChosen == 0)  {
					for (int j = 1; j <= curBrVerts_Index; j++) {
						remove_edge(curBranchVerts[j-1], curBranchVerts[j], msTree_buffer);
					}
				}

			}
		}
	
	msTree = msTree_buffer;

	if (debug)
		std::cerr<<"The Graph Pruning is done with respect to the threshold= %f\n" << length_Threshold;
    
    if (debug)
	{
     std::vector<pairE> retLines;
	 typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	 Edge_iter ei, ei_end; 
	 for (tie(ei, ei_end) = edges(msTree); ei != ei_end; ++ei)
	 {
		pairE ne( (int)source(*ei, msTree)-1, (int)target(*ei, msTree)-1 );
		retLines.push_back( ne );
     }

	 vtkFileHandler * fhdl = new vtkFileHandler();
	 fhdl->SetNodes(&nodes);
	 fhdl->SetLines(&retLines);
	 fhdl->Write("PruningMST.vtk");
	 delete fhdl;
	}
	
	return msTree;
     
}


}