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
// morphGraphPrune.cpp : Morphological process on graph 
// Function: Pruning on Minimum Spanning Tree for trivia branches
//           Removing branches with length below certain threshold
// Author: Xiaosong Yuan, RPI, modified by Xiao L


#include "MST.h"


Graph morphGraphPrune(Graph msTree, int num_nodes, struct VoxelPosition *vertexPos, float length_Threshold) 
{

	Graph msTree_buffer(num_nodes+1);
    // copy to buffer for next process
	msTree_buffer = msTree; 
	graph_traits<Graph>::out_edge_iterator  outei, outedge_end;

	Edge_iter   ei, ei_end;
	Vertex_iter vi, vend;
	int curBranchVerts[2000];
	int indVert, indVert_last;
    int branchChosen;
	int prunetimes = 1; 
	  
	for (int i=0; i< prunetimes; i++) {
		for(boost::tie(vi, vend) = vertices(msTree); vi != vend; ++vi) 
		 { // for all vertex in the graph
			int curBrVerts_Index = 0;
			if (out_degree(*vi, msTree) == 1)  
			   { // if it is a leaf tip
				 for (boost::tie(outei, outedge_end) = out_edges(*vi, msTree); outei != outedge_end; ++outei)
				   {
					curBranchVerts[curBrVerts_Index] = source(*outei, msTree);
					curBrVerts_Index++;
					curBranchVerts[curBrVerts_Index] = target(*outei, msTree);
				   } // end for

				  while (out_degree(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree) == 2) 
				   {
					//curVert = vertex(curBranchVerts[curBrVerts_Index], msTree);
					for (boost::tie(outei, outedge_end) = out_edges(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree);
																							outei != outedge_end; ++outei)
					  {
						if (target(*outei, msTree) == (unsigned int)curBranchVerts[curBrVerts_Index-1])
							continue;
						curBranchVerts[curBrVerts_Index+1] = target(*outei, msTree);
					  } //end for
					  curBrVerts_Index++;
				   } // end while
				  
				branchChosen = 1;
                float length_leaf = 0;  // must be reset the value,by xiaoliang
				for (int j = 0; j <= curBrVerts_Index; j++)
				{
					indVert = curBranchVerts[j];
					if(indVert > num_nodes)
				    {
					 std::cerr << "not a node\n";
					 indVert = indVert_last;
				    }

					if (j==0) 
					 {
						indVert_last = indVert;
					  } // end if
					float length_edge = (vertexPos[indVert].x-vertexPos[indVert_last].x)*(vertexPos[indVert].x-vertexPos[indVert_last].x);
					length_edge+= (vertexPos[indVert].y-vertexPos[indVert_last].y)*(vertexPos[indVert].y-vertexPos[indVert_last].y);
					length_edge+= (vertexPos[indVert].z-vertexPos[indVert_last].z)*(vertexPos[indVert].z-vertexPos[indVert_last].z);
					length_edge = sqrt(length_edge);
					length_leaf += length_edge;
					indVert_last = indVert;	
				} //end for
	          /*
				mahalanobis_dist = 0;
				x1 = meanDensityBranch - 30.81;  // minus mean_feature from matlab 
				x2 = length_leaf - 12.90;
				if (x2>10)  x2=10; //keep long dendrite 
				x3 = meanVesselBranch - 215.53;
				mahalanobis_dist = x1*x1*0.0136+ 2*x1*x2*(-0.0030)+ 2*x1*x3*(-0.0008)+ x2*x2*0.0286+ 2*x2*x3*0.0012
								+x3*x3*0.0005;
				if (mahalanobis_dist > 1.96*1.96)  branchChosen = 0;   // At 5% level of significance
	         */
				//if (meanDensityBranch < 150)     branchChosen = 0;
				//if (length_leaf > 25)            branchChosen = 1;
              
				if (length_leaf < length_Threshold)    
				{  
					branchChosen = 0;     // Parameter length_leaf to choose -> 2, 3, 5
				}

				// Remove the branch based on MDL critierion
				if (branchChosen == 0) 
				{
					for (int j = 1; j <= curBrVerts_Index; j++)
					{
						remove_edge(curBranchVerts[j-1], curBranchVerts[j], msTree_buffer);
					}// end for
				}// end if
			}//if (out_degree(*vi, msTree) == 1)  
		}// end for(boost::tie(vi, vend) = vertices(msTree); vi != vend; ++vi) 
		msTree = msTree_buffer;
	}  // end of prunetimes

	std::cout << "The Graph Pruning is done with respect to the threshold " << length_Threshold << std::endl;
	return msTree;

}
