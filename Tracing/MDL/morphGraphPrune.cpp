// morphGraphPrune.cpp : Morphological process on graph 
// Function: Pruning on Minimum Spanning Tree for trivia branches
//           Removing branches with length below certain threshold
// Author: Xiaosong Yuan, RPI
// Date: May 30th, 2007

#include "MinSpanTree.h"
#include "morphGraphPrune.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <iostream>
#include <fstream>
//#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <boost/config.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list <vecS, vecS, undirectedS, no_property, property <edge_weight_t, float> >  Graph;
typedef graph_traits < Graph >::edge_descriptor  Edge;
typedef graph_traits < Graph >::vertex_descriptor  Vertex;

typedef graph_traits < Graph >::edge_iterator Edge_iter; 
typedef graph_traits < Graph >::vertex_iterator Vertex_iter;

struct  VoxelPosition
{
	float x;
	float y;
	float z;
};

// Arguments of function: degree_nodes[], tree graph, vertexPos[]
// Output to the same Graph: msTree (after remove_edge()).

Graph morphGraphPrune(Graph msTree, int num_nodes, struct VoxelPosition *vertexPos, float length_leaf) 
{
	//struct VoxelPosition *vertexPos;
	Graph msTree_buffer(num_nodes+1);
	//int num_edge_MST = 0;
	//int num_vertex_MST = num_vertices(msTree); 
	//int num_edge_MST = num_edges(msTree); //not work

	msTree_buffer = msTree; // copy to buffer for next process

	Edge_iter   ei, ei_end;
	Vertex_iter vi, vend;
	graph_traits<Graph>::out_edge_iterator  outei, outedge_end;
	int curBranchVerts[2000];
	int curBrVerts_Index = 0;
	float length_edge;
	int indVert, indVert_last;
	int num_leaves;
	int i, j;
    int branchChosen;


	// Consider all leaves
	//Vertex curVert;
	num_leaves = 0;
	int prunetimes = 1; 
	  
	for (i=0; i< prunetimes; i++) {
		for(boost::tie(vi, vend) = vertices(msTree); vi != vend; ++vi)  { // for all vertex in the graph
			curBrVerts_Index = 0;
			if (out_degree(*vi, msTree) == 1)  { // if it is a leaf tip
				for (boost::tie(outei, outedge_end) = out_edges(*vi, msTree); outei != outedge_end; ++outei) {
					curBranchVerts[curBrVerts_Index] = source(*outei, msTree);
					curBrVerts_Index++;
					curBranchVerts[curBrVerts_Index] = target(*outei, msTree);
				}
				while (out_degree(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree) == 2) {
					//curVert = vertex(curBranchVerts[curBrVerts_Index], msTree);
					for (boost::tie(outei, outedge_end) = out_edges(vertex(curBranchVerts[curBrVerts_Index], msTree), msTree);
																							outei != outedge_end; ++outei) {
						if (target(*outei, msTree) == (unsigned int)curBranchVerts[curBrVerts_Index-1])
							continue;
						curBranchVerts[curBrVerts_Index+1] = target(*outei, msTree);
					}
					curBrVerts_Index++;
				}
				branchChosen = 1;
				// Evaluate with MDL if the branch is chosen
				if (i == prunetimes-1) num_leaves++;
				for (j = 0; j <= curBrVerts_Index; j++) {
					indVert = curBranchVerts[j];

					if (j==0) {
						indVert_last = indVert;
					}
					length_edge = (vertexPos[indVert].x-vertexPos[indVert_last].x)*(vertexPos[indVert].x-vertexPos[indVert_last].x);
					length_edge+= (vertexPos[indVert].y-vertexPos[indVert_last].y)*(vertexPos[indVert].y-vertexPos[indVert_last].y);
					length_edge+= (vertexPos[indVert].z-vertexPos[indVert_last].z)*(vertexPos[indVert].z-vertexPos[indVert_last].z);
					length_edge = sqrt(length_edge);
					length_leaf += length_edge;
					indVert_last = indVert;

				}
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

				if (length_leaf < 4)      branchChosen = 0;     // Parameter length_leaf to choose -> 2, 3, 5

				// Remove the branch based on MDL critierion
				if (branchChosen == 0)  {
					for (j = 1; j <= curBrVerts_Index; j++) {
						remove_edge(curBranchVerts[j-1], curBranchVerts[j], msTree_buffer);
					}
				}

			}
		}
		msTree = msTree_buffer;
	}  // end of prunetimes


	return msTree;

}
