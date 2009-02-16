// Author: Xiaosong Yuan, RPI
// Date: May 30th, 2007
// Head file

#include <iostream>
//#include "morphGraphPrune.h"
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
typedef std::pair<int, int>  E;

typedef graph_traits < Graph >::edge_iterator Edge_iter;  
typedef graph_traits < Graph >::vertex_iterator Vertex_iter;




// Head file: MinSpanTree.h

/*
struct  VoxelPosition
{
	float x;
	float y;
	float z;
};
*/
