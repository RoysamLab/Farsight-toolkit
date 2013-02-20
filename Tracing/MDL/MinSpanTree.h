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

// Author: Xiaosong Yuan, RPI
// Date: May 30th, 2007
// Head file

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <iostream>
#include <fstream>
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

