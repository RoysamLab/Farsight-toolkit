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

#ifndef _EDGE_H_
#define _EDGE_H_

#include "vertex.h"
#include "edge.h"

Edge::Edge(Vertex *v, Face *f) {
  vertex = v;
  face = f;
  next = NULL;
  opposite = NULL;
}

Edge::~Edge() { 
//	clearOpposite();
  if (opposite != NULL)
    opposite->opposite = NULL;
}

void Edge::extract_func(Edge *e, int &a, int &b, int &c, int &d) {
  a = (*e)[0]->getIndex();
  b = (*e)[1]->getIndex();
  c = 0;
  d = 0;
}

void Edge::Print() { 
  printf ("EDGE %d -> %d\n", (*this)[0]->getIndex(), (*this)[1]->getIndex()); 
}

#endif
