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
