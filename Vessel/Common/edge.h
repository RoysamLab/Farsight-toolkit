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

#ifndef EDGE_H
#define EDGE_H

#include <limits.h>
#include <stdio.h>
#include <assert.h>

class Vertex;
class Face;

// ===================================================================
// half-edge data structure

class Edge { 

public:

  // ========================
  // CONSTRUCTORS & DESTRUCTOR
  Edge(Vertex *v, Face *f);
  ~Edge();
//void ClearOpposite();
  // here's the hash function to use for edges so they
  // can be efficiently accessed within the Bag data structure
  static void extract_func(Edge *e, int &a, int &b, int &c, int &d);
 
  // =========
  // ACCESSORS
  Vertex* getVertex() const { assert (vertex != NULL); return vertex; }
  Edge* getNext() const { assert (next != NULL); return next; }
  Edge* getPrev() const { assert (next != NULL); assert(next->getNext()!=NULL); return next->getNext();}
  Face* getFace() const { 
		if(face == NULL ) 
		for(;;); 
	  assert (face != NULL);
	   return face; }
  Edge* getOpposite() const {
    // warning!  the opposite edge might be NULL!
    return opposite; }
	void setVertex( Vertex *v ){ vertex =v;	}
  Vertex* operator[](int i) const { 
    if (i==0) return getVertex();
    if (i==1) return getNext()->getVertex();
    assert(0);
  }

  // =========
  // MODIFIERS
  void setOpposite(Edge *e) {
		if(opposite!=NULL)
			{
				printf("opp ! null\n");
			//	for(;;);
			}
    //assert (opposite == NULL); 
    		if(e==NULL)
			{
				printf("e == null\n");
				for(;;);
			}
    assert (e != NULL);
    		if(e->opposite!=NULL)
			{
				printf("e->opp ! null\n");
				for(;;);
			}
    assert (e->opposite == NULL);
    opposite = e; 
    e->opposite = this; 
  }
  void clearOpposite() { 
    if (opposite == NULL) return; 
    assert (opposite->opposite == this); 
    opposite->opposite = NULL;
    opposite = NULL; 
  }
  void setNext(Edge *e) {
    assert (next == NULL);
    assert (e != NULL);
    assert (face == e->face);
    next = e;
  }

  void Print();
  
private:

  Edge(const Edge&) { assert(0); }
  //Edge& operator=(const Edge&) { assert(0); }

  

  // ==============
  // REPRESENTATION
  // in the half edge data adjacency data structure, the edge stores everything!
  Vertex *vertex;
  Face *face;
  Edge *opposite;
  Edge *next;
};

// ===================================================================

#endif
