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

#ifndef _FACE_H
#define _FACE_H

#include <limits.h>
#include "boundingbox.h"
#include "edge.h"
#include "vertex.h"
#include "vectors.h"
// ===========================================================

class Face {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR

  Face(const Vec3f &c, const Vec3f &e) {
    edge = NULL;
    color = c; 
    emit = e; 
    mark = false;
  }
  ~Face() {}

  // here's the hash function to use for faces so they
  // can be efficiently accessed within the Bag data structure
  static void extract_func(Face *t, int &a, int &b, int &c, int &d) {
    a = (*t)[0]->getIndex(); 
    b = (*t)[1]->getIndex(); 
    c = (*t)[2]->getIndex();
    d = 0;
   // d = (*t)[3]->getIndex();
  }

  // =========
  // ACCESSORS
  Vertex* operator[](int i) const { 
    assert (edge != NULL);
    if (i==0) return edge->getVertex();
    if (i==1) return edge->getNext()->getVertex();
    if (i==2) return edge->getNext()->getNext()->getVertex();
   // if (i==3) return edge->getNext()->getNext()->getNext()->getVertex();
    assert(0);
  }
  Edge* getEdge() const { 
    assert (edge != NULL);
    return edge; 
  }
  void setMark(){mark = true;}
  bool getMark(){return mark;}
  void unMark(){mark = false;}

  // =========
  // MODIFIERS
  void setEdge(Edge *e) {
    assert (edge == NULL);
    edge = e;
  }

  // ===================
  // RADIOSITY ACCESSORS
  Vec3f getColor() const { return color; }
  Vec3f getEmit() const { return emit; }
  void setColor(Vec3f c) { color = c; }
  void setEmit(Vec3f e ) { emit = e; }
  float getArea() const;
  Vec3f getNormal() const{
		
		Vec3f v1 = (*edge)[0]->get() - (*edge)[1]->get();
		Vec3f v2 = (*(edge->getNext()))[1]->get() - (*(edge->getNext()))[0]->get();
		Vec3f res;
		Vec3f::Cross3(res,v1,v2);
		return res;
		
		}
  Vec3f getCenter() const 
{

	Edge *e = edge;
	Vec3f center = Vec3f(0,0,0);
	for(int co = 0; co <3; co++)
	{
		center = center + e->getVertex()->get();
		e = e->getNext();
	}
	center = center/3.0;
	return center;
}
  int getRadiosityPatchIndex() const { return radiosity_patch_index; }

  // ===================
  // RADIOSITY MODIFIERS
  void setRadiosityPatchIndex(int i) { radiosity_patch_index = i; }

  // NOTE: If you want to modify a face, remove it from the mesh,
  // delete it, create a new copy with the changes, and re-add it.
  // This will ensure the edges get updated appropriately.
Face * opposite;

protected:

  // don't use this constructor
  Face& operator = (const Face &f) { assert(0); }
  
  // ==============
  // REPRESENTATION
  Edge *edge;

  
  // radiosity variables
  Vec3f emit;   // the energy emitted per unit of patch area
  Vec3f color;  // the color (diffuse reflectivity, BRDF)
  int radiosity_patch_index;  // an awkward pointer to this patch in the Radiosity patch array
    bool mark;
};

// ===========================================================

#endif
