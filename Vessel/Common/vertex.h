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

#ifndef _VERTEX_H
#define _VERTEX_H

#include <stdio.h>
#include <assert.h>

#include "vectors.h"

class Vertex;
//#pragma warning(disable : 4244)
// ==========================================================

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const Vec3f &pos) : position(pos) { index = i; mark = false; }
  virtual ~Vertex() { }
  
  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x(); }
  double y() const { return position.y(); }
  double z() const { return position.z(); }
  const Vec3f& get() const { return position; }
  bool getMark() { return mark;}
  void negateMark() { mark = !mark;}
  void print(){printf(" x %f y %f z %f ",position.x(),position.y(),position.z());}
  // =========
  // MODIFIERS
  void set(Vec3f v) { position = v; }
  void set(double x, double y, double z) { position.Set(x,y,z); }
  void setNewPos(Vec3f v) { newpos = v;}
  void updatePos(){position = newpos;}
  Vec3f getNewPos(){return newpos;}
   Vec3f inc;
   int inccounts;
private:

  // don't use these constructors
  Vertex() { assert(0); }
  //Vertex& operator=(const Vertex&) { assert(0); }
  Vertex(const Vertex&) { assert(0); }
  
  // ==============
  // REPRESENTATION
  Vec3f position;
  Vec3f newpos;
//  Vec3f inc;
  
  bool mark;
  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure
  int index;  

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation.

};

// ==========================================================

#endif

