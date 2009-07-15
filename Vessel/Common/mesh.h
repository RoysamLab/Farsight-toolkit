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

#ifndef MESH_H
#define MESH_H
#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 
#include "vectors.h"
#include "array.h"
#include "bag.h"
#include "boundingbox.h"
#include "argparser.h"

#include <vector>

class Vertex;
class Edge;
class Face;
class VertexParent;

struct trio{
	int a, b,c;
	char annotation[512];
	int type;
	int edit_number;
};
// ======================================================================
// ======================================================================

class Mesh {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Mesh();
  virtual ~Mesh();
  void Load(const char *input_file);
    
  // ========
  // VERTICES
  int numVertices() const { return vertices->Count(); }
  Vertex* addVertex(const Vec3f &pos);
  // this creates a relationship between 3 vertices (2 parents, 1 child)
  void setParentsChild(Vertex *p1, Vertex *p2, Vertex *child);
  // this accessor will find a child vertex (if it exists) when given
  // two parent vertices
  Vertex* getChildVertex(Vertex *p1, Vertex *p2) const;
  // look up vertex by index from original .obj file
  Vertex* getVertex(int i) const {
    assert (i >= 0 && i < numVertices());
    Vertex *v = (*vertices)[i];
    assert (v != NULL);
    return v; }

  // =====
  // EDGES
  int numEdges() const { return edges->Count(); }
  // this efficiently looks for an edge with the given vertices, using a hash table
  Edge* getEdge(Vertex *a, Vertex *b) const;

  // =========
  // FACES
  int numFaces() const { return faces->Count(); }
  void addFace(Vertex *a, Vertex *b, Vertex *c, const Vec3f &color, const Vec3f &emit);
  Face * addFace2(Vertex *a, Vertex *b, Vertex *c, const Vec3f &color, const Vec3f &emit);
  void removeFace(Face *f);
  void removeFace2(Face *f);
  // ===============
  // OTHER ACCESSORS
  BoundingBox* getBoundingBox() const { return bbox; }
  Bag<Face*>* getFaces() const { return faces; }
  Bag<Edge*>* getEdges() const { return edges; }
  Array<Vertex*>* getVertices() const {return vertices;}
  // ===============
  // OTHER FUNCTIONS
  void PaintWireframe();
  void Subdivision();
  void PaintPoints();
  void PaintVotes(int);
  void PaintConfidence();
  void PaintCenterLines();
  void texture_list_init();
  void set_pic_name(char *);
 GLuint *texturexy, *textureyz, *texturexz;
 GLuint *texture;
 int rwidth, rlength, rdepth;
 GLuint list;
 GLuint box;
  vector<trio> undo_operations;
  vector<trio> redo_operations;
private:

  // helper functions
  Vertex* AddEdgeVertex(Vertex *a, Vertex *b);
  Vertex* AddMidVertex(Vertex *a, Vertex *b, Vertex *c);
  char pic_name[1024];
  // ==============
  // REPRESENTATION

  Array<Vertex*> *vertices;
  Bag<Edge*> *edges;
  Bag<Face*> *faces;
  BoundingBox *bbox;
  Bag<VertexParent*> *vertex_parents;

};

// ======================================================================
// ======================================================================

#endif




