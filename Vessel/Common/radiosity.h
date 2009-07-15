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

#ifndef _RADIOSITY_H_
#define _RADIOSITY_H_

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 
#include "array.h"
#include "vectors.h"
#include "argparser.h"

class Mesh;
class Face;
class Vertex;


class Radiosity {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Radiosity(Mesh *m, ArgParser *args);
  ~Radiosity();
  void Reset();
  void Initialize();
  void Cleanup();

  // =========
  // ACCESSORS
  Mesh* getMesh() const { return mesh; }
  Face* getPatch(int i) const {
    assert (i >= 0 && i < n);
    return patches[i]; }
  float getFormFactor(int i, int j) const {
    // F_i,j radiant energy leaving i arriving at j
    assert (i >= 0 && i < n);
    assert (j >= 0 && j < n);
    return 0; // CHANGED
    return formfactors[i*n+j]; }
  float getArea(int i) const {
    assert (i >= 0 && i < n);
    return area[i]; }
  Vec3f getUndistributed(int i) const {
    assert (i >= 0 && i < n);
    return undistributed[i]; }
  Vec3f getAbsorbed(int i) const {
    assert (i >= 0 && i < n);
    return absorbed[i]; }
  Vec3f getRadiance(int i) const {
    assert (i >= 0 && i < n);
    return radiance[i]; }

  // =========
  // MODIFIERS
  bool Iterate();
  void setFormFactor(int i, int j, float value) { 
    assert (i >= 0 && i < n);
    assert (j >= 0 && j < n);
    formfactors[i*n+j] = value; }
  void normalizeFormFactors(int i) {
    float sum = 0;
    int j;
    for (j = 0; j < n; j++) {
      sum += getFormFactor(i,j); }
    assert(sum > 0);
    for (j = 0; j < n; j++) {
      setFormFactor(i,j,getFormFactor(i,j)/sum); } }
  void setArea(int i, float value) {
    assert (i >= 0 && i < n);
    area[i] = value; }
  void setUndistributed(int i, Vec3f value) { 
    assert (i >= 0 && i < n);
    undistributed[i] = value; }
  void findMaxUndistributed();
  void setAbsorbed(int i, Vec3f value) { 
    assert (i >= 0 && i < n);
    absorbed[i] = value; }
  void setRadiance(int i, Vec3f value) { 
    assert (i >= 0 && i < n);
    radiance[i] = value; }
	
  // =====
  // PAINT
  void Paint(ArgParser *args);
  void PaintSelection(ArgParser *args);
  Vec3f whichVisualization(enum RENDER_MODE mode, Face *f, int i);
  void insertInterpolatedColor(int index, Face *f, Vertex *v);
  void insertColor(Vec3f v);
	Face * getFaceFromID(int);
  // ===============
  // HELPER FUNCIONS
  void computeFormFactors(int i);

private:

  // ==============
  // REPRESENTATION
  Mesh *mesh;
  ArgParser *args;

  // num patches
  int n; 
  Face* *patches;

  // a nxn matrix
  // F_i,j radiant energy leaving i arriving at j
  float *formfactors;

  // length n vectors
  float *area;
  Vec3f *undistributed; // energy per unit area
  Vec3f *absorbed;      // energy per unit area
  Vec3f *radiance;      // energy per unit area

  int max_undistributed_patch;  // the patch with the most undistributed energy
  float total_undistributed;    // the total amount of undistributed light
  float total_area;             // the total area of the scene
  bool volume_rendering_setup_complete;
};

#endif
