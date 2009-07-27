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

#ifndef _CAMERA_H_
#define _CAMERA_H_

#include <assert.h>
#include "vectors.h"

// ====================================================================

class Camera {

public:
  // CONSTRUCTOR & DESTRUCTOR
  Camera(Vec3f c, Vec3f d, Vec3f u);
 virtual ~Camera() {}

  // GL NAVIGATION
  virtual void glInit(int w, int h) = 0;
  virtual  void glInit_forselection(int w,int h,int x, int y)=0;
  void glPlaceCamera(void);
  virtual void dollyCamera(float dist) = 0;
  virtual void truckCamera(float dx, float dy) = 0;
  virtual void rotateCamera(float rx, float ry) = 0;
  virtual void Print() = 0;
  virtual float getAngle()=0;;
protected:
  Camera() { assert(0); } // don't use

  // HELPER FUNCTIONS
  const Vec3f getHorizontal() const {
    Vec3f answer;
    Vec3f::Cross3(answer, direction, up);
    answer.Normalize();
    return answer; }
  const Vec3f getScreenUp() const {
    Vec3f answer;
    Vec3f::Cross3(answer, getHorizontal(), direction);
    return answer; }

  // REPRESENTATION
  Vec3f center;
  Vec3f direction;
  Vec3f up;
};

// ====================================================================

class PerspectiveCamera : public Camera {

public:
  // CONSTRUCTOR & DESTRUCTOR
  PerspectiveCamera(Vec3f c, Vec3f d, Vec3f u, float a);
  virtual ~PerspectiveCamera(void) { }

  // GL NAVIGATION
  void glInit(int w, int h);
  void glInit_forselection(int w,int h,int x, int y);
  void dollyCamera(float dist);
  void truckCamera(float dx, float dy);
  void rotateCamera(float rx, float ry);
  void Print() {
    printf ("PerspectiveCamera {\n");
    printf ("    center    ");
    center.Write(stdout);
    printf ("    direction ");
    direction.Write(stdout);
    printf ("    up        ");
    up.Write(stdout);
    printf ("    angle      %f\n", angle);
    printf ("}\n");
  }    
  float getAngle()
  {
	  return angle;
  }

private:
  PerspectiveCamera() { assert(0); } // don't use

  // REPRESENTATION
  float angle;
  Vec3f lowerLeft;
  Vec3f xAxis;
  Vec3f yAxis;
};

// ====================================================================

#endif
