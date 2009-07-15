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

#include <math.h>
#include "camera.h"
#include "matrix.h"
#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 

// ====================================================================
// ====================================================================
// CONSTRUCTORS

#ifndef M_PI
#define M_PI (2*acos(double(0)))
#endif
Camera::Camera(Vec3f c, Vec3f d, Vec3f u) {
  center = c;
  direction = d;
  up = u;
  // normalize the vectors
  up.Normalize();
  direction.Normalize();
}

PerspectiveCamera::PerspectiveCamera(Vec3f c, Vec3f d, Vec3f u, float a) : Camera(c,d,u) {
  angle = a;

  Vec3f screenCenter = center + direction;
  float screenHeight = tan(angle/2.0);

  lowerLeft = screenCenter - (getScreenUp() * screenHeight) - (getHorizontal() * screenHeight);
  xAxis = getHorizontal() * 2 * screenHeight;
  yAxis = getScreenUp() * 2 * screenHeight;
}

// ====================================================================
// ====================================================================
// GL INIT
// Create a camera with the appropriate dimensions that
// crops the screen in the narrowest dimension.

void PerspectiveCamera::glInit(int w, int h) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  float aspect = float(w)/float(h);
  float asp_angle = angle * 180/M_PI;
  if (aspect > 1) asp_angle /= aspect;

  gluPerspective(asp_angle, aspect, 1, 100.0);
}
void PerspectiveCamera::glInit_forselection(int w,int h,int x, int y){
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  gluPickMatrix(x,viewport[3]-y-1,1,1,viewport); 

  float aspect = float(w)/float(h);
  float asp_angle = angle * 180/M_PI;
  if (aspect > 1) asp_angle /= aspect;

  gluPerspective(asp_angle, aspect, 1, 100.0);
}
// ====================================================================
// ====================================================================
// GL PLACE CAMERA
// Place a camera within an OpenGL scene

void Camera::glPlaceCamera(void) {
  Vec3f lookAt = center + direction;
  gluLookAt(center.x(), center.y(), center.z(),
	    lookAt.x(), lookAt.y(), lookAt.z(),
            up.x(), up.y(), up.z());
}

// ====================================================================
// dollyCamera, truckCamera, and RotateCamera
//
// Asumptions:
//  - up is really up (i.e., it hasn't been changed
//    to point to "screen up")
//  - up and direction are normalized
// Special considerations:
//  - If your constructor precomputes any vectors for
//    use in 'generateRay', you will likely to recompute those
//    values at athe end of the these three routines
// ====================================================================

// ====================================================================
// dollyCamera: Move camera along the direction vector
// ====================================================================

void PerspectiveCamera::dollyCamera(float dist) {
  center += direction*dist;
  Vec3f screenCenter = center + direction;
  float screenHeight = tan(angle/2.0);
  lowerLeft = screenCenter - (getScreenUp() * screenHeight) - (getHorizontal() * screenHeight);
}

// ====================================================================
// truckCamera: Translate camera perpendicular to the direction vector
// ====================================================================

void PerspectiveCamera::truckCamera(float dx, float dy) {
  center += getHorizontal()*dx + getScreenUp()*dy;
  Vec3f screenCenter = center + direction;
  float screenHeight = tan(angle/2.0);
  lowerLeft = screenCenter - (getScreenUp() * screenHeight) - (getHorizontal() * screenHeight);
}

// ====================================================================
// rotateCamera: Rotate around the up and horizontal vectors
// ====================================================================

void PerspectiveCamera::rotateCamera(float rx, float ry) {
  // Don't let the model flip upside-down (There is a singularity
  // at the poles when 'up' and 'direction' are aligned)
  float tiltAngle = acos(up.Dot3(direction));
  if (tiltAngle-ry > 3.13)
    ry = tiltAngle - 3.13;
  else if (tiltAngle-ry < 0.01)
    ry = tiltAngle - 0.01;
    
// ARUN CHANGED THIS!
  Matrix rotMat = Matrix::MakeAxisRotation(up, rx);
  rotMat *= Matrix::MakeAxisRotation(getHorizontal(), ry);

  rotMat.Transform(center);
  rotMat.TransformDirection(direction);
  direction.Normalize();
  Vec3f screenCenter = center + direction;
  float screenHeight = tan(angle/2.0);

  lowerLeft = screenCenter - (getScreenUp() * screenHeight) - (getHorizontal() * screenHeight);
  xAxis = getHorizontal() * 2 * screenHeight;
  yAxis = getScreenUp() * 2 * screenHeight;
}

// ====================================================================
// ====================================================================
