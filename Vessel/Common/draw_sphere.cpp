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

#include "draw_sphere.h"
#include <math.h>

/*********** Draw a sphere (based on an OpenGL example) ***************/

#define X .525731112119133606
#define Z .850650808352039932

void Normalize(float v[3]) {
  GLfloat d = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (d==0.0) // Error: Zero length vector
    return;
  v[0] /= d;  v[1] /= d;  v[2] /= d;
}

void DrawTriangle(float *v1, float *v2, float *v3) {
  glBegin(GL_POLYGON);
    glNormal3fv(v1); glVertex3fv(v1);
    glNormal3fv(v2); glVertex3fv(v2);
    glNormal3fv(v3); glVertex3fv(v3);
  glEnd();
}

void Subdivide(float *v1, float *v2, float *v3, int depth) {
  GLfloat v12[3], v23[3], v31[3];
  GLint i;

  if (depth==0) {
    DrawTriangle(v1, v2, v3);
    return;
  }
  for (i=0; i<3; i++) {
    v12[i] = v1[i]+v2[i];
    v23[i] = v2[i]+v3[i];
    v31[i] = v3[i]+v1[i];
  }
  Normalize(v12);
  Normalize(v23);
  Normalize(v31);
  Subdivide(v1, v12, v31, depth-1);
  Subdivide(v2, v23, v12, depth-1);
  Subdivide(v3, v31, v23, depth-1);
  Subdivide(v12, v23, v31, depth-1);
}

void DrawSphere(float radius, int detail) {
  static GLfloat vdata[12][3] = {
    {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
    {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
    {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
  };

  static GLint tindices[20][3] = {
    {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
    {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
    {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
    {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
  };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glScalef(radius, radius, radius);

  for (int i=0; i<20; i++) {
    Subdivide(&vdata[tindices[i][0]][0],
              &vdata[tindices[i][1]][0],
              &vdata[tindices[i][2]][0], detail);
  }

  glPopMatrix();
}
/**********************************************************************/

