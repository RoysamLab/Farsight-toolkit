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

#include "boundingbox.h"
#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 

// ====================================================================
// ====================================================================

void BoundingBox::paint() const {
  // draw a wireframe box to represent the boundingbox
  
  glColor3f(1,1,1);
  glLineWidth(1);
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);

  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(min.x(),max.y(),min.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(min.x(),max.y(),min.z());

  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(min.x(),max.y(),min.z());
  glVertex3f(min.x(),max.y(),max.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(max.x(),max.y(),max.z());

  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(min.x(),max.y(),max.z());
  glVertex3f(max.x(),max.y(),max.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(max.x(),max.y(),max.z());
  glVertex3f(min.x(),max.y(),max.z());

  glEnd();
  glEnable(GL_LIGHTING);	   
}

// ====================================================================
// ====================================================================
