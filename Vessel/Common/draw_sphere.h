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

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif 
//#include <GL/glx.h>
#include <stdio.h>
#include <sys/types.h>
//#include <unistd.h>

void Normalize(float v[3]);
void DrawTriangle(float *v1, float *v2, float *v3);
void Subdivide(float *v1, float *v2, float *v3, int depth);
void DrawSphere(float radius, int detail);
