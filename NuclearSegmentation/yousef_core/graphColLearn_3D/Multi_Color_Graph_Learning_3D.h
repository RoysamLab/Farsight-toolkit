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

#ifndef Multi_Color_Graph_Learning_3D_H
#define Multi_Color_Graph_Learning_3D_H

#include <vector>
#include <iostream>
#include <math.h>

float* multiColGraphLearning(float* X_vals, int* labs_vals, int* color_im, int r, int c, int z, int *NC, int refinemetRange);
//added by Yousef on 11/3/2008
void distToEdge(int *** edge_im, int R, int C, int Z);
//////////////////////////////

#endif
