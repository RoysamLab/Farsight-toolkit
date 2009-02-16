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
