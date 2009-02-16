#ifndef CELL_BINARIZATION_H
#define CELL_BINARIZATION_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<limits>
#include<exception>
#include <algorithm>

#include "new_graph.h"
#include "itkMinErrorThresholdImageFilter.h"

int Cell_Binarization_2D(unsigned char* imgIn, int *imgOut, int R, int C, int shd); 
int Cell_Binarization_3D(unsigned char *imgIn, int* imgOut, int R, int C, int Z, int shd); 
double compute_poisson_prob(double intensity, double alpha);
void Seg_GC_Full_2D(unsigned char* IM, int r, int c, double alpha_F, double alpha_B, double P_I, int* num_nodes, int* num_edges, int* Seg_out);
void Seg_GC_Full_3D(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, int* Seg_out);
void CompMixPoss(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int shiftDown);
void CompMixPoss3D(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown);
void MinErrorThresholding(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown, int *imgOut);
void Seg_GC_Full_3D_Blocks(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, int* Seg_out, int* imBlock);

#endif
