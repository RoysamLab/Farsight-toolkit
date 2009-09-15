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


//added by Yousef on 9/3/2009
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageFileWriter.h"
////////////////////

int Cell_Binarization_2D(unsigned char* imgIn, unsigned short *imgOut, int R, int C, int shd); 
int Cell_Binarization_3D(unsigned char *imgIn, unsigned short* imgOut, int R, int C, int Z, int shd, int div); 
double compute_poisson_prob(double intensity, double alpha);
void Seg_GC_Full_2D(unsigned char* IM, int r, int c, double alpha_F, double alpha_B, double P_I, int* num_nodes, int* num_edges, unsigned short* Seg_out);
void Seg_GC_Full_3D(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, unsigned short* Seg_out);
void CompMixPoss(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int shiftDown);
void CompMixPoss3D(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown);
void MinErrorThresholding(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown, unsigned short *imgOut);
void Seg_GC_Full_3D_Blocks(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, unsigned short* Seg_out, int* imBlock);
//void subtractGradientImage(unsigned char* IM, int r, int c, int z, int sampl_ratio);
#endif
