/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef CELL_BINARIZATION_H
#define CELL_BINARIZATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <exception>
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
int Neuron_Binarization_3D(unsigned char *imgIn, unsigned short* imgOut, int R, int C, int Z, int shd, int div);
double compute_poisson_prob(double intensity, double alpha);
void Seg_GC_Full_2D(unsigned char* IM, int r, int c, double alpha_F, double alpha_B, double P_I, int* num_nodes, int* num_edges, unsigned short* Seg_out);
void Seg_GC_Full_3D(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, unsigned short* Seg_out);
void CompMixPoss(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int shiftDown);
void CompMixPoss3D(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown);
void MinErrorThresholding(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown, unsigned short *imgOut);
void Seg_GC_Full_3D_Blocks(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, unsigned short* Seg_out, int* imBlock);
void threeLevelMinErrorThresh(unsigned char* im, float* Alpha1, float* Alpha2, float* Alpha3, float* P_I1, float* P_I2, int r, int c, int z);
void Seg_GC_Full_2D_Three_Level(unsigned char* IM,int r, int c,double alpha_F,double alpha_B1, double alpha_B2, double P_I1, double P_I2,int* num_nodes, int* num_edges, unsigned short* Seg_out);
#endif
