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

#ifndef SEEDSDETECTION_2D_H
#define SEEDSDETECTION_2D_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
//#include "itkApproximateSignedDistanceMapImageFilter.h"
#include <itkDanielssonDistanceMapImageFilter.h>

//typedef    float     InputPixelType;
//typedef itk::Image< InputPixelType,  2 >   InputImageType;

//int detect_seeds(itk::SmartPointer<InputImageType>, int , int , const double, float*);

//double get_maximum(double** A, int r1, int r2, int c1, int c2);

//void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, int* im_bin);

//int distMap(itk::SmartPointer<InputImageType> im, int r, int c, float* IMG);

int detectSeeds2D( float* IM, float* IM_out, unsigned short* IM_bin, int r, int c, double* sigma_min, double* sigma_max, double* scale, unsigned short* bImg, bool autoParamEstimation);

#endif