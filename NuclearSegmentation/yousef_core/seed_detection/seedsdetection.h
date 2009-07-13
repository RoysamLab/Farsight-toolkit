#ifndef SEEDSDETECTION_2D_H
#define SEEDSDETECTION_2D_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include "itkImage.h"
#include "itkSimpleFilterWatcher.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
//#include "itkApproximateSignedDistanceMapImageFilter.h"
#include <itkDanielssonDistanceMapImageFilter.h>

typedef    float     InputPixelType;
typedef itk::Image< InputPixelType,  2 >   InputImageType;

using namespace std;

int detect_seeds(itk::SmartPointer<InputImageType>, int , int , const double, float*);

double get_maximum(double** A, int r1, int r2, int c1, int c2);

void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, int* im_bin);

int distMap(itk::SmartPointer<InputImageType> im, int r, int c, float* IMG);

int detectSeeds2D( float* IM, float* IM_out, int* IM_bin, int r, int c, double sigma_min, double sigma_max, double scale, int* bImg);

#endif