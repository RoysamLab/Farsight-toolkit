#ifndef ADAPTIVE_BINARIZATION_H
#define ADAPTIVE_BINARIZATION_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<limits>
#include<exception>
#include <algorithm>

#include "new_graph.h"
#include "itkMinErrorThresholdImageFilter.h"


#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
////////////////////

typedef unsigned char inputPixelType;
typedef unsigned char UCHARPixelType;
typedef itk::Image< inputPixelType, 3 >         inputImageType;
typedef itk::Image< UCHARPixelType, 3 >         binaryImageType;
typedef itk::ImageFileReader<inputImageType> ReaderType;
typedef itk::ImageFileWriter<binaryImageType> WriterType;

double computePoissonProb(int intensity, double alpha);
binaryImageType::Pointer Adaptive_Binarization(inputImageType::Pointer _inputImage);
//void subtractGradientImage(unsigned char* IM, int r, int c, int z, int sampl_ratio);

#endif
