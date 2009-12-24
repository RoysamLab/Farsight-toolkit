/*
 *  ftkPreprocess.h
 *  Farsight
 *
 *  Created by RAGHAV on 12/3/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */



#ifndef __ftkPreprocess_h
#define __ftkPreprocess_h


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "ftkImage.h"
#include "itkImage.h"
#include <vector>
//ITK Preprocessing includes

#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkOpeningByReconstructionImageFilter.h"
#include "itkClosingByReconstructionImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkLaplacianImageFilter.h"
#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkMinMaxCurvatureFlowImageFilter.h"



class ftkPreprocess
{
public:
	ftkPreprocess();
	
	typedef unsigned char   InpPixelType;
	typedef float			FloatPixelType;
	typedef double			DoublePixelType;
	
	typedef itk::Image<InpPixelType, 3> InpImageType;		
    typedef itk::Image<FloatPixelType, 3> FloatImageType;
	typedef itk::Image<DoublePixelType, 3> DoubleImageType;	
	
	ftk::Image::Pointer CurvAnisotropicDiffusion(void);
	ftk::Image::Pointer GrayscaleClose(void);
	ftk::Image::Pointer GrayscaleOpen(void);
	ftk::Image::Pointer GrayscaleDilate(void);
	ftk::Image::Pointer GrayscaleErode(void);
	ftk::Image::Pointer SigmoidFilter(void);
	ftk::Image::Pointer MedianFilter(void);
	ftk::Image::Pointer GADiffusion(void);
	ftk::Image::Pointer Resample(void);
	ftk::Image::Pointer OpeningbyReconstruction(void);
    ftk::Image::Pointer ClosingbyReconstruction(void);	
	ftk::Image::Pointer MeanFilter(void);
	ftk::Image::Pointer LaplacianFilter(void);
	ftk::Image::Pointer ThreeDSmoothingRGFilter(void);
	ftk::Image::Pointer NormalizeImage(void);
	ftk::Image::Pointer ShiftScale(void);
	ftk::Image::Pointer SobelEdgeDetection(void);
	ftk::Image::Pointer CurvatureFlow(void);
	ftk::Image::Pointer MinMaxCurvatureFlow(void);
	ftk::Image::Pointer VesselFilter(void);// 
	
	ftk::Image::Pointer myImg;
	std::vector <double> filterParams;
	//private:

};



#endif
