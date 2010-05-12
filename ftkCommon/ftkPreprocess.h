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
#ifndef __ftkPreprocess_h
#define __ftkPreprocess_h


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

// NOTE THAT THIS CLASS IS WRITTEN FOR 8-BIT IMAGE PROCESSING

#include "ftkImage.h"
#include "itkImage.h"
#include "ftkObject.h"
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
#include "itkInvertIntensityImageFilter.h"
#include "itkPolylineMask2DImageFilter.h"
#include "itkPolyLineParametricPath.h"
#include "itkExtractImageFilter.h"


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
	
	void MaskTest(void);

	void MaskImage(std::vector< ftk::Object::Point > roiPoints); //skip
	void InvertIntensity(void); //done
	void CurvAnisotropicDiffusion(void); //done
	void GrayscaleClose(void); //done
	void GrayscaleOpen(void);  //done
	void GrayscaleDilate(void);
	void GrayscaleErode(void);
	void SigmoidFilter(void); 
	void MedianFilter(void);   //done
	void GADiffusion(void);		//done
	void Resample(void);		//skip
	void OpeningbyReconstruction(void); //skip
    void ClosingbyReconstruction(void);	//skip
	void MeanFilter(void);
	void LaplacianFilter(void);
	void ThreeDSmoothingRGFilter(void);
	void NormalizeImage(void);
	void ShiftScale(void);
	void SobelEdgeDetection(void);
	void CurvatureFlow(void);
	void MinMaxCurvatureFlow(void);
	void VesselFilter(void);// 
	
	ftk::Image::Pointer myImg;
	int channelNumber;
	std::vector <double> filterParams;
	//private:

};



#endif
