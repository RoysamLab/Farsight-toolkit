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
/**************************************************************************  
 *  Volume dataset processing:
 *  Author: Xiaosong Yuan, RPI
 *  Modified on Sep. 29, 2005  
 *  Adapted Jan. 2010 by Isaac Abbott 
 *
 *  Input parameters
 *  1. sizeExpand   
 *  2. preproess          
 *************************************************************************/
#ifndef __mdlVolumeProcess_h
#define __mdlVolumeProcess_h

#include <string>
#include <vector>
#include <iostream>

#include "itkImage.h"
//#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkOtsuThresholdImageFilter.h" 
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkShiftScaleImageFilter.h"
//#include "itkNormalizeImageFilter.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

namespace mdl
{

class VolumeProcess
{
public:
	typedef unsigned char PixelType;
    static const unsigned int Dimension = 3;
    typedef itk::Image< PixelType, Dimension > ImageType;

	VolumeProcess();
	void SetUseCAD(bool inp = true){ useCAD = inp; };
	void SetOverwrite(bool inp = true){ overwriteInput = inp; };
	void SetDebug(bool inp = true){ debug = inp; };
	void SetInput(ImageType::Pointer inImage);
	bool Update();
	ImageType::Pointer GetOutput();

private:
	static const unsigned char m_NumberOfHistogramBins = 128;

	//Parameters
	bool useCAD;			//Use curvature anisotropic diffusion
	bool overwriteInput;	//Replace the original image with processed one
	bool debug;				//If debug is true, process in steps and print stuff
	
	//Images
	ImageType::Pointer m_inputImage;
	ImageType::Pointer m_outputImage;

	//Functions:
	double getItkOtsuThreshold(ImageType::Pointer img);
	double getXiaoLiangOtsuThreshold(ImageType::Pointer img);
	void dialateImage(ImageType::Pointer img, int iterations, int border);

};

}  // end namespace mdl

#endif
