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
#ifndef __ftkMultipleImageHandler_h
#define __ftkMultipleImageHandler_h

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

// NOTE THAT THIS CLASS IS WRITTEN FOR 8-BIT IMAGE PROCESSING
#include "itkImage.h"
#include <itkImageIOBase.h>
#include <itkLightObject.h>
#include <itkObjectFactory.h>
#include <itkSmartPointer.h>
#include "itkRGBPixel.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"

#include <vector>
#include <string>

#include <math.h>

#include "ftkUtils.h"
#include "ftkProjectManager.h"

namespace ftk
{

class MultipleImageHandler
{
	typedef unsigned char   UCharPixelType;
	typedef unsigned short	UShortPixelType;
	typedef float			FloatPixelType;
	typedef double			DoublePixelType;
	typedef itk::RGBPixel<UCharPixelType> RGBPixelType;
	
	typedef itk::Image<UCharPixelType, 2> UCharImageType2D;
	typedef itk::Image<UCharPixelType, 3> UCharImageType3D;
	typedef itk::Image<UShortPixelType, 2> UShortImageType2D;
	typedef itk::Image<UShortPixelType, 3> UShortImageType3D;
	typedef itk::Image<RGBPixelType, 2> RGBImageType2D;
	typedef itk::Image<RGBPixelType, 3> RGBImageType3D;

	typedef std::vector<std::string> StrVector;
	typedef std::pair<int, int> PairType;
	typedef std::vector< PairType > PairVector;

public:
	MultipleImageHandler();
	void SetOutputDirectory(std::string dir){ outputDirectory = dir; };
	void SetOutputBase(std::string base){ outputBaseString = base; };

	//These two functions take in image series and re-partion it into 3D blocks.
	void SeriesToBlocks(std::string seriesFormat, int startIndex, int endIndex, int dx, int dy, int dz, int color=-1);
	void SeriesToBlocks(StrVector inFiles, int dx, int dy, int dz, int color=-1);

	//Take in an image series and create a maximum projection image along the z dimension
	UCharImageType2D::Pointer SeriesProjection(std::string seriesFormat, int startIndex, int endIndex, std::string outName = "");
	UCharImageType2D::Pointer SeriesProjection(StrVector inFiles, std::string outName = "");
	UCharImageType2D::Pointer ImageProjection(std::string inFile, std::string outName = "");

	//Use these functions to just load a specific region in a 3D image:
	UCharImageType3D::Pointer ExtractRegion(StrVector inFiles, UCharImageType3D::RegionType region, bool rescale = false, std::string fname = "");
	UCharImageType3D::Pointer ExtractRegionColor(StrVector inFiles, UCharImageType3D::RegionType region, int color=-1, std::string fname = "");

	//These functions help to divide up the image into a region:
	UCharImageType3D::RegionType CreateRegion(PairType xPair, PairType yPair, PairType zPair);
	PairVector CreateOutputRegions(int size, int divs);
	
private:
	//functions:
	UCharImageType2D::Pointer ReadImage2D(std::string filename);

	std::string outputDirectory;
	std::string outputBaseString;

};

} // end namespace

#endif
