#ifndef _COLORSEGMENTATION_H_
#define _COLORSEGMENTATION_H_

//ITK Includes
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkMeanImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//STL and Other Includes
#include <iostream>
#include <algorithm>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "dhColors.h"
#include "dhHistogram.h"
#include "dhSeedGrid.h"
#include "dhClassifiers.h"

typedef float FloatPixelType;

typedef itk::RGBPixel< dh::RGBType > RGBPixelType;
typedef itk::RGBPixel< dh::RLIType > RLIPixelType;

typedef itk::Image< RGBPixelType, 3 > RGBImageType;
typedef itk::Image< RLIPixelType, 3 > RLIImageType;
typedef itk::Image< unsigned char, 3> UcharImageType;

enum Pixel_Class { UNKNOWN, RED_CELL, BLUE_CELL, BKGD_FIELD };

class ColorSegmentation
{
public:
	//Constructor
	ColorSegmentation(RGBImageType::Pointer input);
	~ColorSegmentation(){if(hist) delete hist;};//Destructor

	void SetTesting(bool t = true){ TESTING = t; };		//default is false
	void SetIgnoreBackground(bool i = true){ IGNORE_BACKGROUND = i; }; //default is false
	void SetLightBackground(bool d = true){ LIGHT_BACKGROUND = d; }; //default is false
	void SetGenerateProjections(bool d = true){ GEN_PROJ = d; }; //default is false

	//Methods:
	void InvertRGBImage(RGBImageType::Pointer img);
	void SmoothRGBImage(RGBImageType::Pointer img);
	void MaskBackgroundFromInput();

	void TransformToRLI();			//First step
	void ComputeBinary(int num_bins, int num_in_fg);
	void FindArchetypalColors();	//Compute Archetypal Colors
	void SetArchetypalColors(dh::RLI r, dh::RLI b, dh::RLI w);
	void SetArchetypalColors(dh::_RGB r, dh::_RGB b, dh::_RGB w);
	void ComputeClassWeights2();
	void ComputeClassWeights();		//Get Grayscales Based On Distances From Atypes
	void VoteBasedOnWeights();

	//A few parameters:
	bool IGNORE_BACKGROUND;
	bool LIGHT_BACKGROUND;
	bool TESTING;

	bool GEN_PROJ; //Generate Projection Images:

protected:
	//Image Pointers
	RGBImageType::Pointer rgb_input;
	RLIImageType::Pointer rli_image;
	UcharImageType::Pointer bin_image;

	dh::Histogram * hist;

	UcharImageType::Pointer red_weights;
	UcharImageType::Pointer blue_weights;

	//Intermediate values
	// These actually contain RLI values:
	dh::RLI archTypRED, archTypBLUE, archTypBACK; //1->Red-ish 2->Blue-ish

private:
	void GenerateProjection(int dir, std::string outFilename);
	void GenerateATColors(std::string outFilename);
	void GenerateColors(dh::_RGB c1, dh::_RGB c2, dh::_RGB c3, std::string outFilename);
};

#endif
