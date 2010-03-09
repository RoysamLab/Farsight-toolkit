#ifndef _COLORSEGMENTATION_H_
#define _COLORSEGMENTATION_H_

//ITK Includes
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
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

#include "dhEvalState.h"
#include "RLI_Slice.h"
#include "dhColors.h"
#include "dhHistogram.h"
#include "RGB_Atype.h"
#include "dhSeedGrid.h"
#include "dhClassifiers.h"

typedef unsigned char UcharPixelType;
typedef unsigned short UshrtPixelType;
typedef float FloatPixelType;
typedef itk::RGBPixel< UcharPixelType > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 > RGBImageType;
typedef itk::Image< UcharPixelType, 3 > UcharImageType;
typedef itk::Image< UshrtPixelType, 3 > UshrtImageType;
typedef itk::Image< FloatPixelType, 3 > FloatImageType;

enum Pixel_Class { UNKNOWN, RED_CELL, BLUE_CELL, BKGD_FIELD };

class ColorSegmentation
{
	//Flags
	int foreground_dark, number_of_bins, number_in_foreground, bin_done;

	//Archetypal Colors
	//RGB arch_typ1, arch_typ2, bkgrnd_typ; //1-> Red-ish 2-> Blue-ish

	//Image Pointers
	RGBImageType::Pointer rgb_input;
	UcharImageType::Pointer intial_binarization;
	UcharImageType::Pointer red_image;
	UcharImageType::Pointer lime_image;
	UcharImageType::Pointer intensity_image;
	UcharImageType::Pointer red_wts;
	UcharImageType::Pointer blue_wts;
public:
	//Constructor
	ColorSegmentation(RGBImageType::Pointer input, int fore_ground_dark,int number_bins, int number_in_fg);

	//Destructor
	~ColorSegmentation();

	//Run Initial Binarization
	void RunInitialBinarization();
	//Get Binary
	UcharImageType::Pointer get_binary();

	dh::RGBHistogram * hist;
	dh::RGB_Atype at;
	dh::_RGB arch_typ1, arch_typ2, bkgrnd_typ; //1-> Red-ish 2-> Blue-ish

	//Compute Archetypal Colors
	void FindArchetypalColors();

	//Set Previously computed Atypes
	//void SetAtypes( RGB a_typ1, RGB a_typ2, RGB bk_typ ){
	//	arch_typ1 = a_typ1; arch_typ2 = a_typ2; bkgrnd_typ = bk_typ;
	//}

	//Get Grayscales Based On Distances From Atypes
	void ComputeClassWeights();

	//A few parameters:
	bool RLI_MODE;
	bool GEN_PROJ;

private:
	void go_best_dir( dh::_RGB& ma, bool& moved, const dh::_RGB& sa, const int r = 1, int res = 1 );
};

#endif
