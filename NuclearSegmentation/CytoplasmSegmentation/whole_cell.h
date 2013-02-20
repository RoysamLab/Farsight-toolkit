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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef _CYTOPLASM_SEG_
#define _CYTOPLASM_SEG_

//ITK Includes
#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryThresholdImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
//#include "itkImageFileWriter.h" 

//LOCAL INCLUDES
#include "NuclearSegmentation/yousef_core/cell_binarization/cell_binarization.h"
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"

//Standard C-Library Includes
#include <iostream>
#include <limits.h>
#include <float.h>

//int Cell_Binarization_2D ( unsigned char* imgIn, int *imgOut, int R, int C, int shd );
int gradient_enhanced_distance_map( float *INP_IM_2D, float *GRAD_IMW, int size1, int size2 );

typedef int IntPixelType;
typedef itk::Image< unsigned char, 2 > UcharImageType;
typedef itk::Image< unsigned short int, 2 > UShortImageType;
typedef itk::Image< int, 2 > IntImageType;
typedef itk::Image< float, 2 > FltImageType;
typedef itk::Image< unsigned long int, 2 > ULongImageType;

class WholeCellSeg{
	unsigned short *bin_Image;
	//Binarization Parameters
	int shift_bin, num_levels, num_levels_incl, draw_real_bounds, draw_synth_bounds;
	//Scaling and mem_scaling
	int radius_of_synth_bounds, remove_small_objs, scaling, mem_scaling,use_mem_img;

	//Flags
	int nuc_im_set, cyt_im_set, mem_im_set, bin_done, seg_done;

	//Image Pointers
	UShortImageType::Pointer nuclab_inp;
	UShortImageType::Pointer nuclab_inp_cpy;
	UShortImageType::Pointer cyt_im_inp;
	UShortImageType::Pointer mem_im_inp;
	UShortImageType::Pointer bin_im_out;
	UShortImageType::Pointer intermediate_bin_im_out;
	UShortImageType::Pointer seg_im_out;

public:
	//Function to set parameters
	void set_parameters ( int *parameters );
	void set_radius ( int parameter );

	//Set the input images, they have to be 16-bit grayscale images passed as ITK smart pointers
	void set_nuc_img( UShortImageType::Pointer nuc_im_16 ){
		nuclab_inp = nuc_im_16;
		nuc_im_set = 1;
	}
	void set_cyt_img( UShortImageType::Pointer cyt_im_16 ){
		cyt_im_inp = cyt_im_16;
		cyt_im_set = 1;
		draw_real_bounds=1; draw_synth_bounds=0;

	}
	void set_mem_img( UShortImageType::Pointer mem_im_16 ){
		mem_im_inp = mem_im_16;
		mem_im_set = 1;
		use_mem_img = 1;
	}

	void RealBoundaries();
	void BinarizationForRealBounds();
	void RemoveSmallObjs();
	void SyntheticBoundaries();

	//Two steps in the segmentation
	//1-> Run Binarization
	void RunBinarization();

	//2-> Create distance map and call marker based watershed
	void RunSegmentation();

	//Return pointer to binary after binarization
	UShortImageType::Pointer getBinPointer();

	//Return pointer to final segmenation result
	UShortImageType::Pointer getSegPointer();

	//Constructor
	WholeCellSeg();

	//Destructor
	~WholeCellSeg();
};

#endif
