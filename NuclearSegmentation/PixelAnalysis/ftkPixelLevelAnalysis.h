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
  Module:    $RCSfile: ftkPixelLevelAnalysisRules.h,v $
  Language:  C++
  Date:      $Date: 2008/11/27 1:00:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/

#ifndef __ftkPixelLevelAnalysisRules_h
#define __ftkPixelLevelAnalysisRules_h

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>

#include "ftkImage.h"
#include "ftkCommon/ftkUtils.h"
#include "Nuclear_Association/ftkNuclearAssociationRules.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

namespace ftk
{

typedef itk::Image< unsigned short, 3 > UShortImageType;

class PixelLevelAnalysis{
	public:
		PixelLevelAnalysis(){
			OutputFilename.clear();
		}
		~PixelLevelAnalysis() { }

		void SetInputs( std::string ROIImageName, std::string TargetImageName,
			std::string output_filenames, int radius, int mode, int erode_radius );
		bool RunAnalysis1();
		bool RunAnalysis2();
		bool RunAnalysis3();

	private:
		std::string       ROIImageName;
		std::string    ROIBinImageName;
		std::string    TargetImageName;
		std::string TargetBinImageName;
		std::string  	OutputFilename;
		UShortImageType::Pointer   ROIImagePtr;
		UShortImageType::Pointer TargetImagePtr;
		void WriteInitialOutputs();
		void WriteOutputImage(std::string OutName, UShortImageType::Pointer OutPtr);
		void CleanUpSaltNPepperNThinStructs( UShortImageType::Pointer InputImage );
		unsigned short uns_zero, uns_max;
		int pixel_distance;
		int pixelMode;
		int erodeRadius;
};

}
#endif
