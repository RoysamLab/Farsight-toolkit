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
  Module:    $RCSfile: ftkPixelLevelAnalysisRules.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/27 1:00:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/

#ifndef __ftkPixelLevelAnalysisRules_cxx
#define __ftkPixelLevelAnalysisRules_cxx

#include "ftkPixelLevelAnalysis.h"

void ftk::PixelLevelAnalysis::SetInputs( std::string ROIImageNames, std::string TargetImageNames, std::string output_filenames, int radius, int mode ){

	pixelMode = mode;
	OutputFilename = output_filenames;

	pixel_distance = radius;

	typedef itk::ImageFileReader< UShortImageType >	FileReaderType;

	FileReaderType::Pointer ROIReader    = FileReaderType::New();
	FileReaderType::Pointer TargetReader = FileReaderType::New();

	ROIReader   ->SetFileName( ROIImageNames );
	TargetReader->SetFileName( TargetImageNames );

	ROIImagePtr   = UShortImageType::New();
	TargetImagePtr = UShortImageType::New();

	ROIReader   ->Update();
	TargetReader->Update();

	unsigned short uns_zero, uns_max;
	uns_zero = 0;
	uns_max  = itk::NumericTraits<unsigned short>::max();
/*
	typedef itk::RescaleIntensityImageFilter< UShortImageType, UShortImageType > RescaleUSUSType;
	RescaleUSUSType::Pointer RescaleUSUS1 = RescaleUSUSType::New();
	RescaleUSUSType::Pointer RescaleUSUS2 = RescaleUSUSType::New();
	RescaleUSUS1->SetOutputMaximum( uns_max  );
	RescaleUSUS2->SetOutputMaximum( uns_max  );
	RescaleUSUS1->SetOutputMinimum( uns_zero );
	RescaleUSUS2->SetOutputMinimum( uns_zero );
	RescaleUSUS1->SetInput( ROIReader   ->GetOutput() );
	RescaleUSUS2->SetInput( TargetReader->GetOutput() );
	RescaleUSUS1->Update();
	RescaleUSUS2->Update();
*/
	ROIImagePtr    = ROIReader->GetOutput(); //RescaleUSUS1
	TargetImagePtr = TargetReader->GetOutput(); //RescaleUSUS2

	if( ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[0] !=
		TargetImagePtr->GetLargestPossibleRegion().GetSize()[0] ||
		ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[1] !=
		TargetImagePtr->GetLargestPossibleRegion().GetSize()[1] ||
		ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[2] !=
		TargetImagePtr->GetLargestPossibleRegion().GetSize()[2] )
		std::cout << "WARNING: The input images must be of the same size!"
				  << std::endl;

	ROIImageName    = ROIImageNames;
	TargetImageName = TargetImageNames;
}


void ftk::PixelLevelAnalysis::WriteInitialOutputs(){

	std::string::size_type idx1 =    ROIImageName.find_last_of( '.' );
	std::string::size_type idx2 = TargetImageName.find_last_of( '.' );
	if ( idx1 == std::string::npos || idx2 == std::string::npos ){
    	std::cerr << "No file extension";
    	return;
	}
	ROIBinImageName    = ROIImageName.substr(0,idx1)    + "_bin.tif";
	TargetBinImageName = TargetImageName.substr(0,idx2) + "_bin.tif";

	std::ofstream output_txt_file( OutputFilename.c_str() , ios::app );
	output_txt_file	<< std::endl;

	output_txt_file	<< "ROI File          : " << ROIImageName	    << std::endl;
	output_txt_file	<< "ROI Binary File   : " << ROIBinImageName    << std::endl;
	output_txt_file	<< "Target File       : " << TargetImageName    << std::endl;
	output_txt_file	<< "Target Binary File: " << TargetBinImageName << std::endl;

	output_txt_file.close();
}

void ftk::PixelLevelAnalysis::WriteOutputImage(std::string OutName, UShortImageType::Pointer OutPtr){
	typedef itk::ImageFileWriter< UShortImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( OutName.c_str() );
	writer->SetInput( OutPtr );
	writer->Update();
}

bool ftk::PixelLevelAnalysis::RunAnalysis1(){
	this->WriteInitialOutputs();
	unsigned short thresh_roi, thresh_target;
	//this->WriteOutputImage( ROIBinImageName,    ROIImagePtr   );
	//this->WriteOutputImage( TargetBinImageName, TargetImagePtr);
	if( pixelMode == 1 ){
		thresh_roi    = returnthresh( ROIImagePtr,    1, 1 );
		thresh_target = returnthresh( TargetImagePtr, 1, 1 );
	}
	if( pixelMode == 4 ){
		thresh_roi    = returnthresh( ROIImagePtr,    2, 2 );
		thresh_target = returnthresh( TargetImagePtr, 2, 2 );
	}

	//Create an image with the atoms set as bright pixels
	UShortImageType::Pointer roi_bin    = UShortImageType::New();
	UShortImageType::Pointer target_bin = UShortImageType::New();
	UShortImageType::PointType origin1;
	origin1[0] = 0;
	origin1[1] = 0;
	origin1[2] = 0;
	roi_bin    ->SetOrigin(origin1);
	target_bin ->SetOrigin(origin1);

	UShortImageType::SizeType size1,size2; 
	size1[0] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[0];
	size1[1] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[1];
	size1[2] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[2];
	size2[0] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[0];
	size2[1] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[1];
	size2[2] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[2];

	UShortImageType::IndexType start;
	start[0] = 0; // first index on X
	start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z
	UShortImageType::RegionType region1,region2;
	region1.SetSize(size1);
	region1.SetIndex(start);
	region2.SetSize(size2);
	region2.SetIndex(start);

	roi_bin   ->SetRegions(region1);
	target_bin->SetRegions(region2);
	roi_bin   ->Allocate();
	target_bin->Allocate();
	roi_bin   ->FillBuffer(0);
	target_bin->FillBuffer(0);
	roi_bin   ->Update();
	target_bin->Update();

	double roi_count, target_count;
	roi_count = 0, target_count = 0;

	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< UShortImageType > IteratorType;
	IteratorType	  pix_buf_roi_bin( roi_bin,               roi_bin->GetRequestedRegion() );
	IteratorType	  pix_buf_tar_bin( target_bin,         target_bin->GetRequestedRegion() );
	ConstIteratorType pix_buf_roi_im ( ROIImagePtr,       ROIImagePtr->GetRequestedRegion() );
	ConstIteratorType pix_buf_tar_im ( TargetImagePtr, TargetImagePtr->GetRequestedRegion() );
	pix_buf_roi_bin.GoToBegin(); pix_buf_tar_bin.GoToBegin(); pix_buf_roi_im.GoToBegin(); pix_buf_tar_im.GoToBegin();
	unsigned short uns_zero, uns_max;
	uns_zero = 0;
	uns_max  = itk::NumericTraits<unsigned short>::max();
	for ( ; !(pix_buf_roi_bin.IsAtEnd() || pix_buf_tar_bin.IsAtEnd() || pix_buf_roi_im.IsAtEnd() || pix_buf_tar_im.IsAtEnd());
			++pix_buf_roi_bin, ++pix_buf_tar_bin, ++pix_buf_roi_im, ++pix_buf_tar_im ){
		if( pix_buf_roi_im.Get() >= thresh_roi ){
			roi_count = roi_count + 1.0;
			if( pix_buf_tar_im.Get() >= thresh_target ){
				target_count = target_count + 1.0;
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_zero);
			}
		} else {
			if( pix_buf_tar_im.Get() >= thresh_target ){
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_zero);
			}
		}
	}

	std::ofstream output_txt_file( OutputFilename.c_str(), ios::app );
	output_txt_file	<< std::endl;
	output_txt_file	<< "The percentage of pixels that are positive in the ROI image\nthat are also positive in the target image is:\n";

	long double percentage_of_pixels;
	percentage_of_pixels = 100.0 * (long double)target_count / (long double)roi_count;
	output_txt_file	<< setprecision (5);
	output_txt_file	<< percentage_of_pixels << std::endl;
	output_txt_file.close();

	WriteOutputImage(ROIBinImageName, roi_bin);
	WriteOutputImage(TargetBinImageName, target_bin);

	return true;
}

bool ftk::PixelLevelAnalysis::RunAnalysis2(){
	unsigned short uns_zero, uns_max;
	uns_zero = 0;
	uns_max  = itk::NumericTraits<unsigned short>::max();

	typedef itk::OtsuThresholdImageFilter	 < UShortImageType, UShortImageType > OtsuThresholdFilterType;
	OtsuThresholdFilterType::Pointer otsufilter = OtsuThresholdFilterType::New();
	otsufilter->SetInput( ROIImagePtr );
	otsufilter->SetOutsideValue(    0    );
	otsufilter->SetInsideValue ( uns_max );
	otsufilter->Update();
	UShortImageType::Pointer roi_bin1 = UShortImageType::New();
	roi_bin1 = otsufilter->GetOutput();

	this->WriteInitialOutputs();
	unsigned short thresh_target;
	if( pixelMode == 2 )
		thresh_target = returnthresh( TargetImagePtr, 1, 1 );
	if( pixelMode == 5 )
		thresh_target = returnthresh( TargetImagePtr, 2, 2 );

	typedef itk::BinaryBallStructuringElement<  unsigned short int, 3 >	StructuringElementType;
	typedef itk::BinaryDilateImageFilter	 < UShortImageType, UShortImageType, StructuringElementType > DilateFilterType;
	DilateFilterType::Pointer shell_image_filter = DilateFilterType::New();
	StructuringElementType shell_individual;
	shell_individual.SetRadius( pixel_distance ); //radius shell
	shell_individual.CreateStructuringElement();
	shell_image_filter->SetKernel( shell_individual );
	shell_image_filter->SetDilateValue( uns_max );
	shell_image_filter->SetInput( roi_bin1 );
	shell_image_filter->Update();
	UShortImageType::Pointer roi_bin2 = UShortImageType::New();
	roi_bin2 = shell_image_filter->GetOutput();

	//Create an image with the atoms set as bright pixels
	UShortImageType::Pointer roi_bin    = UShortImageType::New();
	UShortImageType::Pointer target_bin = UShortImageType::New();
	UShortImageType::PointType origin1;
	origin1[0] = 0;
	origin1[1] = 0;
	origin1[2] = 0;
	roi_bin    ->SetOrigin(origin1);
	target_bin ->SetOrigin(origin1);

	UShortImageType::SizeType size1,size2; 
	size1[0] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[0];
	size1[1] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[1];
	size1[2] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[2];
	size2[0] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[0];
	size2[1] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[1];
	size2[2] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[2];

	UShortImageType::IndexType start;
	start[0] = 0; // first index on X
	start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z
	UShortImageType::RegionType region1,region2;
	region1.SetSize(size1);
	region1.SetIndex(start);
	region2.SetSize(size2);
	region2.SetIndex(start);

	roi_bin   ->SetRegions(region1);
	target_bin->SetRegions(region2);
	roi_bin   ->Allocate();
	target_bin->Allocate();
	roi_bin   ->FillBuffer(0);
	target_bin->FillBuffer(0);
	roi_bin   ->Update();
	target_bin->Update();

	double roi_count, target_count;
	roi_count = 0, target_count = 0;

	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< UShortImageType > IteratorType;
	IteratorType	  pix_buf_roi_bin( roi_bin,        roi_bin->GetRequestedRegion() );
	IteratorType	  pix_buf_tar_bin( target_bin,     target_bin->GetRequestedRegion() );
	ConstIteratorType pix_buf_roi_im1( roi_bin1,       roi_bin1->GetRequestedRegion() );
	ConstIteratorType pix_buf_roi_im2( roi_bin2,       roi_bin2->GetRequestedRegion() );
	ConstIteratorType pix_buf_tar_im ( TargetImagePtr, TargetImagePtr->GetRequestedRegion() );
	pix_buf_roi_bin.GoToBegin(); pix_buf_tar_bin.GoToBegin(); pix_buf_roi_im1.GoToBegin(); pix_buf_roi_im2.GoToBegin(); pix_buf_tar_im.GoToBegin();
	for ( ; !(pix_buf_roi_bin.IsAtEnd() || pix_buf_tar_bin.IsAtEnd() || pix_buf_roi_im1.IsAtEnd() || pix_buf_roi_im2.IsAtEnd() || pix_buf_tar_im.IsAtEnd());
			++pix_buf_roi_bin, ++pix_buf_tar_bin, ++pix_buf_roi_im1, ++pix_buf_roi_im2, ++pix_buf_tar_im ){
		if( pix_buf_roi_im1.Get() != pix_buf_roi_im2.Get() ){
			roi_count = roi_count + 1.0;
			if( pix_buf_tar_im.Get() >= thresh_target ){
				target_count = target_count + 1.0;
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_zero);
			}
		} else {
			if( pix_buf_tar_im.Get() >= thresh_target ){
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_zero);
			}
		}
	}

	std::ofstream output_txt_file( OutputFilename.c_str(), ios::app );
	output_txt_file	<< std::endl;
	output_txt_file	<< "The percentage of pixels that are positive in the ROI image\nthat are also positive in the target image is:\n";

	long double percentage_of_pixels;
	percentage_of_pixels = 100.0 * (long double)target_count / (long double)roi_count;
	output_txt_file	<< setprecision (5);
	output_txt_file	<< percentage_of_pixels << std::endl;
	output_txt_file.close();

	WriteOutputImage(ROIBinImageName, roi_bin);
	WriteOutputImage(TargetBinImageName, target_bin);

	return true;
}

bool ftk::PixelLevelAnalysis::RunAnalysis3(){
	this->WriteInitialOutputs();
	unsigned short thresh_roi, thresh_target;
	//this->WriteOutputImage( ROIBinImageName,    ROIImagePtr   );
	//this->WriteOutputImage( TargetBinImageName, TargetImagePtr);
	if( pixelMode == 3 ){
		thresh_roi    = returnthresh( ROIImagePtr,    1, 1 );
		thresh_target = returnthresh( TargetImagePtr, 1, 1 );
	}
	if( pixelMode == 6 ){
		thresh_roi    = returnthresh( ROIImagePtr,    2, 2 );
		thresh_target = returnthresh( TargetImagePtr, 2, 2 );
	}

	//Create an image with the atoms set as bright pixels
	UShortImageType::Pointer roi_bin    = UShortImageType::New();
	UShortImageType::Pointer target_bin = UShortImageType::New();
	UShortImageType::PointType origin1;
	origin1[0] = 0;
	origin1[1] = 0;
	origin1[2] = 0;
	roi_bin    ->SetOrigin(origin1);
	target_bin ->SetOrigin(origin1);

	UShortImageType::SizeType size1,size2; 
	size1[0] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[0];
	size1[1] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[1];
	size1[2] = ROIImagePtr   ->GetLargestPossibleRegion().GetSize()[2];
	size2[0] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[0];
	size2[1] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[1];
	size2[2] = TargetImagePtr->GetLargestPossibleRegion().GetSize()[2];

	UShortImageType::IndexType start;
	start[0] = 0; // first index on X
	start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z
	UShortImageType::RegionType region1,region2;
	region1.SetSize(size1);
	region1.SetIndex(start);
	region2.SetSize(size2);
	region2.SetIndex(start);

	roi_bin   ->SetRegions(region1);
	target_bin->SetRegions(region2);
	roi_bin   ->Allocate();
	target_bin->Allocate();
	roi_bin   ->FillBuffer(0);
	target_bin->FillBuffer(0);
	roi_bin   ->Update();
	target_bin->Update();

	double roi_count, target_count, independent_target_count;
	roi_count = 0, target_count = 0, independent_target_count = 0;

	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< UShortImageType > IteratorType;
	IteratorType	  pix_buf_roi_bin( roi_bin,               roi_bin->GetRequestedRegion() );
	IteratorType	  pix_buf_tar_bin( target_bin,         target_bin->GetRequestedRegion() );
	ConstIteratorType pix_buf_roi_im ( ROIImagePtr,       ROIImagePtr->GetRequestedRegion() );
	ConstIteratorType pix_buf_tar_im ( TargetImagePtr, TargetImagePtr->GetRequestedRegion() );
	pix_buf_roi_bin.GoToBegin(); pix_buf_tar_bin.GoToBegin(); pix_buf_roi_im.GoToBegin(); pix_buf_tar_im.GoToBegin();
	unsigned short uns_zero, uns_max;
	uns_zero = 0;
	uns_max  = itk::NumericTraits<unsigned short>::max();
	for ( ; !(pix_buf_roi_bin.IsAtEnd() || pix_buf_tar_bin.IsAtEnd() || pix_buf_roi_im.IsAtEnd() || pix_buf_tar_im.IsAtEnd());
			++pix_buf_roi_bin, ++pix_buf_tar_bin, ++pix_buf_roi_im, ++pix_buf_tar_im ){
		if( pix_buf_roi_im.Get() >= thresh_roi ){
			++roi_count;
			if( pix_buf_tar_im.Get() >= thresh_target ){
				++target_count;
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_max );
				pix_buf_tar_bin.Set( uns_zero);
			}
		} else {
			if( pix_buf_tar_im.Get() >= thresh_target ){
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_max );
			}
			else{
				pix_buf_roi_bin.Set( uns_zero);
				pix_buf_tar_bin.Set( uns_zero);
			}
		}
		if( pix_buf_tar_im.Get() >= thresh_target )
			++independent_target_count;
	}

	std::ofstream output_txt_file( OutputFilename.c_str(), ios::app );
	output_txt_file	<< std::endl;
	long double percentage_of_pixels;
	output_txt_file	<< "The percentage of pixels that are positive in the ROI image is:\n";
	percentage_of_pixels = 100.0 * (long double)roi_count / ( (long double)size1[0] * (long double)size1[1] * (long double)size1[2] );
	output_txt_file	<< setprecision (5);
	output_txt_file	<< percentage_of_pixels << std::endl;

	output_txt_file	<< "The percentage of pixels that are positive in the Target image is:\n";
	percentage_of_pixels = 100.0 * (long double)independent_target_count / ( (long double)size1[0] * (long double)size1[1] * (long double)size1[2] );
	output_txt_file	<< setprecision (5);
	output_txt_file	<< percentage_of_pixels << std::endl;

	output_txt_file	<< "The percentage of pixels that are positive in the ROI image\nthat are also positive in the target image is:\n";
	percentage_of_pixels = 100.0 * (long double)target_count / (long double)roi_count;
	output_txt_file	<< setprecision (5);
	output_txt_file	<< percentage_of_pixels << std::endl;
	output_txt_file.close();

	WriteOutputImage(ROIBinImageName, roi_bin);
	WriteOutputImage(TargetBinImageName, target_bin);

	return true;
}


#endif
