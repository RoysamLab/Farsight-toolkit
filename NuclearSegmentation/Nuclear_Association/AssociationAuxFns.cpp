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
#ifndef _ASC_FEAT_AUX_FN_CPP_
#define _ASC_FEAT_AUX_FN_CPP_

#include <math.h>

#include "itkImage.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"

#include "ftkCommon/itkLabelGeometryImageFilter.h"
#include "NuclearSegmentation/CytoplasmSegmentation/whole_cell.h"

#define MM_PI		3.14159265358979323846
#define MM_PI_2		1.57079632679489661923

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 3 > USImageType;

std::vector<float> compute_ec_features( USImageType::Pointer input_image,  USImageType::Pointer inp_labeled, int number_of_rois, unsigned short thresh){

	//Dialate input first
	WholeCellSeg *dialate_filter = new WholeCellSeg;
	typedef itk::ExtractImageFilter< USImageType, UShortImageType > LabelExtractType;
	typedef itk::ExtractImageFilter< UShortImageType, USImageType > LabelExtractType1;
	LabelExtractType::Pointer deFilter = LabelExtractType::New();
	USImageType::RegionType dRegion = inp_labeled->GetLargestPossibleRegion();
	dRegion.SetSize(2,0);
	deFilter->SetExtractionRegion(dRegion);
	deFilter->SetInput( inp_labeled );
	deFilter->Update();
	try{
		deFilter->Update();
	}
	catch( itk::ExceptionObject & excep ){
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	dialate_filter->set_nuc_img( deFilter->GetOutput() );
	dialate_filter->RunSegmentation();
	UShortImageType::Pointer input_lab = dialate_filter->getSegPointer();
	LabelExtractType1::Pointer deFilter1 = LabelExtractType1::New();
	UShortImageType::RegionType dRegion1 = input_lab->GetLargestPossibleRegion();
	deFilter1->SetExtractionRegion(dRegion1);
	deFilter1->SetInput( input_lab );
	deFilter1->Update();
	try{
		deFilter1->Update();
	}
	catch( itk::ExceptionObject & excep ){
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	delete dialate_filter;

	USImageType::Pointer input_labeled = deFilter1->GetOutput();

	std::vector< float > quantified_numbers;
	std::vector< unsigned short > labelsList;

	typedef itk::LabelGeometryImageFilter< USImageType > GeometryFilterType;
	typedef itk::LabelStatisticsImageFilter< USImageType,USImageType > StatisticsFilterType;
	typedef itk::CastImageFilter< USImageType, USImageType > CastUSUSType;
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	typedef GeometryFilterType::LabelIndicesType labelindicestype;

	//int size1 = input_image->GetLargestPossibleRegion().GetSize()[0];
	//int size2 = input_image->GetLargestPossibleRegion().GetSize()[1];

	GeometryFilterType::Pointer geomfilt1 = GeometryFilterType::New();

	StatisticsFilterType::Pointer statsfilt = StatisticsFilterType::New();
	statsfilt->UseHistogramsOn();

	geomfilt1->SetInput( input_labeled );
	geomfilt1->SetCalculatePixelIndices( true );
	geomfilt1->Update();

	statsfilt->SetInput( input_image );
	statsfilt->SetLabelInput( input_labeled );
	statsfilt->Update();
	labelsList = geomfilt1->GetLabels();

	CastUSUSType::Pointer castUSUSfilter = CastUSUSType::New();

	for( unsigned short i=0; (int)i <= labelsList.size(); ++i ){
		std::vector<float> quantified_numbers_cell;
		//assert(quantified_numbers_cell.size() == number_of_rois);
		for( int j=0; j<number_of_rois; ++j ) quantified_numbers_cell.push_back(0);
		double centroid_x = (double)(geomfilt1->GetCentroid(i)[0]);
		double centroid_y = (double)(geomfilt1->GetCentroid(i)[1]);
		if( !statsfilt->GetSum(labelsList[i]) ){
			IteratorType iterator ( input_labeled, input_labeled->GetRequestedRegion() );
			labelindicestype indices1;
			indices1 = geomfilt1->GetPixelIndices(labelsList[i]);
			labelindicestype::iterator itPixind = indices1.begin();
			for( int j=0; j<(int)indices1.size(); ++j, ++itPixind ){
				IteratorType iterator1 ( input_image, input_image->GetRequestedRegion() );
				iterator1.SetIndex( *itPixind );
				if( thresh >= iterator1.Get() )
					continue;
				iterator.SetIndex( *itPixind );
				double x = (double)(iterator.GetIndex()[0]);
				double y = (double)(iterator.GetIndex()[1]);
				double angle = atan2((centroid_y-y),fabs(centroid_x-x));
				if( (centroid_x-x)>0 )
					angle += MM_PI_2;
				else
					angle = MM_PI+MM_PI-(angle+MM_PI_2);
				angle = ((number_of_rois-1)*angle)/(2*MM_PI);
				double angle_fraction[1];
				if( modf( angle, angle_fraction ) > 0.5 )
					angle = ceil( angle );
				else
					angle = floor( angle );
				quantified_numbers_cell[(int)angle] += iterator1.Get();
			}
		}
		for( int j=0; j<number_of_rois; ++j ) quantified_numbers.push_back(quantified_numbers_cell[j]);
	}
	std::vector< float > qfied_num;
	for( int i=0; i<(int)(quantified_numbers.size()/number_of_rois); ++i ){
		int counter=0;
		for( int j=0; j<number_of_rois; ++j ){
			if( quantified_numbers[i*number_of_rois+j] )
				++counter;
		}
		qfied_num.push_back((float)counter);
	}
	return qfied_num;
}

unsigned short returnthresh( USImageType::Pointer input_image, int num_bin_levs, int num_in_fg ){
//Instantiate the different image and filter types that will be used
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	typedef itk::ImageRegionConstIterator< USImageType > ConstIteratorType;
	typedef itk::Statistics::ScalarImageToHistogramGenerator< USImageType > ScalarImageToHistogramGeneratorType;
	typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
	CalculatorType::Pointer calculator = CalculatorType::New();
	scalarImageToHistogramGenerator->SetNumberOfBins( 128 );
	calculator->SetNumberOfThresholds( num_bin_levs );
	scalarImageToHistogramGenerator->SetInput( input_image);
	scalarImageToHistogramGenerator->Compute();
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	USPixelType thresh;

	for(int i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<USPixelType>(*itNum));

	return thresh;

}

#endif