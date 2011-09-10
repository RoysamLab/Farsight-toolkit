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
#include "itkExtractImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

#include "itkLabelGeometryImageFilter.h"
#include "NuclearSegmentation/CytoplasmSegmentation/whole_cell.h"

#define MM_PI		3.14159265358979323846
#define MM_PI_2		1.57079632679489661923

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 3 > USImageType;
typedef float FloatPixelType;
typedef itk::Image< FloatPixelType, 3 > FloatImageType;

std::vector<float> compute_ec_features( USImageType::Pointer input_image,  USImageType::Pointer inp_labeled, int number_of_rois, unsigned short thresh, int surr_dist){

	//Dialate input first
	WholeCellSeg *dialate_filter = new WholeCellSeg;
	typedef itk::ExtractImageFilter< USImageType, UShortImageType > LabelExtractType;
	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	LabelExtractType::Pointer deFilter = LabelExtractType::New();
	USImageType::RegionType dRegion = inp_labeled->GetLargestPossibleRegion();
	dRegion.SetSize(2,0);
	deFilter->SetExtractionRegion(dRegion);
	deFilter->SetInput( inp_labeled );
	try{
		deFilter->Update();
	}
	catch( itk::ExceptionObject & excep ){
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	dialate_filter->set_nuc_img( deFilter->GetOutput() );
	dialate_filter->set_radius( surr_dist );
	dialate_filter->RunSegmentation();
	UShortImageType::Pointer input_lab = dialate_filter->getSegPointer();

	USImageType::Pointer input_labeled = USImageType::New();
	USImageType::PointType origint;
	origint[0] = 0;
	origint[1] = 0;
	origint[2] = 0;
	input_labeled->SetOrigin( origint );
	USImageType::IndexType startt;
	startt[0] = 0;  // first index on X
	startt[1] = 0;  // first index on Y
	startt[2] = 0;  // first index on Z
	USImageType::SizeType  sizet;
	sizet[0] = inp_labeled->GetLargestPossibleRegion().GetSize()[0];  // size along X
	sizet[1] = inp_labeled->GetLargestPossibleRegion().GetSize()[1];  // size along Y
	sizet[2] = inp_labeled->GetLargestPossibleRegion().GetSize()[2];  // size along Z
	USImageType::RegionType regiont;
	regiont.SetSize( sizet );
	regiont.SetIndex( startt );
	input_labeled->SetRegions( regiont );
	input_labeled->Allocate();
	input_labeled->FillBuffer(0);
	input_labeled->Update();

	ConstIteratorType pix_buf1( input_lab, input_lab->GetRequestedRegion() );
	IteratorType iterator2 ( input_labeled, input_labeled->GetRequestedRegion() );
	iterator2.GoToBegin();
	for ( pix_buf1.GoToBegin(); !pix_buf1.IsAtEnd(); ++pix_buf1 ){
		iterator2.Set( pix_buf1.Get() );
		++iterator2;
	}

	std::vector< float > quantified_numbers;
	std::vector< unsigned short > labelsList;

	typedef itk::LabelGeometryImageFilter< USImageType > GeometryFilterType;
	typedef GeometryFilterType::LabelIndicesType labelindicestype;

	//int size1 = input_image->GetLargestPossibleRegion().GetSize()[0];
	//int size2 = input_image->GetLargestPossibleRegion().GetSize()[1];

	GeometryFilterType::Pointer geomfilt1 = GeometryFilterType::New();

	geomfilt1->SetInput( input_labeled );
	geomfilt1->SetCalculatePixelIndices( true );
	geomfilt1->Update();
	labelsList = geomfilt1->GetLabels();

	bool zp=false;
	for( unsigned short i=0; (int)i < labelsList.size(); ++i ){
		if( labelsList[i] == 0 ){ zp=true; continue; }
		std::vector<float> quantified_numbers_cell;
		//assert(quantified_numbers_cell.size() == number_of_rois);
		for( int j=0; j<number_of_rois; ++j ) quantified_numbers_cell.push_back((float)0.0);
		double centroid_x = (double)(geomfilt1->GetCentroid(labelsList[i])[0]);
		double centroid_y = (double)(geomfilt1->GetCentroid(labelsList[i])[1]);
		labelindicestype indices1;
		indices1 = geomfilt1->GetPixelIndices(labelsList[i]);
		for( labelindicestype::iterator itPixind = indices1.begin(); itPixind!=indices1.end(); ++itPixind ){
			IteratorType iterator1 ( input_image, input_image->GetRequestedRegion() );
			iterator1.SetIndex( *itPixind );
			if( iterator1.Get() < thresh )
				continue;
			double x = (double)(iterator1.GetIndex()[0]);
			double y = (double)(iterator1.GetIndex()[1]);
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
		for( int j=0; j<number_of_rois; ++j ) quantified_numbers.push_back(quantified_numbers_cell[j]);
	}
	std::vector< float > qfied_num;
	int qnum_sz = zp? (int)(labelsList.size()-1) : (int)(labelsList.size());
	for( int i=0; i<qnum_sz; ++i ){
		int counter=0;
		for( int j=0; j<number_of_rois; ++j ){
			if( quantified_numbers[(i*number_of_rois+j)] > 1 )
				++counter;
		}
		qfied_num.push_back((float)counter);
	}
	return qfied_num;
}

unsigned short returnthresh( USImageType::Pointer input_image, int num_bin_levs, int num_in_fg ){
	//Instantiate the different image and filter types that will be used
	typedef itk::ImageRegionConstIterator< USImageType > ConstIteratorType;
	//typedef itk::Statistics::ScalarImageToHistogramGenerator< USImageType > ScalarImageToHistogramGeneratorType;
	//typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::Statistics::Histogram< FloatPixelType > HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	std::cout<<"Starting threshold computation\n";

	//Create a temporary histogram container:
	const int numBins = 256;
	double tempHist[numBins];
	for(int i=0; i<numBins; ++i)
	{
		tempHist[i] = 0;
	}

	//Populate the histogram (assume pixel type is actually uchar):
	ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		USPixelType pix = it.Get();
		if(pix <= 255)
		{
			++tempHist[pix];
		}
	}
	
	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<unsigned short>::max();
	double max = 0.0;
	for(int i=0; i<numBins; ++i)
	{
		if( tempHist[i] > max )
			max = tempHist[i];
	}

	double scaleFactor = 1;
	if(max >= floatIntegerMax)
	{
		scaleFactor = floatIntegerMax / max;
	}

	HistogramType::Pointer histogram = HistogramType::New() ;
	// initialize histogram
	HistogramType::SizeType size;
	HistogramType::MeasurementVectorType lowerBound;
	HistogramType::MeasurementVectorType upperBound;

	lowerBound.SetSize(1);
	upperBound.SetSize(1);
	size.SetSize(1);

	lowerBound.Fill(0.0);
	upperBound.Fill(255.0);
	size.Fill(numBins);

	histogram->SetMeasurementVectorSize(1);
	histogram->Initialize(size, lowerBound, upperBound ) ;

	int i=0;
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter )
	{
		float norm_freq = (float)(tempHist[i] * scaleFactor);
		iter.SetFrequency(norm_freq);
		++i;
	}

	std::cout<<"Histogram computed\n";

	//ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
	//scalarImageToHistogramGenerator->SetNumberOfBins( 256 );
	//scalarImageToHistogramGenerator->SetInput( input_image);
	//scalarImageToHistogramGenerator->Compute();
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bin_levs );
	//calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
	calculator->SetInputHistogram( histogram );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	FloatPixelType thresh;

	for(int i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<FloatPixelType>(*itNum));

	std::cout<<"Threshold computed: "<<thresh<<std::endl;

	return (USPixelType)(thresh+0.5);

}

#endif
