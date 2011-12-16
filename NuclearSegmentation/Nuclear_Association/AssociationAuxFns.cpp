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
#include "itkSignedDanielssonDistanceMapImageFilter.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_real.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/vnl_double_3x3.h"

#include "itkLabelGeometryImageFilter.h"
#include "NuclearSegmentation/CytoplasmSegmentation/whole_cell.h"

#define MM_PI		3.14159265358979323846
#define MM_PI_2		1.57079632679489661923

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 3 > USImageType;
typedef float FloatPixelType;
typedef itk::Image< FloatPixelType, 3 > FloatImageType;

std::vector<float> compute_ec_features( USImageType::Pointer input_image,  USImageType::Pointer inp_labeled, int number_of_rois, unsigned short thresh, int surr_dist){

	std::vector< float > qfied_num;
	std::vector< unsigned short > labelsList;
	std::vector< double > quantified_numbers_cell;

	typedef itk::ExtractImageFilter< USImageType, UShortImageType > LabelExtractType;
	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	typedef itk::SignedDanielssonDistanceMapImageFilter<FloatImageType, FloatImageType > DTFilter;
	typedef itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorTypeFloat;
	typedef itk::LabelGeometryImageFilter< USImageType > GeometryFilterType;
	typedef GeometryFilterType::LabelIndicesType labelindicestype;

	int sz_x, sz_y, sz_z;
	sz_x = input_image->GetLargestPossibleRegion().GetSize()[0];
	sz_y = input_image->GetLargestPossibleRegion().GetSize()[1];
	sz_z = input_image->GetLargestPossibleRegion().GetSize()[2];

	if( sz_x==1 || sz_y==1 || sz_z==1 ){
		LabelExtractType::Pointer deFilter = LabelExtractType::New();
		USImageType::RegionType dRegion = inp_labeled->GetLargestPossibleRegion();
		dRegion.SetSize(2,0);
		deFilter->SetExtractionRegion(dRegion);
		deFilter->SetInput( inp_labeled );
		deFilter->SetDirectionCollapseToIdentity();
		try{
			deFilter->Update();
		}
		catch( itk::ExceptionObject & excep ){
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
		//Dialate input first
		WholeCellSeg *dialate_filter = new WholeCellSeg;
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
		int qnum_sz = zp? (int)(labelsList.size()-1) : (int)(labelsList.size());
		for( int i=0; i<qnum_sz; ++i ){
			int counter=0;
			for( int j=0; j<number_of_rois; ++j ){
				if( quantified_numbers[(i*number_of_rois+j)] > 1 )
					++counter;
			}
			qfied_num.push_back((float)counter);
		}
	}
	else{
		GeometryFilterType::Pointer geomfilt1 = GeometryFilterType::New();
		geomfilt1->SetInput( inp_labeled );
		geomfilt1->SetCalculatePixelIndices( true );
		geomfilt1->Update();
		labelsList = geomfilt1->GetLabels();

		USImageType::SizeType  sizee;
		sizee[0] = inp_labeled->GetLargestPossibleRegion().GetSize()[0];  // size along X
		sizee[1] = inp_labeled->GetLargestPossibleRegion().GetSize()[1];  // size along Y
		sizee[2] = inp_labeled->GetLargestPossibleRegion().GetSize()[2];  // size along Z

		std::vector<double> quantified_numbers_cell;
		for( int j=0; j<(number_of_rois*(int)labelsList.size()); ++j ) quantified_numbers_cell.push_back((double)0.0);
		for( int i=0; i<labelsList.size(); ++i ){
			//Get label indices
			labelindicestype indices1;
			indices1 = geomfilt1->GetPixelIndices(labelsList[i]);

			//Create vnl array 3xN( label indicies )
			vnl_matrix<double> B(3,(int)indices1.size());

			//Get Centroid
			double centroid_x = (double)(geomfilt1->GetCentroid(labelsList[i])[0]);
			double centroid_y = (double)(geomfilt1->GetCentroid(labelsList[i])[1]);
			double centroid_z = (double)(geomfilt1->GetCentroid(labelsList[i])[2]);

			//Create an image with bounding box + 2 * outside distance + 2
			//and get distance map for the label
			GeometryFilterType::BoundingBoxType boundbox = geomfilt1->GetBoundingBox(labelsList[i]);
			int ssz_x = boundbox[1]-boundbox[0]+2*surr_dist+2;
			int ssz_y = boundbox[3]-boundbox[2]+2*surr_dist+2;
			int ssz_z = boundbox[5]-boundbox[4]+2*surr_dist+2;
			FloatImageType::Pointer inp_lab = FloatImageType::New();
			FloatImageType::PointType origint; origint[0] = 0; origint[1] = 0; origint[2] = 0;
			inp_lab->SetOrigin( origint );
			FloatImageType::IndexType startt;
			startt[0] = 0;  // first index on X
			startt[1] = 0;  // first index on Y
			startt[2] = 0;  // first index on Z
			FloatImageType::SizeType  sizet;
			sizet[0] = ssz_x;  // size along X
			sizet[1] = ssz_y;  // size along Y
			sizet[2] = ssz_z;  // size along Z
			FloatImageType::RegionType regiont;
			regiont.SetSize( sizet );
			regiont.SetIndex( startt );
			inp_lab->SetRegions( regiont );
			inp_lab->Allocate();
			inp_lab->FillBuffer(0.0);
			inp_lab->Update();
			IteratorTypeFloat iterator444 ( inp_lab, inp_lab->GetRequestedRegion() );

			//Populate matrix with deviations from the centroid for principal axes and
			//at the same time set up distance-transform computation
			int ind=0;
			for( labelindicestype::iterator itPixind = indices1.begin(); itPixind!=indices1.end(); ++itPixind, ++ind){
				IteratorType iterator3( input_image, input_image->GetRequestedRegion() );
				iterator3.SetIndex( *itPixind );
				B(0,ind) = iterator3.GetIndex()[0]-centroid_x;
				B(1,ind) = iterator3.GetIndex()[1]-centroid_y;
				B(2,ind) = iterator3.GetIndex()[2]-centroid_z;
				FloatImageType::IndexType cur_in;
				cur_in[0] = iterator3.GetIndex()[0]-boundbox[0]+1+surr_dist;
				cur_in[1] = iterator3.GetIndex()[1]-boundbox[2]+1+surr_dist;
				cur_in[2] = iterator3.GetIndex()[2]-boundbox[4]+1+surr_dist;
				iterator444.SetIndex( cur_in );
				iterator444.Set( 255.0 );
			}
			//Compute distance transform for the current object
			DTFilter::Pointer dt_obj= DTFilter::New() ;
			dt_obj->SetInput( inp_lab );
			dt_obj->SquaredDistanceOff();
			dt_obj->InsideIsPositiveOff();
			try{
				dt_obj->Update() ;
			} catch( itk::ExceptionObject & err ){
				std::cout << "Error in Distance Transform: " << err << std::endl;
			}
			FloatImageType::Pointer dist_im = dt_obj->GetOutput();

			//Use KLT to compute pricipal axes
			vnl_matrix<double> B_transp((int)indices1.size(),3);
			B_transp = B.transpose();
			vnl_matrix<double>  COV(3,3);
			COV = B * B_transp;
			double norm = 1.0/(double)indices1.size();
			COV = COV * norm;
			//Eigen decomposition
			vnl_real_eigensystem Eyegun( COV );
			vnl_matrix<vcl_complex<double> > EVals = Eyegun.D;
			double Eval1 = vnl_real(EVals)(0,0); double Eval2 = vnl_real(EVals)(1,1); double Eval3 = vnl_real(EVals)(2,2);
			vnl_double_3x3 EVectMat = Eyegun.Vreal;
			double V1[3],V2[3],EP_norm[3];
			if( Eval1 >= Eval3 && Eval2 >= Eval3 ){
				if( Eval1 >= Eval2 ){
					V1[0] = EVectMat(0,0); V1[1] = EVectMat(1,0); V1[2] = EVectMat(2,0);
					V2[0] = EVectMat(0,1); V2[1] = EVectMat(1,1); V2[2] = EVectMat(2,1);
				} else {
					V2[0] = EVectMat(0,0); V2[1] = EVectMat(1,0); V2[2] = EVectMat(2,0);
					V1[0] = EVectMat(0,1); V1[1] = EVectMat(1,1); V1[2] = EVectMat(2,1);
				}
			} else if( Eval1 >= Eval2 && Eval3 >= Eval2 ) {
				if( Eval1 >= Eval3 ){
					V1[0] = EVectMat(0,0); V1[1] = EVectMat(1,0); V1[2] = EVectMat(2,0);
					V2[0] = EVectMat(0,2); V2[1] = EVectMat(1,2); V2[2] = EVectMat(2,2);
				} else {
					V2[0] = EVectMat(0,0); V2[1] = EVectMat(1,0); V2[2] = EVectMat(2,0);
					V1[0] = EVectMat(0,2); V1[1] = EVectMat(1,2); V1[2] = EVectMat(2,2);
				}
			} else {
				if( Eval2 >= Eval3 ){
					V1[0] = EVectMat(0,1); V1[1] = EVectMat(1,1); V1[2] = EVectMat(2,1);
					V2[0] = EVectMat(0,2); V2[1] = EVectMat(1,2); V2[2] = EVectMat(2,2);
				} else {
					V2[0] = EVectMat(0,1); V2[1] = EVectMat(1,1); V2[2] = EVectMat(2,1);
					V1[0] = EVectMat(0,2); V1[1] = EVectMat(1,2); V1[2] = EVectMat(2,2);
				}
			}
			double n_sum = sqrt( V1[0]*V1[0]+V1[1]*V1[1]+V1[2]*V1[2] );
			V1[0] /= n_sum; V1[1] /= n_sum; V1[2] /= n_sum;
			n_sum = sqrt( V2[0]*V2[0]+V2[1]*V2[1]+V2[2]*V2[2] );
			V2[0] /= n_sum; V2[1] /= n_sum; V2[2] /= n_sum;
			//Get the normal to the plane formed by the biggest two EVs
			EP_norm[0] = V1[1]*V2[2]-V1[2]*V2[1];
			EP_norm[1] = V1[2]*V2[0]-V1[0]*V2[2];
			EP_norm[2] = V1[0]*V2[1]-V1[1]*V2[0];
			//Reassign V2 so that it is orthogonal to both EP_norm and V1
			V2[0] = V1[1]*EP_norm[2]-V1[2]*EP_norm[1];
			V2[1] = V1[2]*EP_norm[0]-V1[0]*EP_norm[2];
			V2[2] = V1[0]*EP_norm[1]-V1[1]*EP_norm[0];
			//Now we have the point normal form; EP_norm is the normal and
			//centroid_x, centroid_y, centroid_z is the normal
			//The equation to the plane is EP_norm[0](x-centroid_x)+EP_norm[1](y-centroid_y)+EP_norm[2](z-centroid_z)=0
			double dee = (centroid_x*EP_norm[0]+centroid_y*EP_norm[1]+centroid_z*EP_norm[2])*(-1.00);

			//Iterate through and assign values to each region
			typedef itk::ImageRegionConstIterator< FloatImageType > ConstIteratorTypeFloat;
			ConstIteratorTypeFloat pix_buf2( dist_im, dist_im->GetRequestedRegion() );
			IteratorType iterator44( input_image, input_image->GetRequestedRegion() );

			for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2 ){
				//Use pixels that are only within the defined radius from the nucleus
				if( pix_buf2.Get() <= (float)surr_dist ){
					USImageType::IndexType cur_in;
					double n_vec[3];
					cur_in[0] = pix_buf2.GetIndex()[0]+boundbox[0]-1-surr_dist;
					cur_in[1] = pix_buf2.GetIndex()[1]+boundbox[2]-1-surr_dist;
					cur_in[2] = pix_buf2.GetIndex()[2]+boundbox[4]-1-surr_dist;
					if( cur_in[0] < 0 || cur_in[1] < 0 || cur_in[2] < 0 ) continue;
					if( cur_in[0] >= sizee[0] || cur_in[1] >= sizee[1] || cur_in[2] >= sizee[2] ) continue;
					iterator44.SetIndex( cur_in );
					if( iterator44.Get() < thresh ) continue;

					//The projection of the point on the plane formed by the fist two major axes
					double xxx, yyy, zzz;
					xxx = (double)cur_in[0] - EP_norm[0]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
									/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					yyy = (double)cur_in[1] - EP_norm[1]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
									/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					zzz = (double)cur_in[2] - EP_norm[2]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
									/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					//The vector from the centroid to the projected point
					n_vec[0] = centroid_x-xxx;
					n_vec[1] = centroid_y-yyy;
					n_vec[2] = centroid_z-zzz;
					n_sum = sqrt( n_vec[0]*n_vec[0] + n_vec[1]*n_vec[1] + n_vec[2]*n_vec[2] );
					n_vec[0] /= n_sum; n_vec[1] /= n_sum; n_vec[2] /= n_sum;
					//n_vec is the normalized vect in the direction of the projected point
					//V1 is the largest eigenvector
					//Get the dot and cross product between the two
					double doooot, crooos,fin_est_angle;
					doooot = n_vec[0]*V1[0]+n_vec[1]*V1[1]+n_vec[2]*V1[2];
					crooos = n_vec[0]*V2[0]+n_vec[1]*V2[1]+n_vec[2]*V2[2];
					//if( crooos>0 )
					//	crooos = sqrt(1-(crooos*crooos));
					//else
					//	crooos = (-1.00)*(sqrt(1-(crooos*crooos)));//

					fin_est_angle = atan2( crooos, doooot );
					unsigned short bin_num, bin_low_num;
					bin_low_num = ceil((double)number_of_rois/2);
					if( fin_est_angle<0 )
						bin_num = floor(fin_est_angle/MM_PI*bin_low_num);
					else
						bin_num = bin_low_num+floor(abs(fin_est_angle)/MM_PI*(number_of_rois-bin_low_num));
					if( bin_num >= number_of_rois ) bin_num = number_of_rois-1;
					quantified_numbers_cell.at((i*number_of_rois+bin_num)) += iterator44.Get();
				}
			}
		}
		double meean = 0, std_deev = 0, thrrrr;
		for( uint64_t i=0; i<quantified_numbers_cell.size(); ++i )
			if( quantified_numbers_cell.at(i)>1 )
				meean += quantified_numbers_cell.at(i)/((double)quantified_numbers_cell.size());
		for( uint64_t i=0; i<quantified_numbers_cell.size(); ++i )
			if( quantified_numbers_cell.at(i)>1 )
				std_deev += ((meean-quantified_numbers_cell.at(i))*(meean-quantified_numbers_cell.at(i)))/
					    ((double)quantified_numbers_cell.size());
		std_deev = sqrt( std_deev );
		thrrrr = meean - 2 * std_deev;
		for( uint64_t i=0; i<(uint64_t)labelsList.size(); ++i ){
			int count = 0;
			for( uint64_t j=0; j<number_of_rois; ++j )
				if( quantified_numbers_cell.at((i*number_of_rois+j)) > thrrrr )
					++count;
			qfied_num.at(i) = count;
		}
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
		tempHist[i] = 0;

	//Populate the histogram (assume pixel type is actually uchar):
	ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		USPixelType pix = it.Get();
		if(pix <= 255)
			++tempHist[pix];
	}
	
	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<unsigned short>::max();
	double max = 0.0;
	for(int i=0; i<numBins; ++i)
		if( tempHist[i] > max )
			max = tempHist[i];

	double scaleFactor = 1;
	if(max >= floatIntegerMax)
		scaleFactor = floatIntegerMax / max;

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
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter ){
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
