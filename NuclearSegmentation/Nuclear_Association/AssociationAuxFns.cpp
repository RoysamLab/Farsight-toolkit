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

#define _USE_MATH_DEFINES
#include <math.h>

#include <algorithm>

#include "itkImage.h"
#include "itkIntTypes.h"
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

#include "ftkNuclearAssociationRules.h"

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 3 > USImageType;
typedef float FloatPixelType;
typedef itk::Image< FloatPixelType, 3 > FloatImageType;

std::vector<float> compute_ec_features( USImageType::Pointer input_image,  USImageType::Pointer inp_labeled, int number_of_rois, unsigned short thresh, int surr_dist, int inside_dist ){

	std::vector< float > qfied_num;
	std::vector< USImageType::PixelType > labelsList;
	std::vector< double > quantified_numbers_cell;

	typedef itk::ExtractImageFilter< USImageType, UShortImageType > LabelExtractType;
	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	typedef itk::SignedDanielssonDistanceMapImageFilter<FloatImageType, FloatImageType > DTFilter;
	typedef itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorTypeFloat;
	typedef itk::LabelGeometryImageFilter< USImageType > GeometryFilterType;
	typedef GeometryFilterType::LabelIndicesType labelindicestype;

	itk::SizeValueType sz_x, sz_y, sz_z;
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

		GeometryFilterType::Pointer geomfilt1 = GeometryFilterType::New();

		geomfilt1->SetInput( input_labeled );
		geomfilt1->SetCalculatePixelIndices( true );
		geomfilt1->Update();
		labelsList = geomfilt1->GetLabels();

		bool zp=false;
		for( USImageType::PixelType i=0; i < labelsList.size(); ++i ){
			if( labelsList[i] == 0 ){ zp=true; continue; }
			std::vector<float> quantified_numbers_cell;
			for( unsigned j=0; j<number_of_rois; ++j ) quantified_numbers_cell.push_back((float)0.0);
			double centroid_x = geomfilt1->GetCentroid(labelsList[i])[0];
			double centroid_y = geomfilt1->GetCentroid(labelsList[i])[1];
			labelindicestype indices1;
			indices1 = geomfilt1->GetPixelIndices(labelsList[i]);
			for( labelindicestype::iterator itPixind = indices1.begin(); itPixind!=indices1.end(); ++itPixind ){
				IteratorType iterator1 ( input_image, input_image->GetRequestedRegion() );
				iterator1.SetIndex( *itPixind );
				if( iterator1.Get() < thresh )
					continue;
				double x = iterator1.GetIndex()[0];
				double y = iterator1.GetIndex()[1];
				double angle = atan2((centroid_y-y),fabs(centroid_x-x));
				if( (centroid_x-x)>0 )
					angle += M_PI_2;
				else
					angle = M_PI+M_PI-(angle+M_PI_2);
				angle = ((number_of_rois-1)*angle)/(2*M_PI);
				double angle_fraction[1];
				unsigned angular_index;
				if( modf( angle, angle_fraction ) > 0.5 )
					angular_index = ceil( angle );
				else
					angular_index = floor( angle );
				
				quantified_numbers_cell[angular_index] += iterator1.Get();
			}
			for( unsigned j=0; j<number_of_rois; ++j ) quantified_numbers.push_back(quantified_numbers_cell[j]);
		}
		unsigned qnum_sz = zp? (labelsList.size()-1) : (labelsList.size());
		for( unsigned i=0; i<qnum_sz; ++i ){
			unsigned counter=0;
			for( unsigned j=0; j<number_of_rois; ++j ){
				if( quantified_numbers[(i*number_of_rois+j)] > (255.0*surr_dist) )
					++counter;
			}
			qfied_num.push_back(counter);
		}
	}
	else{
		GeometryFilterType::Pointer geomfilt1 = GeometryFilterType::New();
		geomfilt1->SetInput( inp_labeled );
		geomfilt1->SetCalculatePixelIndices( true );
		geomfilt1->Update();
		labelsList = geomfilt1->GetLabels();
		std::cout<<std::endl<<"The size is: "<<labelsList.size();
		if( labelsList.size() == 1 )
		{
			qfied_num.clear();
			return qfied_num;
		}
		bool zp=false; USPixelType zero;
		//Check if the background is also included
		for( USPixelType i=0; i<labelsList.size(); ++i ) if( labelsList[i] == 0 ){ zp=true; zero = i; }

		USImageType::SizeType  sizee;
		sizee[0] = inp_labeled->GetLargestPossibleRegion().GetSize()[0];  // size along X
		sizee[1] = inp_labeled->GetLargestPossibleRegion().GetSize()[1];  // size along Y
		sizee[2] = inp_labeled->GetLargestPossibleRegion().GetSize()[2];  // size along Z

		itk::SizeValueType roi_list_size = zp ?
					((itk::SizeValueType)number_of_rois*(labelsList.size()-1)*2) : 
					((itk::SizeValueType)number_of_rois*labelsList.size()*2);
		std::vector<double> quantified_numbers_cell((roi_list_size),0.0
			);
		std::cout<<"Bounding boxes computed"<<std::endl;

#ifdef _OPENMP
itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
#pragma omp parallel for
#if _OPENMP < 200805L
		for( int i=0; i<labelsList.size(); ++i )
#else
		for( USImageType::PixelType i=0; i<labelsList.size(); ++i )
#endif
#else
		for( USImageType::PixelType i=0; i<labelsList.size(); ++i )
#endif
		{
			itk::SizeValueType ind;
			if( zp && (zero==i) ) continue;
			if( zp && (i>zero)  ) ind = i-1;
			else ind = i;
			//Get label indices
			labelindicestype indices1;
			double centroid_x,centroid_y,centroid_z;
			GeometryFilterType::BoundingBoxType boundbox;
			#pragma omp critical
			{
				indices1 = geomfilt1->GetPixelIndices(labelsList[i]);
				//Get Centroid
				centroid_x = geomfilt1->GetCentroid(labelsList[i])[0];
				centroid_y = geomfilt1->GetCentroid(labelsList[i])[1];
				centroid_z = geomfilt1->GetCentroid(labelsList[i])[2];
				//Create an image with bounding box + 2 * outside distance + 2
				//and get distance map for the label
				boundbox = geomfilt1->GetBoundingBox(labelsList[i]);
			}
			//Create vnl array 3xN( label indicies )
			vnl_matrix<double> B(3,indices1.size());

			FloatImageType::Pointer inp_lab = FloatImageType::New();
			FloatImageType::PointType origint; origint[0] = 0; origint[1] = 0; origint[2] = 0;
			inp_lab->SetOrigin( origint );
			FloatImageType::IndexType startt;
			startt[0] = 0;  // first index on X
			startt[1] = 0;  // first index on Y
			startt[2] = 0;  // first index on Z
			FloatImageType::SizeType  sizet;
			sizet[0] = boundbox[1]-boundbox[0]+2*surr_dist+2;  // size along X
			sizet[1] = boundbox[3]-boundbox[2]+2*surr_dist+2;  // size along Y
			sizet[2] = boundbox[5]-boundbox[4]+2*surr_dist+2;  // size along Z
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
			itk::SizeValueType ind1=0;
			for( labelindicestype::iterator itPixind = indices1.begin(); itPixind!=indices1.end(); ++itPixind ){
				IteratorType iterator3( input_image, input_image->GetRequestedRegion() );
				iterator3.SetIndex( *itPixind );
				B(0,(ind1)) = iterator3.GetIndex()[0]-centroid_x;
				B(1,(ind1)) = iterator3.GetIndex()[1]-centroid_y;
				B(2,(ind1)) = iterator3.GetIndex()[2]-centroid_z;
				FloatImageType::IndexType cur_in;
				cur_in[0] = iterator3.GetIndex()[0]-boundbox[0]+1+surr_dist;
				cur_in[1] = iterator3.GetIndex()[1]-boundbox[2]+1+surr_dist;
				cur_in[2] = iterator3.GetIndex()[2]-boundbox[4]+1+surr_dist;
				iterator444.SetIndex( cur_in );
				iterator444.Set( 255.0 );
				++ind1;
			}

			//Compute distance transform for the current object
			DTFilter::Pointer dt_obj= DTFilter::New() ;
			dt_obj->SetInput( inp_lab );
			dt_obj->SquaredDistanceOff();
			dt_obj->InsideIsPositiveOff();
			try{
				dt_obj->Update() ;
			} catch( itk::ExceptionObject & err ){
				std::cerr << "Error in Distance Transform: " << err << std::endl;
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
			//centroid_x, centroid_y, centroid_z is the point
			//The equation to the plane is EP_norm[0](x-centroid_x)+EP_norm[1](y-centroid_y)+EP_norm[2](z-centroid_z)=0
			double dee = (centroid_x*EP_norm[0]+centroid_y*EP_norm[1]+centroid_z*EP_norm[2])*(-1.00);

			//Iterate through and assign values to each region
			typedef itk::ImageRegionConstIterator< FloatImageType > ConstIteratorTypeFloat;
			ConstIteratorTypeFloat pix_buf2( dist_im, dist_im->GetRequestedRegion() );
			IteratorType iterator44( input_image, input_image->GetRequestedRegion() );

			for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2 ){
				//Use pixels that are only within the defined radius from the nucleus
				double current_distance = pix_buf2.Get();
				if( (current_distance <= (double)surr_dist) && (current_distance>=(-1*inside_dist)) ){
					USImageType::IndexType cur_in;//,cur_in_cpy;
					double n_vec[3];
					cur_in[0] = pix_buf2.GetIndex()[0]+boundbox[0]-1-surr_dist;
					cur_in[1] = pix_buf2.GetIndex()[1]+boundbox[2]-1-surr_dist;
					cur_in[2] = pix_buf2.GetIndex()[2]+boundbox[4]-1-surr_dist;
					//cur_in_cpy[0] = cur_in[0]; cur_in_cpy[1] = cur_in[1]; cur_in_cpy[2] = cur_in[2];
					if( cur_in[0] < 0 || cur_in[1] < 0 || cur_in[2] < 0 ) continue;
					if( cur_in[0] >= sizee[0] || cur_in[1] >= sizee[1] || cur_in[2] >= sizee[2] ) continue;
					iterator44.SetIndex( cur_in );
					USImageType::PixelType pixel_intensity;
					pixel_intensity = iterator44.Get();
					if( pixel_intensity < thresh ) continue;

					//The projection of the point on the plane formed by the fist two major axes
					double xxx, yyy, zzz;
					xxx = cur_in[0] - EP_norm[0]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
									/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					yyy = cur_in[1] - EP_norm[1]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
									/(EP_norm[0]*EP_norm[0]+EP_norm[1]*EP_norm[1]+EP_norm[2]*EP_norm[2]));
					zzz = cur_in[2] - EP_norm[2]*((EP_norm[0]*cur_in[0]+EP_norm[1]*cur_in[1]+EP_norm[2]*cur_in[2]+dee)
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

					fin_est_angle = atan2( crooos, doooot );
					USPixelType bin_num;
					//Compute bin num
					if( fin_est_angle<0 )
						fin_est_angle += (2*M_PI);
					bin_num = floor(fin_est_angle*number_of_rois/(2*M_PI));

					//Check which side of the plane the point lies on
					double v_norm = (cur_in[0]-centroid_x)*(cur_in[0]-centroid_x)
									+(cur_in[1]-centroid_y)*(cur_in[1]-centroid_y)
									+(cur_in[2]-centroid_z)*(cur_in[2]-centroid_z);
					v_norm = sqrt( v_norm );
					double doot   = (cur_in[0]-centroid_x)*EP_norm[0]/v_norm + (cur_in[1]-centroid_y)*EP_norm[1]/v_norm + (cur_in[2]-centroid_z)*EP_norm[2]/v_norm;

					if( doot<0 )
						bin_num += number_of_rois;
					quantified_numbers_cell.at((ind*(2*number_of_rois)+bin_num)) += pixel_intensity;
				}
			}
		}
		number_of_rois = number_of_rois*2;
		if( labelsList.size() == 1 )
		{
			qfied_num.clear();
			return qfied_num;
		}
		std::vector<double> quantified_numbers_cell_cpy(roi_list_size);
		std::copy(quantified_numbers_cell.begin(), quantified_numbers_cell.end(), quantified_numbers_cell_cpy.begin() );
		//Run k-means
		//Most of the code is adapted from mul/mbl/mbl_k_means.cxx
		std::sort(quantified_numbers_cell.begin(), quantified_numbers_cell.end());
	
		bool skipKmeans;
		double Positive_thresh;
		if( quantified_numbers_cell.at(0) == quantified_numbers_cell.at(quantified_numbers_cell.size()-1) )
			skipKmeans = true;

		if( !skipKmeans ){
			std::cout<<"Starting k-means\n";
			std::cout<<"First:"<<quantified_numbers_cell.at(0)<<" Last:"<<quantified_numbers_cell.at(quantified_numbers_cell.size()-1)<<std::endl;
			unsigned k = 2;
			//Vectors and matrices for k-means
			std::vector< USImageType::PixelType > partition( roi_list_size, 0 );
			std::vector< double > sums   ( k, 0.0 );
			std::vector< double > centers( k, 0.0 );
			std::vector< USImageType::PixelType > nNearest( k, 0 );

			//Use the elements that are evenly spaced to get the intial centers
			for( unsigned i=0; i<k; ++i ){
				double index = ((double)(i)*roi_list_size)/(k+1);
				centers.at(i) = quantified_numbers_cell.at((itk::SizeValueType)index);
				bool duplicated;
				std::cout<<"Initializing centers\n"<<std::flush;
				if(i){
					if( centers.at((i-1)) == centers.at(i) ){
						duplicated = true;
						itk::SizeValueType ind=i+1;
						while( centers.at((i-1))==quantified_numbers_cell.at(ind) )
							++ind;
						centers.at(i) = quantified_numbers_cell.at(ind);
						sums.at(i)    = quantified_numbers_cell.at(ind);
					}
				}
				if(!duplicated)
					sums.at(i) = quantified_numbers_cell.at((i+1)/(k+1));
				++nNearest[i];
			}

			bool changed = true;
			std::cout<<"Waiting for kmeans to converge\n"<<std::flush;
			while(changed){
				changed = false;
				for(itk::SizeValueType i=0; i<roi_list_size; ++i){
					unsigned bestCentre = 0;
					double bestDist = fabs((centers.at(0)-quantified_numbers_cell.at(i)));
					for(unsigned j=1; j<k; ++j){
						double dist = fabs((centers.at(j)-quantified_numbers_cell.at(i)));
						if( dist < bestDist ){
							bestDist = dist; bestCentre = j;
						}
					}
					sums[bestCentre] += quantified_numbers_cell.at(i);
					++ nNearest[bestCentre];
					if( bestCentre != partition.at(i) ){
						changed = true;
						partition.at(i) = bestCentre;
					}
				}
				for( unsigned j=0; j<k; ++j) centers.at(j)  = sums.at(j)/nNearest.at(j);
				for( unsigned j=0; j<k; ++j) sums.at(j)     = 0;
				for( unsigned j=0; j<k; ++j) nNearest.at(j) = 0;
			}
			for( unsigned i=0; i<k; ++i )
				std::cout<<"Center "<<i<<" "<<centers.at(i)<<"\n";

			Positive_thresh = ((centers.at(0)+centers.at(1))/2) < (255.0*thresh)?
						 ((centers.at(0)+centers.at(1))/2) : (255.0*thresh); //Arbitrary upper thresh
		}
		else
			Positive_thresh = 255.0*thresh;

		std::cout<<"Positive_thresh "<<Positive_thresh<<"\n";

		std::cout<<"Done k-means\n";
		itk::SizeValueType ind = 0;
		for( USPixelType i=0; i<labelsList.size(); ++i ){
			if( zp && (zero==i) ) continue;
			int num_positive_rois = 0;
			for( unsigned j=0; j<number_of_rois; ++j ){
				itk::SizeValueType index_of_roi = ind*number_of_rois+j;
				if( quantified_numbers_cell_cpy.at(index_of_roi)>Positive_thresh )
					++num_positive_rois;
			}
			qfied_num.push_back(num_positive_rois);
			++ind;
		}
	}
	std::cout<<"Done surroundedness\n"<<std::flush;
	return qfied_num;
}

USImageType::PixelType returnthresh( itk::SmartPointer<USImageType> input_image,
			     int num_bin_levs, int num_in_fg ){
	//Instantiate the different image and filter types that will be used
	typedef itk::ImageRegionConstIterator< USImageType > ConstIteratorType;
	typedef itk::Statistics::Histogram< float > HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	std::cout<<"Starting threshold computation\n";

	//Create a temporary histogram container:
	const int numBins = itk::NumericTraits<USImageType::PixelType>::max();
	double *tempHist;
	tempHist = (double*) malloc( sizeof(double) * numBins );
	for(USImageType::PixelType i=0; i<numBins; ++i)
		tempHist[i] = 0;

	USImageType::PixelType maxval = itk::NumericTraits<USImageType::PixelType>::ZeroValue();
	USImageType::PixelType minval = itk::NumericTraits<USImageType::PixelType>::max();
	//Populate the histogram (assume pixel type is actually is some integer type):
	ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		USImageType::PixelType pix = it.Get();
		++tempHist[pix];
		if( pix > maxval ) maxval = pix;
		if( pix < minval ) minval = pix;
	}
	//return max of type if there is no variation in the staining
	if( (maxval-minval)<3 ) return itk::NumericTraits<USImageType::PixelType>::max(); 
	const USImageType::PixelType numBinsPresent = maxval+1;
	
	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<USImageType::PixelType>::max();
	double max = 0.0;
	for(USImageType::PixelType i=0; i<numBinsPresent; ++i)
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
	upperBound.Fill((double)maxval);
	size.Fill(numBinsPresent);

	histogram->SetMeasurementVectorSize(1);
	histogram->Initialize(size, lowerBound, upperBound ) ;

	USImageType::PixelType i=0;
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter ){
		float norm_freq = (float)(tempHist[i] * scaleFactor);
		iter.SetFrequency(norm_freq);
		++i;
	}

	std::cout<<"Histogram computed\n";

	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bin_levs );
	calculator->SetInputHistogram( histogram );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	float thresh;

	for(USImageType::PixelType i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<float>(*itNum));

	std::cout<<"Threshold computed: "<<thresh<<std::endl;

	return (USImageType::PixelType)(thresh+0.5);
}

#endif
