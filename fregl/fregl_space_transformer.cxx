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

#include "fregl_space_transformer.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include "itkImageRegionIterator.h"
#include <tinyxml/tinyxml.h>
#include <itkExtractImageFilter.h>
#include "itkMedianImageFunction.h"

static std::string ToString(double val);

template <class TPixel>
fregl_space_transformer< TPixel >::
fregl_space_transformer( )
{
	anchor_ = -1;
	anchor_set_ = false;

	//spacing_ is default to isotropic.
	spacing_[0] = 1;
	spacing_[1] = 1;
	spacing_[2] = 1;
}

	template <class TPixel>
fregl_space_transformer< TPixel >::fregl_space_transformer( typename fregl_joint_register< TPixel >::Pointer joint_reg )
{
	joint_register_ = joint_reg;
	anchor_ = -1;
	anchor_set_ = false;

	//spacing_ is default to isotropic.
	spacing_[0] = 1;
	spacing_[1] = 1;
	spacing_[2] = 1;
}

	template <class TPixel>
void fregl_space_transformer< TPixel >::set_anchor(std::string const &anchor_name, bool in_anchor, bool overlap_only, bool space_set)
{
	std::vector<std::string> const & image_names = joint_register_->image_names();

	int anchor_index = -1;
	for (unsigned i = 0; i < image_names.size(); i++) {
		std::cout << anchor_name << " " << image_names[i] << std::endl;
		if (anchor_name.compare(image_names[i]) == 0) 
			anchor_index = i;
	}

	if (anchor_index<0) {
		std::cerr<<"To_image is not found!!!"<<std::endl; 
		return ;
	}

	anchor_ = anchor_index;
	anchor_set_ = true;

	image_id_indices_.clear();
	image_origins_.clear();
	std::vector<SizeType> const & image_sizes = joint_register_->image_sizes();

	// Get valid image indices for mapping 
	for (unsigned int i = 0; i<joint_register_->number_of_images(); i++) {
		// If overalpping w/ the anchor is required, the image is only
		// considered if overlap w/ the anchor
		if (overlap_only && !joint_register_->is_overlapped(anchor_, i)) 
			continue;
		// If the image and the anchor do not belong to the same subgraph,
		// the image should not be considered.
		if (!joint_register_->in_same_subgraph(anchor_index, i))
			continue;

		image_id_indices_.push_back(i);
	}

	weight_images_2D_.resize(image_id_indices_.size(), 0);
	normalized_weight_images_2D_.resize(image_id_indices_.size(), 0);

	if ( in_anchor ) {
		origin_[0] = 0;
		origin_[1] = 0;
		origin_[2] = 0;
		roi_origin_[0] = 0;
		roi_origin_[1] = 0;
		roi_origin_[2] = 0;

		image_size_ = image_sizes[anchor_];
		roi_size_ = image_sizes[anchor_];
		std::cout<<"Origin = "<<origin_[0]<<","<<origin_[1]<<","<<origin_[2]<<std::endl;
		std::cout<<"Image_size = "<<image_size_[0]<<" x "<<image_size_[1]<<" x "<<image_size_[2]<<std::endl;

		return;  
	}

	if (space_set) {
		roi_origin_ = origin_;
		roi_size_ = image_size_;
		std::cout<<"Origin = "<<origin_[0]<<","<<origin_[1]<<","<<origin_[2]<<std::endl;
		std::cout<<"Image_size = "<<image_size_[0]<<" x "<<image_size_[1]<<" x "<<image_size_[2]<<std::endl;
		return;
	}

	// Compute the global space to determine the origin. Since the
	// transforms are affine, each bounding box of the transformed image
	// is determined by the 8 corners of the original image.
	PointType min_pt, max_pt;
	bool min_set = false, max_set = false; 

	for (unsigned int i = 0; i<joint_register_->number_of_images(); i++) {
		if (overlap_only && !joint_register_->is_overlapped(anchor_, i)) 
			continue;
		if (!joint_register_->in_same_subgraph(anchor_index, i))
			continue;

		TransformTypePointer xform = joint_register_->get_transform(i, anchor_);
		// Only images which can be transformed to the anchor space will
		// be considered
		if ( !xform ) continue;
		SizeType size = image_sizes[i];
		for (unsigned int sx = 0; sx < size[0]; sx+= size[0]-1)
			for (unsigned int sy = 0; sy < size[1]; sy+= size[1]-1)
				for (unsigned int sz = 0; sz < size[2]; sz+= size[2]-1) {
					PointType pt, xformed_pt;
					pt[0] = sx;
					pt[1] = sy;
					pt[2] = sz;
					//std::cout<<"Point = "<<pt[0]<<","<<pt[1]<<","<<pt[2]<<std::endl;

					xformed_pt = xform->TransformPoint(pt);
					//std::cout<<"xformed_Point = "<<xformed_pt[0]<<","<<xformed_pt[1]<<","<<xformed_pt[2]<<std::endl;
					if (sx == 0 && sy == 0 && sz==0){
						image_origins_.push_back(xformed_pt);
					}
					if (!min_set) {
						min_set = true;
						min_pt = xformed_pt;
					}
					else {
						if (xformed_pt[0] < min_pt[0]) min_pt[0] = xformed_pt[0];
						if (xformed_pt[1] < min_pt[1]) min_pt[1] = xformed_pt[1];
						if (xformed_pt[2] < min_pt[2]) min_pt[2] = xformed_pt[2];
					}
					if (!max_set) {
						max_pt = xformed_pt;
						max_set = true;
					}
					else {
						if (xformed_pt[0] > max_pt[0]) max_pt[0] = xformed_pt[0];
						if (xformed_pt[1] > max_pt[1]) max_pt[1] = xformed_pt[1];
						if (xformed_pt[2] > max_pt[2]) max_pt[2] = xformed_pt[2];
					}
				}
	}
	origin_ = min_pt;
	// Making the origin a integer position
	origin_[0] = floor(origin_[0]);
	origin_[1] = floor(origin_[1]);
	origin_[2] = floor(origin_[2]);
	roi_origin_[0] = floor(origin_[0]);
	roi_origin_[1] = floor(origin_[1]);
	roi_origin_[2] = floor(origin_[2]);
	//std::cout<<"Origin = "<<origin_[0]<<", "<<origin_[1]<<", "<<origin_[2]<<std::endl;

	image_size_[0] =  ceil(max_pt[0]-origin_[0] + 1);
	image_size_[1] =  ceil(max_pt[1]-origin_[1] + 1);
	image_size_[2] =  ceil(max_pt[2]-origin_[2] + 1);
	roi_size_[0] =  ceil(max_pt[0]-origin_[0] + 1);
	roi_size_[1] =  ceil(max_pt[1]-origin_[1] + 1);
	roi_size_[2] =  ceil(max_pt[2]-origin_[2] + 1);

	std::cout<<"Origin = "<<origin_[0]<<","<<origin_[1]<<","<<origin_[2]<<std::endl;
	std::cout<<"Image_size = "<<image_size_[0]<<" x "<<image_size_[1]<<" x "<<image_size_[2]<<std::endl;

}

//  Set the roi for processing
template <class TPixel>
void 
fregl_space_transformer< TPixel >::
set_roi(PointType s_origin, SizeType s_size) {

	roi_origin_[0] = s_origin[0];
	roi_origin_[1] = s_origin[1];
	roi_origin_[2] = s_origin[2];
	roi_size_[0] =  s_size[0];
	roi_size_[1] =  s_size[1];
	roi_size_[2] =  s_size[2];

	// Now make sure origin and size are in range or fix
	if (roi_origin_[0] < origin_[0]) roi_origin_[0] = origin_[0];
	if (roi_origin_[1] < origin_[1]) roi_origin_[1] = origin_[1];
	if (roi_origin_[2] < origin_[2]) roi_origin_[2] = origin_[2];
	if (roi_origin_[0] > origin_[0] + image_size_[0]) roi_origin_[0] = origin_[0] + image_size_[0];
	if (roi_origin_[1] > origin_[1] + image_size_[1]) roi_origin_[1] = origin_[1] + image_size_[1];
	if (roi_origin_[2] > origin_[2] + image_size_[2]) roi_origin_[2] = origin_[2] + image_size_[2];

	// Make sure size is in range or fix
	if (roi_origin_[0] + roi_size_[0] > origin_[0] + image_size_[0])
		roi_size_[0] = (origin_[0] + image_size_[0]) - roi_origin_[0];
	if (roi_origin_[1] + roi_size_[1] > origin_[1] + image_size_[1])
		roi_size_[1] = (origin_[1] + image_size_[1]) - roi_origin_[1];
	if (roi_origin_[1] + roi_size_[1] > origin_[1] + image_size_[1])
		roi_size_[1] = (origin_[1] + image_size_[1]) - roi_origin_[1];
	//	std::cout<<"ROI Origin = "<<roi_origin_[0]<<","<<roi_origin_[1]<<","<<roi_origin_[2]<<std::endl;
	//	std::cout<<"ROI Size = "<<roi_size_[0]<<" x "<<roi_size_[1]<<" x "<<roi_size_[2]<<std::endl;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::in_image(PointType loc, int image_index, 
		PointType& xformed_loc) const
{
	bool failed = false;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) return false;

	/*
	   if (inverse_xforms_.empty()) {
	   std::cerr<<" Inverse xforms not yet set! "<<std::endl;
	   failed = true;
	   }
	   */
	xformed_loc = inverse_xform->TransformPoint(loc);
	SizeType size = joint_register_->image_size(index);

	//Check if the point is in range
	if (xformed_loc[0]<0 || xformed_loc[0] >= size[0])
		failed = true;
	else if (xformed_loc[1]<0 || xformed_loc[1] >= size[1]) failed = true;
	else if (xformed_loc[2]<0 || xformed_loc[2] >= size[2]) failed = true;

	return !failed;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::
in_image(vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 3 > & xformed_loc) const
{
	PointType pt_loc, xformed_pt_loc;
	pt_loc[0] = loc[0];
	pt_loc[1] = loc[1];
	pt_loc[2] = loc[2];
	bool in_range = in_image(pt_loc, image_index, xformed_pt_loc);
	xformed_loc[0] = xformed_pt_loc[0];
	xformed_loc[1] = xformed_pt_loc[1];
	xformed_loc[2] = xformed_pt_loc[2];

	return in_range;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::
image_in_roi(int image_index) const {

	int index = image_id_indices_[image_index];
	TransformTypePointer xform = joint_register_->get_transform(index, anchor_);
	// If the image is not even in anchor space then return false
	if (!xform) return false;
	std::vector<SizeType> const & image_sizes = joint_register_->image_sizes();

	// Get the origin and size of the image in anchor space
	SizeType i_size = image_sizes[image_index];
	PointType i_origin = image_origins_[image_index];

	//Set the far corner of the image space
	PointType i_end = i_origin;
	i_end[0] += i_size[0];
	i_end[1] += i_size[1];
	i_end[2] += i_size[2];

	// Set the far corner of the roi
	PointType roi_end = roi_origin_;
	roi_end[0] += roi_size_[0];
	roi_end[1] += roi_size_[1];
	roi_end[2] += roi_size_[2];

	// Little complicated test to see if any piece of the image is in the anchor space.
	if ((i_end[0] >= roi_origin_[0] && i_end[1] >= roi_origin_[1] && i_end[2] >= roi_origin_[2]) && 
			(i_origin[0] <= roi_end[0] && i_origin[1] <= roi_end[1] && i_origin[2] <= roi_end[2])) return true;
	else return false;

}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::
in_image_2d(vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 2 > & xformed_loc) const
{
	bool failed = false;
	PointType pt_loc, xformed_pt_loc;
	pt_loc[0] = loc[0];
	pt_loc[1] = loc[1];
	pt_loc[2] = loc[2];
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) return false;

	xformed_pt_loc = inverse_xform->TransformPoint(pt_loc);
	SizeType size = joint_register_->image_size(index);

	//Check if the point is in range
	if (xformed_pt_loc[0]<0 || xformed_pt_loc[0] >= size[0])
		failed = true;
	else if (xformed_pt_loc[1]<0 || xformed_pt_loc[1] >= size[1]) failed = true;

	xformed_loc[0] = xformed_pt_loc[0];
	xformed_loc[1] = xformed_pt_loc[1];

	return !failed;
}

template <class TPixel>
bool
fregl_space_transformer< TPixel >::
in_range(PointType loc, int from_index, int to_index,
		PointType& xformed_loc) const
{
	bool failed = false;
	int to = image_id_indices_[to_index];
	int from = image_id_indices_[from_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(from,to);
	if ( !inverse_xform ) return false;

	xformed_loc = inverse_xform->TransformPoint(loc);
	SizeType size = joint_register_->image_size(to);

	//Check if the point is in range
	if (xformed_loc[0]<0 || xformed_loc[0] >= size[0])
		failed = true;
	else if (xformed_loc[1]<0 || xformed_loc[1] >= size[1]) failed = true;
	else if (xformed_loc[2]<0 || xformed_loc[2] >= size[2]) failed = true;

	return !failed;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::
in_range(vnl_vector_fixed< float, 3 > loc, int from_index, int to_index, 
		vnl_vector_fixed< float, 3 >& xformed_loc) const
{
	PointType pt_loc, xformed_pt_loc;
	pt_loc[0] = loc[0];
	pt_loc[1] = loc[1];
	pt_loc[2] = loc[2];
	bool is_in_range = in_range(pt_loc, from_index, to_index, xformed_pt_loc);
	xformed_loc[0] = xformed_pt_loc[0];
	xformed_loc[1] = xformed_pt_loc[1];
	xformed_loc[2] = xformed_pt_loc[2];

	return is_in_range;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::in_anchor(PointType loc, int image_index, 
		PointType& xformed_loc) const
{
	bool failed = false;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(index, anchor_);
	if ( !inverse_xform ) return false;

	/*
	   if (inverse_xforms_.empty()) {
	   std::cerr<<" Inverse xforms not yet set! "<<std::endl;
	   failed = true;
	   }
	   */
	xformed_loc = inverse_xform->TransformPoint(loc);
	SizeType size = joint_register_->image_size(anchor_);

	//Check if the point is in range
	if (xformed_loc[0]<0 || xformed_loc[0] >= size[0])
		failed = true;
	else if (xformed_loc[1]<0 || xformed_loc[1] >= size[1]) failed = true;
	else if (xformed_loc[2]<0 || xformed_loc[2] >= size[2]) failed = true;

	return !failed;
}

template <class TPixel>
bool 
fregl_space_transformer< TPixel >::
in_anchor(vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 3 > & xformed_loc) const
{
	PointType pt_loc, xformed_pt_loc;
	pt_loc[0] = loc[0];
	pt_loc[1] = loc[1];
	pt_loc[2] = loc[2];
	bool in_range = in_anchor(pt_loc, image_index, xformed_pt_loc);
	xformed_loc[0] = xformed_pt_loc[0];
	xformed_loc[1] = xformed_pt_loc[1];
	xformed_loc[2] = xformed_pt_loc[2];

	return in_range;
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageTypePointer 
fregl_space_transformer< TPixel >::
transform_image(ImageTypePointer in_image, int image_index, int background, bool use_NN_interpolator, bool normalizeImage, ImageTypePointer background_image, double sigma, double median)
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}

	ImageTypePointer normalize_image = NULL;
	if(normalizeImage)
	{
		if(background_image)
		{
			std::cout<< "Use background image."<<std::endl;
			intensity_normalize(in_image, background_image, normalize_image,  sigma, median);
		}
		else
		{
			std::cout<< "Generate background image."<<std::endl;
			intensity_normalize(in_image, normalize_image, sigma, median);
		}
	}
	/*
	   ImageTypePointer image = ImageType::New();

	   ImageType::IndexType start;
	   start[0] = 0;
	   start[1] = 0;
	   start[2] = 0;

	   ImageType::RegionType region;
	   region.SetSize( image_size_ );
	   region.SetIndex( start );

	   image->SetRegions( region );
	   image->Allocate();

	   image->SetOrigin( origin );
	   */

	// Set the resampler to generate the transformed image
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) 
		return NULL;

	std::cout<<"Transform image "<<joint_register_->image_names()[index]<<std::endl;

	typename TransformType::ParametersType params(12);
	params = inverse_xform->GetParameters();
	/*
	   std::cout<<"Inverse parameters = "<<params[0]<<", "<<params[1]<<", "
	   <<params[2]<<", "<<params[3]<<", "<<params[4]<<", "<<params[5]
	   <<", "<<params[6]<<", "<<params[7]<<", "<<params[8]<<", "
	   <<params[9]<<", "<<params[10]<<", "<<params[11]<<std::endl;
	   */
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	if (use_NN_interpolator) {
		typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpoType;
		typename NNInterpoType::Pointer nn_interpolator = NNInterpoType::New();
		resampler->SetInterpolator(nn_interpolator);
	}
	if(normalize_image)
	{
		std::cout<< "Set normalized image!"<<std::endl;
		resampler->SetInput(normalize_image);
	}
	else
	{
		resampler->SetInput(in_image);
	}
	resampler->SetTransform( inverse_xform );
	resampler->SetSize( image_size_ );
	resampler->SetOutputOrigin( origin_ );
	resampler->SetOutputSpacing(spacing_);
	resampler->SetDefaultPixelValue( background );
	try {
		resampler->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
		return NULL;
	}

	return resampler->GetOutput();
}

template <class TPixel>
int
fregl_space_transformer< TPixel >::
transform_image_fast(ImageTypePointer in_image, ImageTypePointer& out_image, PointType& offsetNew, int image_index, int background, bool use_NN_interpolator, bool normalizeImage, ImageTypePointer background_image, double sigma, double median)
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return 1;
	}

	ImageTypePointer normalize_image = NULL;
	if(normalizeImage)
	{
		if(background_image)
		{
			std::cout<< "Use background image."<<std::endl;
			intensity_normalize(in_image, background_image, normalize_image, sigma, median);
		}
		else
		{
			std::cout<< "Generate background image."<<std::endl;
			intensity_normalize(in_image, normalize_image, sigma, median);
		}
	}
	// Set the resampler to generate the transformed image
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) 
	{
		std::cerr<<"FIXME: NO TRANSFORMATION"<<std::endl;
		return 1;
	}

	std::cout<<"Transform image "<<joint_register_->image_names()[index]<<std::endl;

	typename TransformType::ParametersType params(12);
	params = inverse_xform->GetParameters();
	/*
	   std::cout<<"Inverse parameters = "<<params[0]<<", "<<params[1]<<", "
	   <<params[2]<<", "<<params[3]<<", "<<params[4]<<", "<<params[5]
	   <<", "<<params[6]<<", "<<params[7]<<", "<<params[8]<<", "
	   <<params[9]<<", "<<params[10]<<", "<<params[11]<<std::endl;
	   */
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	if (use_NN_interpolator) {
		typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpoType;
		typename NNInterpoType::Pointer nn_interpolator = NNInterpoType::New();
		resampler->SetInterpolator(nn_interpolator);
	}
	
	typename ImageType::SizeType image_size_fast;
	itk::uint64_t maxsize = sqrt(double((in_image->GetLargestPossibleRegion().GetSize()[0])*(in_image->GetLargestPossibleRegion().GetSize()[0])+(in_image->GetLargestPossibleRegion().GetSize()[1]*in_image->GetLargestPossibleRegion().GetSize()[1])+(in_image->GetLargestPossibleRegion().GetSize()[2]*in_image->GetLargestPossibleRegion().GetSize()[2])))+10; // Nicolas: Not sure if this is the maximum size, I have to recheck this
	image_size_fast[0] = maxsize*spacing_[0];
	image_size_fast[1] = maxsize*spacing_[1];
	image_size_fast[2] = maxsize*spacing_[2];
	
	itk::Vector<double, 3> oldOffset = inverse_xform->GetOffset();
	std::cout <<  "IMAGE SIZE: " <<image_size_fast <<std::endl;
 	std::cout << "OFFFSET: " << oldOffset<< std::endl;
 	std::cout << "ORIGIN: " << origin_<<std::endl;
	
	PointType origin_fast;
	origin_fast[0] = -(int)oldOffset[0];
	origin_fast[1] = -(int)oldOffset[1];
	origin_fast[2] = -(int)oldOffset[2];
	
 	std::cout << "ORIGIN FAST: " << origin_fast<< std::endl;
	
	if(normalize_image)
	{
		std::cout<< "Set normalized image!"<<std::endl;
		resampler->SetInput(normalize_image);
	}
	else
	{
		resampler->SetInput(in_image);
	}
	resampler->SetTransform( inverse_xform );
	resampler->SetSize( image_size_fast );
	resampler->SetOutputOrigin( origin_fast );
	resampler->SetOutputSpacing(spacing_);
	resampler->SetDefaultPixelValue( background );
	
	try {
		resampler->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
// 		return NULL;
	}
	out_image = resampler->GetOutput();

	offsetNew[0] = -origin_[0]-(int)oldOffset[0];
	offsetNew[1] = -origin_[1]-(int)oldOffset[1];
	offsetNew[2] = -origin_[2]-(int)oldOffset[2];
    return 0;
}

template <class TPixel>
bool fregl_space_transformer< TPixel >::
intensity_normalize(ImageTypePointer in_image, ImageTypePointer& out_image, double sigma, double median)
{
	std::cout<< "Intensity normalizing 1..."<<std::endl;
	int width = in_image->GetLargestPossibleRegion().GetSize()[0];
	int height = in_image->GetLargestPossibleRegion().GetSize()[1];
	int depth = in_image->GetLargestPossibleRegion().GetSize()[2];

	typedef itk::CastImageFilter<ImageType, FloatImageType> CasterImageToFloatFilter;
	typename CasterImageToFloatFilter::Pointer castFilter = CasterImageToFloatFilter::New();
	castFilter->SetInput(in_image);
	try
	{
		castFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in CasterImageToFloatFilter: " << err << std::endl; 
		return false;
	}
	FloatImageTypePointer image = castFilter->GetOutput();
	FloatImageType2DPointer backgroundImage = fregl_util< InputPixelType >::fregl_util_min_projection(image, sigma);
	
	//FloatImageType2DPointer probImage1 = ExtractSlice(image, depth - 1);
	//FloatImageType2DPointer probImage2 = ExtractSlice(image, depth - 2);

	//typedef itk::AddImageFilter< FloatImageType2D, FloatImageType2D, FloatImageType2D> AddImageFilterType;
	//typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
	//addFilter->SetInput1(probImage1);
	//addFilter->SetInput2(probImage2);

	//typedef itk::DivideImageFilter <FloatImageType2D, FloatImageType2D, FloatImageType2D> DivideImageFilterType;
	//typename DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New();
	//divideImageFilter->SetInput1(addFilter->GetOutput());
	//divideImageFilter->SetConstant2(2);

	//typedef itk::DiscreteGaussianImageFilter< FloatImageType2D, FloatImageType2D >  GaussinafilterType;
	//typename GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	//gaussianFilter->SetInput( divideImageFilter->GetOutput());
	//gaussianFilter->SetVariance(sigma);
	//try
	//{
	//	gaussianFilter->Update();
	//}
	//catch( itk::ExceptionObject &err)
	//{
	//	std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
	//	return false;
	//}

	typedef itk::ImageRegionConstIterator< FloatImageType2D > Prob2DConstIteratorType;
	typedef itk::ImageRegionIterator< FloatImageType>  ProbIteratorType;

	for(int i = 0; i < depth; i++)
	{
		Prob2DConstIteratorType backgroundIt( backgroundImage, backgroundImage->GetLargestPossibleRegion());
		backgroundIt.GoToBegin();
		typename FloatImageType::IndexType desiredStart;
		desiredStart[0] = 0;
		desiredStart[1] = 0;
		desiredStart[2] = i;
 
		typename FloatImageType::SizeType desiredSize;
		desiredSize[0] = width;
		desiredSize[1] = height;
		desiredSize[2] = 1;
 
		typename FloatImageType::RegionType desiredRegion(desiredStart, desiredSize);
		ProbIteratorType imageIt( image, desiredRegion);
		imageIt.GoToBegin();
		
		while( !backgroundIt.IsAtEnd() && !imageIt.IsAtEnd())
		{
			double val = (double)imageIt.Get() - (double)backgroundIt.Get();
			if( val < 1e-6)
			{
				val = 0;
			}
			imageIt.Set( val);
			++imageIt;
			++backgroundIt;
		}
	}

	FloatImageType2DPointer minprojAfterImage = fregl_util< InputPixelType >::fregl_util_min_projection(image, 0);
	// rescale
	typedef itk::StatisticsImageFilter<FloatImageType2D> StatisticsImageFilterType;
	typename StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	statisticsImageFilter->SetInput(minprojAfterImage);
	statisticsImageFilter->Update();
	double meanAfterDevide = statisticsImageFilter->GetMean();

	typedef itk::MultiplyImageFilter<FloatImageType, FloatImageType> MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
	multiplyImageFilter->SetInput(image);
	multiplyImageFilter->SetConstant(median / meanAfterDevide);
	std::cout<< "Multiply by: "<<median / meanAfterDevide<<std::endl;
	
	typedef itk::CastImageFilter<FloatImageType, ImageType> CasterFloatToImageFilter;
	typename CasterFloatToImageFilter::Pointer castFilter2 = CasterFloatToImageFilter::New();
	castFilter2->SetInput(multiplyImageFilter->GetOutput());
	try
	{
		castFilter2->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in CasterFloatToImageFilter: " << err << std::endl; 
		return false;
	}
	out_image = castFilter2->GetOutput();
	return true;
}

template <class TPixel>
bool fregl_space_transformer< TPixel >::
intensity_normalize(ImageTypePointer in_image, ImageTypePointer background_image, ImageTypePointer& out_image, double sigma, double magnify)
{
	std::cout<< "Intensity normalizing 2..."<<std::endl;
	int width = in_image->GetLargestPossibleRegion().GetSize()[0];
	int height = in_image->GetLargestPossibleRegion().GetSize()[1];
	int depth = in_image->GetLargestPossibleRegion().GetSize()[2];

	typedef itk::CastImageFilter<ImageType, FloatImageType> CasterImageToFloatFilter;
	typename CasterImageToFloatFilter::Pointer castFilter = CasterImageToFloatFilter::New();
	castFilter->SetInput(in_image);
	try
	{
		castFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in CasterImageToFloatFilter: " << err << std::endl; 
		return false;
	}
	FloatImageTypePointer image = castFilter->GetOutput();

	typename CasterImageToFloatFilter::Pointer castFilter1 = CasterImageToFloatFilter::New();
	castFilter1->SetInput(background_image);
	try
	{
		castFilter1->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in CasterImageToFloatFilter: " << err << std::endl; 
		return false;
	}
	FloatImageTypePointer backimage = castFilter1->GetOutput();
	FloatImageType2DPointer backimage2D = ExtractSlice(backimage, 0);

	typedef itk::DiscreteGaussianImageFilter< FloatImageType2D, FloatImageType2D >  GaussinafilterType;
	typename GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	gaussianFilter->SetInput( backimage2D);
	gaussianFilter->SetVariance(sigma);
	try
	{
		gaussianFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
		return false;
	}

	FloatImageType2DPointer backgroundImage = gaussianFilter->GetOutput();

	typedef itk::ImageRegionConstIterator< FloatImageType2D > Prob2DConstIteratorType;
	typedef itk::ImageRegionIterator< FloatImageType>  ProbIteratorType;

	for(int i = 0; i < depth; i++)
	{
		Prob2DConstIteratorType backgroundIt( backgroundImage, backgroundImage->GetLargestPossibleRegion());
		backgroundIt.GoToBegin();
		typename FloatImageType::IndexType desiredStart;
		desiredStart[0] = 0;
		desiredStart[1] = 0;
		desiredStart[2] = i;
 
		typename FloatImageType::SizeType desiredSize;
		desiredSize[0] = width;
		desiredSize[1] = height;
		desiredSize[2] = 1;
 
		typename FloatImageType::RegionType desiredRegion(desiredStart, desiredSize);
		ProbIteratorType imageIt( image, desiredRegion);
		imageIt.GoToBegin();
		
		while( !backgroundIt.IsAtEnd() && !imageIt.IsAtEnd())
		{
			double val = imageIt.Get();
			if( backgroundIt.Get() > 1e-6)
			{
				val = (double)imageIt.Get() / (double)backgroundIt.Get();
			}
			else
			{
				val = 1;
			}
			imageIt.Set( val);
			++imageIt;
			++backgroundIt;
		}
	}

	FloatImageType2DPointer probImageAfterDevide = ExtractSlice(image, depth - 1);

	typedef itk::MultiplyImageFilter<FloatImageType, FloatImageType> MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
	multiplyImageFilter->SetInput(image);
	multiplyImageFilter->SetConstant(magnify);
	
	typedef itk::CastImageFilter<FloatImageType, ImageType> CasterFloatToImageFilter;
	typename CasterFloatToImageFilter::Pointer castFilter2 = CasterFloatToImageFilter::New();
	castFilter2->SetInput(multiplyImageFilter->GetOutput());
	try
	{
		castFilter2->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in CasterFloatToImageFilter: " << err << std::endl; 
		return false;
	}
	out_image = castFilter2->GetOutput();
	return true;
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::FloatImageType2DPointer 
fregl_space_transformer< TPixel >::
ExtractSlice(FloatImageTypePointer image, int sliceId)
{
	typename FloatImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = sliceId;
	typename FloatImageType::SizeType size;
	size[0] = image->GetLargestPossibleRegion().GetSize()[0];
	size[1] = image->GetLargestPossibleRegion().GetSize()[1];
	size[2] = 0;
	typename FloatImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	typedef itk::ExtractImageFilter< FloatImageType, FloatImageType2D> ExtractFilterType2D;
	typename ExtractFilterType2D::Pointer extractFilter = ExtractFilterType2D::New();
	extractFilter->SetInput(image);
	extractFilter->SetExtractionRegion(desiredRegion);
#if ITK_VERSION_MAJOR >= 4
	extractFilter->SetDirectionCollapseToIdentity(); // This is required.
#endif
	extractFilter->Update();
	FloatImageType2DPointer slice = extractFilter->GetOutput();
	return slice;
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageTypePointer 
fregl_space_transformer< TPixel >::
transform_image_roi(ImageTypePointer in_image, int image_index, int background, bool use_NN_interpolator ) const
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}


	// Set the resampler to generate the transformed image
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) 
		return NULL;

	//	std::cout<<"Transform image in ROI "<<joint_register_->image_names()[image_index]<<std::endl;

	typename TransformType::ParametersType params(12);
	params = inverse_xform->GetParameters();
	/*
	   std::cout<<"Inverse parameters = "<<params[0]<<", "<<params[1]<<", "
	   <<params[2]<<", "<<params[3]<<", "<<params[4]<<", "<<params[5]
	   <<", "<<params[6]<<", "<<params[7]<<", "<<params[8]<<", "
	   <<params[9]<<", "<<params[10]<<", "<<params[11]<<std::endl;
	   */
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	if (use_NN_interpolator) {
		typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpoType;
		typename NNInterpoType::Pointer nn_interpolator = NNInterpoType::New();
		resampler->SetInterpolator(nn_interpolator);
	}
	resampler->SetInput(in_image);
	resampler->SetTransform( inverse_xform );
	resampler->SetSize( roi_size_ );
	resampler->SetOutputOrigin( roi_origin_ );
	resampler->SetOutputSpacing(spacing_);
	resampler->SetDefaultPixelValue( background );
	try {
		resampler->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
		return NULL;
	}

	return resampler->GetOutput();
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageTypePointer 
fregl_space_transformer< TPixel >::
transform_image_whole(ImageTypePointer in_image, int image_index, int background, bool use_NN_interpolator ) const
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}


	// Set the resampler to generate the transformed image
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	int index = image_id_indices_[image_index];
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if ( !inverse_xform ) 
		return NULL;

	//	std::cout<<"Transform Whole Image "<<joint_register_->image_names()[image_index]<<std::endl;

	typename TransformType::ParametersType params(12);
	params = inverse_xform->GetParameters();
	/*
	   std::cout<<"Inverse parameters = "<<params[0]<<", "<<params[1]<<", "
	   <<params[2]<<", "<<params[3]<<", "<<params[4]<<", "<<params[5]
	   <<", "<<params[6]<<", "<<params[7]<<", "<<params[8]<<", "
	   <<params[9]<<", "<<params[10]<<", "<<params[11]<<std::endl;
	   */
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	if (use_NN_interpolator) {
		typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpoType;
		typename NNInterpoType::Pointer nn_interpolator = NNInterpoType::New();
		resampler->SetInterpolator(nn_interpolator);
	}
	// Get the origin and size of the image in anchor space
	std::vector<SizeType> const & image_sizes = joint_register_->image_sizes();
	SizeType i_size = image_sizes[image_index];
	PointType i_origin = image_origins_[image_index];

	resampler->SetInput(in_image);
	resampler->SetTransform( inverse_xform );
	resampler->SetSize( i_size );
	resampler->SetOutputOrigin( i_origin );
	resampler->SetOutputSpacing(spacing_);
	resampler->SetDefaultPixelValue( background );
	try {
		resampler->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
		return NULL;
	}

	return resampler->GetOutput();
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageTypePointer 
fregl_space_transformer< TPixel >::
transform_image_weighted(ImageTypePointer image, int image_index, int background, bool use_NN_interpolator ) const
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}

	int index = image_id_indices_[image_index];
	std::cout<<"Computing photo-bleached weighted image for "<<joint_register_->image_names()[index]<<std::endl;
	TransformTypePointer inverse_xform = joint_register_->get_transform(anchor_,index);
	if (!inverse_xform) return NULL;

	// Set the resampler to generate the transformed image
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpoType;
	typedef itk::NearestNeighborInterpolateImageFunction<ImageType2D, double> NNInterpoType2D;
	typedef itk::LinearInterpolateImageFunction<ImageType2D, double> LInterpoType2D;
	typedef itk::InterpolateImageFunction<ImageType2D, double> InterpoType2D;

	typename TransformType::ParametersType params(12);

	params = inverse_xform->GetParameters();
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	if (use_NN_interpolator) {
		typename NNInterpoType::Pointer nn_interpolator = NNInterpoType::New();
		resampler->SetInterpolator(nn_interpolator);
	}
	resampler->SetInput(image);
	resampler->SetTransform( inverse_xform );
	resampler->SetSize( image_size_ );
	resampler->SetOutputOrigin( origin_ );
	resampler->SetOutputSpacing(spacing_);
	resampler->SetDefaultPixelValue( background );
	try {
		resampler->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
		return NULL;
	}
	ImageTypePointer xformed_image = resampler->GetOutput();

	// Set the interpolator
	typename InterpoType2D::Pointer interpolator;
	if ( use_NN_interpolator ) 
		interpolator = NNInterpoType2D::New();
	else interpolator = LInterpoType2D::New();
	interpolator->SetInputImage( normalized_weight_images_2D_[image_index] );

	// weight the xformed image using the photo-bleaching factors. 
	PointType pt, xformed_pt;
	IndexType pix;
	typename ImageType2D::PointType xformed_pt2d ;
	typename ImageType2D::IndexType pix2d;
	typedef itk::ImageRegionIterator< ImageType > RegionIterator;
	RegionIterator outputIt( xformed_image, xformed_image->GetRequestedRegion() );
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		if ( outputIt.Get() == 0 ) continue;

		pix = outputIt.GetIndex();
		xformed_image->TransformIndexToPhysicalPoint( pix, pt );
		if ( in_image(pt, image_index, xformed_pt) ) {
			xformed_pt2d[0] = xformed_pt[0];
			xformed_pt2d[1] = xformed_pt[1];
			normalized_weight_images_2D_[image_index]->TransformPhysicalPointToIndex( xformed_pt2d, pix2d);

			// use linear interpolation
			double weight = 0;
			if( interpolator->IsInsideBuffer(pix2d) ) {
				weight = interpolator->EvaluateAtContinuousIndex(pix2d);
				outputIt.Set( int(outputIt.Get()*weight/255) );
			}
		}
	}

	return xformed_image;
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageTypePointer 
fregl_space_transformer< TPixel >::
compute_weighted_image() const
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}

	std::cout<<"Computing the weight image"<<std::endl;
	ImageTypePointer weight_image = ImageType::New();
	typename ImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	typename ImageType::RegionType region;
	region.SetSize(image_size_);
	region.SetIndex( start );
	weight_image->SetRegions( region );
	weight_image->Allocate();
	//FillBuffer does not work for large images. Index over-flow.
	//weight_image->FillBuffer( itk::NumericTraits<unsigned char>::Zero );
	typedef itk::ImageRegionIterator< ImageType > RegionIterator;
	RegionIterator outputIt( weight_image, weight_image->GetRequestedRegion() );
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		outputIt.Set(0);
	}
	weight_image->SetSpacing( spacing_ );
	weight_image->SetOrigin( origin_ );
	/*
	// Transform images one by one to the weight space
	PointType loc, xformed_loc;
	IndexType pix;
	unsigned char count;
	std::vector<std::string> const& names = image_names();
	for (unsigned int i = 0; i<image_id_indices_.size(); i++){
	std::cout<<"Computing the weight contribution from image "<<names[i]<<std::endl;
	int index = image_id_indices_[i];
	SizeType size = joint_register_->image_size(index);
	TransformTypePointer xform = joint_register_->get_transform(index, anchor_);
	for (unsigned long xi = 0; xi<size[0]; xi++) 
	for (unsigned long yi = 0; yi<size[1]; yi++) 
	for (unsigned long zi = 0; zi<size[2]; zi++) {
	loc[0] = xi;
	loc[1] = yi;
	loc[2] = zi;
	// no need to check out of range, since the global space
	// should contains all xformed points from the set of
	// images.
	xformed_loc = xform->TransformPoint(loc);
	weight_image->TransformPhysicalPointToIndex( xformed_loc, pix );
	count = weight_image->GetPixel(pix);
	weight_image->SetPixel(pix, count+1);
	} 
	}
	*/
	//typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
	//typedef itk::ImageRegionIterator< ImageType > RegionIterator;
	//RegionIterator outputIt( weight_image, weight_image->GetRequestedRegion() );
	PointType pt, xformed_pt;
	IndexType pix;
	unsigned char count = 0;

	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		pix = outputIt.GetIndex();
		weight_image->TransformIndexToPhysicalPoint( pix, pt );
		count = 0;
		for (unsigned int i = 0; i<image_id_indices_.size(); i++){
			if ( in_image(pt, i, xformed_pt) )
				count++;
		}
		outputIt.Set(count);
	}
	std::cout<<"Computed the weight image"<<std::endl;
	return weight_image;
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageType2DPointer
fregl_space_transformer< TPixel >::
compute_weighted_image_2D() const
{
	if (!anchor_set_) {
		std::cerr<<"Set anchor first"<<std::endl;
		return NULL;
	}

	std::cout<<"Computing the weight image"<<std::endl;
	ImageType2DPointer weight_image = ImageType2D::New();
	typename ImageType2D::IndexType start;
	start[0] = 0;
	start[1] = 0;
	typename ImageType2D::SizeType size;
	size[0] = image_size_[0];
	size[1] = image_size_[1];
	typename ImageType2D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	weight_image->SetRegions( region );
	weight_image->Allocate();
	//FillBuffer does not work for large images. Index over-flow.
	typedef itk::ImageRegionIterator< ImageType2D > RegionIterator;
	RegionIterator outputIt( weight_image, weight_image->GetRequestedRegion() );
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		outputIt.Set(0);
	}
	typename ImageType2D::PointType origin;
	origin[0] = origin_[0];
	origin[1] = origin_[1];
	weight_image->SetOrigin( origin );
	typename ImageType2D::SpacingType spacing;
	spacing[0] = spacing_[0];
	spacing[1] = spacing_[1];
	weight_image->SetSpacing( spacing );


	// Counting the number of images overlapping at each pixel of the
	// middle slice in the global space.
	PointType pt, xformed_pt;
	typename ImageType2D::IndexType pix2d;
	typename ImageType2D::PointType pt2d;
	unsigned char count = 0;
	float z_pos = origin_[2]+(spacing_[2]* image_size_[2]/2.0);//middle of the stack
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		pix2d = outputIt.GetIndex();
		weight_image->TransformIndexToPhysicalPoint( pix2d, pt2d );

		count = 0;
		for (unsigned int i = 0; i<image_id_indices_.size(); i++){
			pt[0] = pt2d[0];
			pt[1] = pt2d[1];
			pt[2] = z_pos;
			if ( in_image(pt, i, xformed_pt) )
				count++;
		}
		outputIt.Set(count);
	}
	std::cout<<"Computed the 2D weight image"<<std::endl;
	return weight_image;
}

template <class TPixel>
std::vector<std::string>
fregl_space_transformer< TPixel >::
image_names() const
{
	std::vector<std::string> names;
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		names.push_back(joint_register_->image_name(image_id_indices_[i]));
	}

	return names;
}

template <class TPixel>
std::vector<typename fregl_space_transformer< TPixel >::SizeType >
fregl_space_transformer< TPixel >::
image_sizes() const
{
	std::vector<SizeType> sizes;
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		sizes.push_back(joint_register_->image_size(image_id_indices_[i]));
	}

	return sizes;
}

template <class TPixel>
std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer> 
fregl_space_transformer< TPixel >::
xforms_to_neighbors(int image_index) const
{
	std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer>
    xforms;
	int from = image_id_indices_[image_index];
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		int to = image_id_indices_[i];
		if ( joint_register_->is_overlapped(from, to)&& from!=to )
			xforms.push_back( joint_register_->get_transform(from, to) );
	}
	return xforms;
}

template <class TPixel>
std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer> 
fregl_space_transformer< TPixel >::
xforms_from_neighbors(int image_index) const
{
	std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer>
    xforms;
	int to = image_id_indices_[image_index];
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		int from = image_id_indices_[i];
		if ( joint_register_->is_overlapped(from, to) && from!=to)
			xforms.push_back( joint_register_->get_transform(from, to) );
	}
	return xforms;
}

template <class TPixel>
std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer> 
fregl_space_transformer< TPixel >::
xforms_to_all(int image_index) const
{
	std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer>
    xforms;
	int from = image_id_indices_[image_index];
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		int to = image_id_indices_[i];   
		if ( from==to ) continue;
		xforms.push_back( joint_register_->get_transform(from, to) );
	}
	return xforms;
}

template <class TPixel>
std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer> 
fregl_space_transformer< TPixel >::
xforms_from_all(int image_index) const
{
	std::vector<typename fregl_space_transformer< TPixel >::TransformTypePointer>
    xforms;
	int to = image_id_indices_[image_index];
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		int from = image_id_indices_[i];   
		if ( from==to ) continue;
		xforms.push_back( joint_register_->get_transform(from, to) );
	}
	return xforms;
}

template <class TPixel>
void 
fregl_space_transformer< TPixel >::
set_individual_weight_map(int index, ImageTypePointer image, float alpha)
{
	typedef itk::ImageRegionIterator< ImageType2D > RegionIterator2D;
	typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;

	//float sigma = 1.0;
	ImageType2DPointer max_image = max_projection( image );
	std::vector<TransformTypePointer> xforms = xforms_from_neighbors( index );
	ImageType2DPointer weight_image = ImageType2D::New();
	weight_image->SetRegions(max_image->GetLargestPossibleRegion());

	/*
	   ImageType2D::IndexType start;
	   start[0] = region_3D.GetIndex()[0];
	   start[1] = region_3D.GetIndex()[1];
	   ImageType2D::SizeType size;
	   size[0] = region_3D.GetSize()[0];
	   size[1] = region_3D.GetSize()[1];

	   ImageType2D::RegionType region;
	   region.SetSize( size );
	   region.SetIndex( start );
	   weight_image->SetRegions( region );
	   */
	weight_image->Allocate();
	typename ImageType2D::RegionType region2D = weight_image->GetRequestedRegion();
	//  ImageType2D::RegionType region3D = image->GetRequestedRegion();

	//FillBuffer does not work for large images. Index over-flow.
	RegionIterator2D outputIt( weight_image, weight_image->GetRequestedRegion() );
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		outputIt.Set(255);
	}

	// Compute the average intensity of the central area of the
	// max_projection. This is computed as the intensity value which
	// corresponds to the peak of the histogram.
	//
	typename ImageType2D::IndexType center_start;
	center_start[0] = region2D.GetIndex()[0]+region2D.GetSize()[0]*0.25;
	center_start[1] = region2D.GetIndex()[1]+region2D.GetSize()[1]*0.25;
	typename ImageType2D::SizeType center_size;
	center_size[0] = region2D.GetSize()[0]*0.5;
	center_size[1] = region2D.GetSize()[1]*0.5;
	typename ImageType2D::RegionType center_region;
	center_region.SetIndex(center_start);
	center_region.SetSize(center_size);
	max_image->SetRequestedRegion(center_region);
	std::vector< long > histogram(256, 0);
	RegionIterator2D centerIt( max_image, max_image->GetRequestedRegion() );
	for ( centerIt.GoToBegin(); !centerIt.IsAtEnd(); ++centerIt) {
		if (centerIt.Get() > 30)
			histogram[centerIt.Get()]++;
	}
	int max = histogram[0];
	int p_center = 0;
	for (unsigned int j = 1; j<histogram.size(); j++) {
		if (histogram[j] > max) {
			max = histogram[j];
			p_center = j;
		}
	}
	if (p_center == 0) p_center = 1;
	std::cout<<"p_center = "<<p_center<<std::endl;

	// For each adjacent images, identify the overlapping volume, and
	// identify the peak of the histogram
	for (unsigned int i = 0; i<xforms.size(); i++) {
		typename TransformType::ParametersType params = xforms[i]->GetParameters();
		double tx = params[9];
		double ty = params[10];
		std::cout<<"tx = "<<tx<<"ty = "<<ty<<std::endl;
		typename ImageType2D::RegionType region_i = max_image->GetLargestPossibleRegion();
		typename ImageType2D::IndexType start_i, end_i; 
		typename ImageType2D::SizeType size_i = region_i.GetSize();

		// Find the bounding box for the overlap area
		start_i[0] = vnl_math_max(int(tx), 0);
		start_i[1] = vnl_math_max(int(ty), 0);
		end_i[0] = vnl_math_min(int(size_i[0]+tx), int(size_i[0]));
		end_i[1] = vnl_math_min(int(size_i[1]+ty), int(size_i[1]));
		size_i[0] = end_i[0] - start_i[0];
		size_i[1] = end_i[1] - start_i[1];
		region_i.SetIndex(start_i);
		region_i.SetSize(size_i);
		max_image->SetRequestedRegion(region_i);
		region_i.Print(std::cout);
		// Compute the histogram info
		histogram.clear();
		histogram.resize(256, 0);
		RegionIterator2D overlapIt( max_image, max_image->GetRequestedRegion() );
		for ( overlapIt.GoToBegin(); !overlapIt.IsAtEnd(); ++overlapIt) {
			if (overlapIt.Get() > 30)
				histogram[overlapIt.Get()]++;
		}

		int max = histogram[0];
		int p_i = 0;
		for (unsigned int j = 1; j<histogram.size(); j++) {
			if (histogram[j] > max) {
				max = histogram[j];
				p_i = j;
			}
		}
		std::cout<<"p_i = "<<p_i<<std::endl;
		weight_image->SetRequestedRegion( region_i );
		float ratio = std::pow(vnl_math_min((float)1.0, (float)p_i/(p_center)), alpha);
		std::cout<<"ratio = "<<ratio<<std::endl;
		RegionIterator2D weightIt( weight_image, weight_image->GetRequestedRegion() );
		for ( weightIt.GoToBegin(); !weightIt.IsAtEnd(); ++weightIt) {
			//weightIt.Set( vnl_math_min(255*ratio, float(weightIt.Get())) );
			weightIt.Set( weightIt.Get()*ratio );
		}    
	}

	// Gaussian smooth the weight_image
	typedef itk::DiscreteGaussianImageFilter< ImageType2D,FloatImageType2D > SmoothingFilterType;
	typedef itk::CastImageFilter< FloatImageType2D, ImageType2D > CastFilterType;
	typename SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
	typename CastFilterType::Pointer caster = CastFilterType::New();
	smoother->SetInput( weight_image );
	smoother->SetVariance(2.0);
	smoother->SetMaximumKernelWidth(15);
	caster->SetInput( smoother->GetOutput() );
	try {
		caster->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
	}
	weight_images_2D_[index]=caster->GetOutput();
}

template <class TPixel>
void 
fregl_space_transformer< TPixel >::
normalize_individual_weight_maps()
{
	typedef itk::ImageRegionConstIterator< ImageType2D > RegionConstIterator2D;
	typedef itk::ImageRegionIterator< ImageType2D > RegionIterator2D;

	std::vector<SizeType> sizes = image_sizes();
	for (unsigned int index = 0; index<normalized_weight_images_2D_.size(); index++) {    
		if (!weight_images_2D_[index]) continue;

		ImageType2DPointer image = ImageType2D::New();
		image->SetRegions( weight_images_2D_[index]->GetLargestPossibleRegion() );
		image->Allocate();
		typename ImageType2D::PointType origin;
		origin[0] = 0;
		origin[1] = 0;
		image->SetOrigin( origin );
		typename ImageType2D::SpacingType spacing;
		spacing[0] = spacing_[0];
		spacing[1] = spacing_[1];
		image->SetSpacing( spacing );

		// Initialized image
		RegionConstIterator2D inputIt( weight_images_2D_[index], weight_images_2D_[index]->GetLargestPossibleRegion() );
		RegionIterator2D outputIt( image, image->GetLargestPossibleRegion() );
		for ( outputIt.GoToBegin(), inputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++inputIt) {
			outputIt.Set( inputIt.Get() );
		}

		// For each pixel, sum up the weights of corresponding pixels in the
		// neighboring images.
		PointType pt, xformed_pt;
		typename ImageType2D::IndexType pix2d;
		typename ImageType2D::PointType pt2d,xformed_pt2d ;
		float z_pos = spacing_[2] * sizes[index][2]/2.0;
		for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
			pix2d = outputIt.GetIndex();
			image->TransformIndexToPhysicalPoint( pix2d, pt2d );
			int total = 0;
			for (unsigned k=0; k<sizes.size(); k++) {
				if (!weight_images_2D_[k]) continue;

				pt[0] = pt2d[0];
				pt[1] = pt2d[1];
				pt[2] = z_pos;
				if ( in_range(pt, index, k, xformed_pt) ) {
					// get the weight from the weight map
					xformed_pt2d[0] = xformed_pt[0];
					xformed_pt2d[1] = xformed_pt[1];
					weight_images_2D_[k]->TransformPhysicalPointToIndex( xformed_pt2d, pix2d);
					total += weight_images_2D_[k]->GetPixel( pix2d );
				}
			}

			// The weight is scaled by 255 to avoid floating point
			outputIt.Set( vnl_math_rnd(outputIt.Get()/(float)total * 255) );
		}
		normalized_weight_images_2D_[index] = image;
	}

	// The memory for the weight maps are released
	weight_images_2D_.clear();
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageType2D::Pointer
fregl_space_transformer< TPixel >::
get_weight_map(int index) const 
{
	assert(normalized_weight_images_2D_[index]);
	return normalized_weight_images_2D_[index];
}

template <class TPixel>
void 
fregl_space_transformer< TPixel >::
write_xml(std::string const & montage_xml, std::string const & montage_directory, 
		std::string const & montage_2d_name, bool overlap_only, 
		bool in_anchor, int channel, int blending, bool use_nn, bool denoised)
{
	TiXmlDocument doc;
	std::string str;

	// TODO: record the channel number

	/* 
	 * Creates a new document, a node and set it as a root node
	 */

	TiXmlElement* root_node = new TiXmlElement( "Montaging" );
	str = ToString( image_id_indices_.size() );
	doc.LinkEndChild( root_node );
	root_node->SetAttribute("number_of_images", str.c_str());
	root_node->SetAttribute("montage_directory", montage_directory.c_str());
	root_node->SetAttribute("montage_2D_proj_image", montage_2d_name.c_str());

	str = joint_register_->image_name(anchor_);
	root_node->SetAttribute("anchor", str.c_str());

	str = ToString( channel );
	root_node->SetAttribute("channel", str.c_str());

	if (overlap_only) 
		root_node->SetAttribute("overlapping_images_only", "yes");
	else
		root_node->SetAttribute("overlapping_images_only", "no");

	if (in_anchor) 
		root_node->SetAttribute("in_anchor_space", "yes");
	else
		root_node->SetAttribute("in_anchor_space", "no");

	switch (blending) {
		case 0: 
			root_node->SetAttribute("blending", "maximum");
			break;
		case 1: 
			root_node->SetAttribute("blending", "averaging");
			break;
		case 2: 
			root_node->SetAttribute("blending", "photo-bleaching-weighted");
			break;
		default:
			root_node->SetAttribute("blending", "unknown");
	}

	if (use_nn) 
		root_node->SetAttribute("nearest_neighbor_interpolation", "yes");
	else
		root_node->SetAttribute("nearest_neighbor_interpolation", "no");

	if (denoised) 
		root_node->SetAttribute("denoised", "yes");
	else
		root_node->SetAttribute("denoised", "no");

	// origin
	TiXmlElement* node = new TiXmlElement("origin");
	root_node->LinkEndChild(node);
	str = ToString(origin_[0]);    
	TiXmlElement* sub_node = new TiXmlElement("x");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);
	str = ToString(origin_[1]);    
	sub_node = new TiXmlElement("y");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);
	str = ToString(origin_[2]);    
	sub_node = new TiXmlElement("z");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);

	//size
	node = new TiXmlElement("size");
	root_node->LinkEndChild(node);
	str = ToString(image_size_[0]);    
	sub_node = new TiXmlElement("x");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);
	str = ToString(image_size_[1]);    
	sub_node = new TiXmlElement("y");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);
	str = ToString(image_size_[2]);    
	sub_node = new TiXmlElement("z");
	sub_node->LinkEndChild( new TiXmlText( str.c_str()) );
	node->LinkEndChild(sub_node);

	// All transformations from the anchor to other images involved in
	// this montaging
	for (unsigned int i = 0; i<image_id_indices_.size(); i++) {
		fregl_reg_record::Pointer reg_rec = joint_register_->get_reg_record(anchor_,image_id_indices_[i]);
		reg_rec->write_xml_node(root_node);
	}    

	/* 
	 * Dumping document to a file
	 */
	doc.SaveFile( montage_xml.c_str() );

}

template <class TPixel>
void 
fregl_space_transformer< TPixel >::
read_xml(std::string const & filename, std::string& montage_directory, 
		std::string& montage_2d_name)
{
	TiXmlDocument doc; /* the resulting document tree */

	//Parse the resource
	if ( !doc.LoadFile( filename.c_str() ) ) {
		vcl_cout<<"Unable to load XML File"<<vcl_endl;
		return ;
	}

	/*Get the root element node */
	TiXmlElement* root_element = doc.FirstChildElement();
	const char* contents = root_element->Value();
	if ( strcmp( contents, "Montaging" ) != 0 ) {
		vcl_cout<<"Incorrect XML root Element: "<<vcl_endl;
		return ;
	}

	//Read the basic information
	//
	montage_directory = root_element->Attribute("montage_directory");
	montage_2d_name = root_element->Attribute("montage_2D_proj_image");
	std::string anchor_image_name = root_element->Attribute("anchor");
	std::string in_anchor_str = root_element->Attribute("in_anchor_space");
	//int channel = atoi(root_element->Attribute("channel")); // not used for now

	bool in_anchor;
	if ( in_anchor_str == "yes") in_anchor = true;
	else in_anchor = false;

	std::string overlap_only_str = root_element->Attribute("overlapping_images_only");

	bool overlap_only;
	if ( overlap_only_str == "yes") overlap_only = true;
	else overlap_only = false;

	//Read the origin, size & transformations
	std::vector<fregl_reg_record::Pointer> reg_records;
	TiXmlElement* cur_node =  root_element->FirstChildElement();
	for ( ; cur_node; cur_node = cur_node->NextSiblingElement() ) {
		const char * value = cur_node->Value();

		//origin
		if ( strcmp(value, "origin") == 0 ) {
			TiXmlElement* sub_node = cur_node->FirstChildElement();
			for (; sub_node; sub_node = sub_node->NextSiblingElement()) {
				if ( strcmp("x", sub_node->Value()) == 0 ) {
					origin_[0] = atoi( sub_node->GetText() );
					continue;
				}
				if ( strcmp("y", sub_node->Value()) == 0 ) {
					origin_[1] = atoi( sub_node->GetText() );
					continue;
				}
				if ( strcmp("z", sub_node->Value()) == 0 ) {
					origin_[2] = atoi( sub_node->GetText() );
					continue;
				}
			}
			continue;
		}

		//size
		if ( strcmp(value, "size") == 0 ) {
			TiXmlElement* sub_node = cur_node->FirstChildElement();
			for (; sub_node; sub_node = sub_node->NextSiblingElement()) {
				if ( strcmp("x", sub_node->Value()) == 0 ) {
					image_size_[0] = atoi( sub_node->GetText() );
					continue;
				}
				if ( strcmp("y", sub_node->Value()) == 0 ) {
					image_size_[1] = atoi( sub_node->GetText() );
					continue;
				}
				if ( strcmp("z", sub_node->Value()) == 0 ) {
					image_size_[2] = atoi( sub_node->GetText() );
					continue;
				}
			}
			continue;
		}

		//transformation
		if ( !strcmp(value, "Transform") ) {
			fregl_reg_record::Pointer reg_rec = new fregl_reg_record();
			reg_rec->read_xml_node(cur_node);
			reg_records.push_back( reg_rec );
			continue;
		}
	}
	joint_register_ = new fregl_joint_register< TPixel >(reg_records);

	//Set anchor_ and image_id_indices_
	this->set_anchor(anchor_image_name, in_anchor, overlap_only, true);
}

static 
	std::string 
ToString(double val)
{
	std::ostringstream strm;
	strm<< val<<std::endl;
	return strm.str();
}

template <class TPixel>
typename fregl_space_transformer< TPixel >::ImageType2DPointer
fregl_space_transformer< TPixel >::
max_projection(ImageTypePointer image, float sigma) const
{
	typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
	typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
	typedef itk::DiscreteGaussianImageFilter< ImageType2D,FloatImageType2D > SmoothingFilterType;
	typedef itk::CastImageFilter< FloatImageType2D, ImageType2D > CastFilterType;

	typename ImageType2D::RegionType region;
	typename ImageType2D::RegionType::SizeType size;
	typename ImageType2D::RegionType::IndexType index;
	typename ImageType::RegionType requestedRegion = image->GetRequestedRegion();
	index[ 0 ] = 0;
	index[ 1 ] = 0;
	size[ 0 ] = requestedRegion.GetSize()[ 0 ];
	size[ 1 ] = requestedRegion.GetSize()[ 1 ];
	region.SetSize( size );
	region.SetIndex( index );
	ImageType2DPointer image2D = ImageType2D::New();
	image2D->SetRegions( region );
	image2D->Allocate();

	//Set the iterator
	SliceIteratorType output3DIt( image, image->GetRequestedRegion() );
	LinearIteratorType output2DIt( image2D, image2D->GetRequestedRegion() );

	unsigned int direction[2];
	direction[0] = 0;
	direction[1] = 1;

	output3DIt.SetFirstDirection(direction[1]);
	output3DIt.SetSecondDirection(direction[0]);
	output2DIt.SetDirection(1 - direction[0]);

	// Initialized the 2D image
	output2DIt.GoToBegin();
	while ( ! output2DIt.IsAtEnd() ) 
	{
		while ( ! output2DIt.IsAtEndOfLine() ) 
		{
			output2DIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
			++output2DIt;
		}
		output2DIt.NextLine();
	}

	// Now do the max projection, 
	output3DIt.GoToBegin();
	output2DIt.GoToBegin();

	while( !output3DIt.IsAtEnd() ) {
		while ( !output3DIt.IsAtEndOfSlice() ) {
			while ( !output3DIt.IsAtEndOfLine() ) {
				output2DIt.Set( vnl_math_max( output3DIt.Get(), output2DIt.Get() ));
				++output3DIt;
				++output2DIt;
			}
			output2DIt.NextLine();
			output3DIt.NextLine();
		}
		output2DIt.GoToBegin();
		output3DIt.NextSlice();
	}

	if (sigma <= 0) return image2D;

	// Perform Gaussian smoothing if the 
	typename SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
	typename CastFilterType::Pointer caster = CastFilterType::New();
	smoother->SetInput( image2D );
	smoother->SetVariance(sigma);
	smoother->SetMaximumKernelWidth(15);
	caster->SetInput( smoother->GetOutput() );
	try {
		caster->Update();
	}
	catch(itk::ExceptionObject& e) {
		vcl_cout << e << vcl_endl;
	}

	return caster->GetOutput();
}

//Explicit Instantiation
template class fregl_space_transformer< unsigned char >;
template class fregl_space_transformer< unsigned short >;
