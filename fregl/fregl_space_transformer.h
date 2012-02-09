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

//: 
// \file
// \brief The class transforms the individual images to a common space
// \author Charlene Tsai
// \date 12/05/2007

#ifndef _fregl_space_transformer_h_
#define _fregl_space_transformer_h_

#include <fregl/fregl_util.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vector>

#include "itkIndex.h"
#include "itkSize.h"
#include "itkPoint.h"
#include "itkImage.h"

#include <fregl/fregl_joint_register.h>
// The space transformer assumes [1,1,1] for voxel spacing. 
class fregl_space_transformer: public vbl_ref_count
{
public:
	typedef vbl_smart_ptr< fregl_space_transformer >  Pointer;
	typedef fregl_util::InputPixelType			InputPixelType;
	typedef fregl_joint_register::TransformType TransformType;
	typedef itk::Image<InputPixelType, 3>        ImageType;
	typedef itk::Image<InputPixelType, 2>        ImageType2D;
	typedef ImageType::PointType                PointType; //physical space
	typedef ImageType::SizeType                 SizeType;
	typedef ImageType::SpacingType              SpacingType;
	typedef ImageType::IndexType                IndexType;
	typedef itk::Image< float, 2 >              FloatImageType2D;

	fregl_space_transformer();
	fregl_space_transformer( fregl_joint_register::Pointer joint_reg );
	~fregl_space_transformer(){}

	//: Set the new anchor image.
	//
	//  xforms_, inverse_xforms_, and anchor_ are updated accordingly
	void set_anchor(std::string const & anchor_name, bool in_anchor = false, bool overlap_only = false, bool space_set = false);
	
	//:  Set the roi for processing
	//   "s_origin" is the new roi origin point
	//   "s_size" is the new roi size
	void set_roi(PointType s_origin, SizeType s_size);

	//: Determine if a point is in the image of given index
	//
	//  "loc" is the location of the point in the anchor image, and the
	//  "xformed_loc" is the transformed location in the image of the
	//  given idnex. The transformed location is always returned, even
	//  if the transformed point is not in the image range.
	bool in_image(PointType loc, int image_index, 
		PointType& xformed_loc) const;

	//: Determine if a point is in the image of given index
	//
	//  "loc" is the location of the point in the anchor image, and the
	//  "xformed_loc" is the transformed location in the image of the
	//  given idnex. The transformed location is always returned, even
	//  if the transformed point is not in the image range.
	bool in_image(vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 3 >& xformed_loc) const;

	//: Determine if an image is in the roi
	//
	//  "image_index" is the joint register index for the image
	//  true is returned if any part of the image is in the roi

	bool image_in_roi(int image_index) const;

	//: Determine if a point is in the 2D projected image of given index
	//
	//  "loc" is the location of the point in the anchor image, and the
	//  "xformed_loc" is the transformed location in the image of the
	//  given idnex. The transformed location is always returned, even
	//  if the transformed point is not in the image range.
	bool in_image_2d(vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 2 >& xformed_loc) const;

	//: Determine if a point in the image of given index is in the anchor
	//
	//  "loc" is the location of the point in the image of given index,
	//  and the "xformed_loc" is the transformed location in the anchor
	//  image. The transformed location is always returned, even if the
	//  transformed point is not in the bounding box of the anchor
	//  image.
	bool in_anchor(PointType loc, int image_index, PointType& xformed_loc) const;

	//: Determine if a point in the image of given index is in the anchor
	//
	//  "loc" is the location of the point in the image of given index,
	//  and the "xformed_loc" is the transformed location in the anchor
	//  image. The transformed location is always returned, even if the
	//  transformed point is not in the bounding box of the anchor
	//  image.
	bool in_anchor( vnl_vector_fixed< float, 3 > loc, int image_index, 
		vnl_vector_fixed< float, 3 >& xformed_loc) const;

	//: Determine if a point is in the "to" image
	//
	//  "loc" is the location of the point in the "from" image, and the
	//  "xformed_loc" is the transformed location in "to" image of the
	//  given idnex. The transformed location is always returned, even
	//  if the transformed point is not in the image range.
	bool in_range(PointType loc, int from_index, int to_index, 
		PointType& xformed_loc) const;

	//: Determine if a point is in the "to" image
	//
	//  "loc" is the location of the point in the "from" image, and the
	//  "xformed_loc" is the transformed location in the "to" image of
	//  the given idnex. The transformed location is always returned,
	//  even if the transformed point is not in the image range.
	bool in_range(vnl_vector_fixed< float, 3 > loc, int from_index, int to_index, 
		vnl_vector_fixed< float, 3 >& xformed_loc) const;

	//: Generate the transformed image using the given transformation
	//
	//  The given image is transformed to the image space of the global
	//  space defined by the anchor image.
	ImageType::Pointer transform_image(ImageType::Pointer in_image, int image_index, int background = 0, bool use_NN_interpolator = false) const;

	//: Generate the transformed image using the given transformation
	//
	//  The given image is transformed to the image space of the global
	//  space defined by the anchor image and the roi.
	ImageType::Pointer transform_image_roi(ImageType::Pointer in_image, int image_index, int background = 0, bool use_NN_interpolator = false) const;

        //: Generate the transformed image using the given transformation
	//
	//  The given image is transformed to the image space of the global
	//  space defined by the anchor image only the size of the actual image read.
	ImageType::Pointer transform_image_whole(ImageType::Pointer in_image, int image_index, int background = 0, bool use_NN_interpolator = false) const;

	//: Generate the photo-bleaching weighted transformed image using the given transformation
	//
	//  The given image is transformed to the image space of the global
	//  space defined by the anchor image.
	ImageType::Pointer transform_image_weighted(ImageType::Pointer in_image, int image_index, int background = 0, bool use_NN_interpolator = false) const;

	//: Generate the image which keeps track of the number of images overlapping at each pixel
	//
	//  The function returns an image of the same dimension as the final image. 
	ImageType::Pointer compute_weighted_image() const;

	//: Generate the image which keeps track of the number of images overlapping at each pixel
	//
	//  The function returns a 2D image of the same x-y dimension as the
	//  final image. The 2D image is taken as the middle slice of the
	//  global volume. It is a good approximation since images sharing
	//  sections contain cells which show up in both images.
	ImageType2D::Pointer compute_weighted_image_2D() const;

	//: Return the names of the set of images 
	std::vector<std::string> image_names() const;

	//: Return the sizes of the set of images 
	std::vector<SizeType> image_sizes() const;

	//: Return the origin of the global space
	//
	//  If on the anchor space is considered, the origin is [0,0,0],
	//  else it is the minimum of the global bounding box
	PointType origin() const { return origin_; }

	//: Return the global image size; 
	SizeType montage_size() const { return image_size_; }

	//: Return the list of xforms from the given image to its neighbors
	//
	//  Neighbors are a set of overlapping images
	std::vector<TransformType::Pointer> xforms_to_neighbors(int image_index) const;

	//: Return the list of xforms to the given image from its neighbors
	//
	//  Neighbors are a set of overlapping images
	std::vector<TransformType::Pointer> xforms_from_neighbors(int image_index) const;

	//: Return the list of xforms from the given image to all other images
	//
	//  The other images are from the set of images in consideration. If
	//  "overlap_only" is set, then the other images are only those
	//  overlap with the anchor image.
	std::vector<TransformType::Pointer> xforms_to_all(int image_index) const;

	//: Return the list of xforms to the given image from all other images
	//
	//  The other images are from the set of images in consideration. If
	//  "overlap_only" is set, then the other images are only those
	//  overlap with the anchor image.
	std::vector<TransformType::Pointer> xforms_from_all(int image_index) const;

	//: Set 2D weight map for the given image
	//
	// This is handle different weighting due to photobleaching
	void set_individual_weight_map(int index, ImageType::Pointer image, float alpha);

	//: Normalized the weight map by the weighted sum
	//
	//  This function should be called only if all weight maps are updated. 
	void normalize_individual_weight_maps();

	//: Return the normalized weight map of given index
	ImageType2D::Pointer get_weight_map(int index) const;

	// IO with xml file
	void write_xml(std::string const & montage_xml, 
		std::string const & montage_directory,
		std::string const & montage_2d_name,
		bool overlap_only, bool in_anchor, int channel = 0,
		int blending=0, bool use_nn=false, bool denoised=false);
	void read_xml(std::string const & filename, std::string& montage_directory, 
		std::string& montage_2d_name);

private:
	ImageType2D::Pointer max_projection(ImageType::Pointer image, float sigma = 0) const;

private: 
	fregl_joint_register::Pointer joint_register_;
	int anchor_;
	//std::vector<TransformType::Pointer> inverse_xforms_; //taking anchor to other spaces
	std::vector<int> image_id_indices_;
	std::vector<ImageType2D::Pointer> weight_images_2D_; //same size as image_id_indices
	std::vector<ImageType2D::Pointer> normalized_weight_images_2D_;
	std::vector<PointType> image_origins_;

	// The physical space is defined by the reference image
	bool anchor_set_;
	PointType origin_;
	SizeType  image_size_;
	PointType roi_origin_;
	SizeType  roi_size_;
	SpacingType spacing_;
};
#endif
