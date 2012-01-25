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
// \brief Utility functions 
// \author Charlene Tsai
// \date 04/14/2008
// 
#ifndef _fregl_util_
#define _fregl_util_

#include "itkImage.h"
#include "vil3d/vil3d_image_view.h"
#include "itkRGBPixel.h"
#include "itkAffineTransform.h"

typedef unsigned char                    InputPixelType;
typedef itk::Image< InputPixelType, 3 >  ImageType;
typedef itk::Image< InputPixelType, 2 >  ImageType2D;
typedef itk::Image< float, 2 >           FloatImageType2D;
typedef itk::RGBPixel< unsigned char >   ColorPixelType;
typedef itk::Image< ColorPixelType, 3 > ColorImageType;
typedef itk::Image< ColorPixelType, 2 > ColorImageType2D;
typedef itk::AffineTransform< double, 3>   TransformType;

// Maximum is taken between the two images
ImageType::Pointer fregl_util_fuse_images(ImageType::Pointer image1, ImageType::Pointer image2);

ImageType::Pointer fregl_util_read_image( std::string const & file_name, bool channel_set = false, int channel = 0, bool denoise = false);

ImageType2D::Pointer fregl_util_max_projection(ImageType::Pointer image, float sigma = 0);

//: Computer the ammount of overlap between two images given the transformation
double fregl_util_overlap(TransformType::Pointer transform, itk::Size<3> size_from, itk::Size<3> size_to);

ColorImageType2D::Pointer fregl_util_max_projection_color(ColorImageType::Pointer image);

//: High frequency image (LoG) is removed from the original image
void fregl_util_reduce_noise(ImageType::Pointer image);

//: convert vil_image_view to itk::Image 
//
//  The code is hacked from ril_itk_convert.h in rpi package
ImageType::Pointer fregl_util_convert_vil_to_itk( vil3d_image_view<InputPixelType> img );

//: extract single channel from lsm image
//ImageType::Pointer fregl_util_lsm_one_channel(std::string filename, int channel);
#endif
