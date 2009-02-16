//: 
// \file
// \brief Utility functions 
// \author Charlene Tsai
// \date 04/14/2008
// 
#ifndef _fregl_util_
#define _fregl_util_

#include "itkImage.h"
#include <vil3d/vil3d_image_view.h>

typedef unsigned char                    InputPixelType;
typedef itk::Image< InputPixelType, 3 >  ImageType;
typedef itk::Image< InputPixelType, 2 >  ImageType2D;
typedef itk::Image< float, 2 >           FloatImageType2D;

// Maximum is taken between the two images
ImageType::Pointer fregl_util_fuse_images(ImageType::Pointer image1, ImageType::Pointer image2);

ImageType::Pointer fregl_util_read_image( std::string const & file_name, bool channel_set = false, int channel = 0, bool denoise = false);

ImageType2D::Pointer fregl_util_max_projection(ImageType::Pointer image, float sigma = 0);

//: High frequency image (LoG) is removed from the original image
void fregl_util_reduce_noise(ImageType::Pointer image);

//: convert vil_image_view to itk::Image 
//
//  The code is hacked from ril_itk_convert.h in rpi package
ImageType::Pointer fregl_util_convert_vil_to_itk( vil3d_image_view<InputPixelType> img );

//: extract single channel from lsm image
//ImageType::Pointer fregl_util_lsm_one_channel(std::string filename, int channel);
#endif
