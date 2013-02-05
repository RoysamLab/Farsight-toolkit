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

#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkTIFFImageIO.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <vtkImageData.h>
#include <ftkImage/ftkImage.h>

template < class TPixel >
class fregl_util
{
public: 
	typedef TPixel							InputPixelType;
	typedef TPixel							GDBICPPixelType;

	typedef itk::Image< InputPixelType, 3 >		ImageType;
	typedef typename ImageType::Pointer		ImageTypePointer;
	typedef itk::Image< InputPixelType, 2 >		ImageType2D;
	typedef typename ImageType2D::Pointer		ImageType2DPointer;

	typedef itk::Image< GDBICPPixelType, 2 >	GDBICPImageType;
	typedef itk::Image< float, 2 >			FloatImageType2D;
	typedef typename FloatImageType2D::Pointer		FloatImageType2DPointer;
	typedef itk::Image< float, 3 >			FloatImageType;
	typedef typename FloatImageType::Pointer		FloatImageTypePointer;
	typedef itk::RGBPixel< unsigned char >		ColorPixelType;
	typedef itk::Image< ColorPixelType, 3 >		ColorImageType;
	typedef itk::Image< ColorPixelType, 2 >		ColorImageType2D;
	typedef itk::AffineTransform< double, 3>	TransformType;

	typedef typename itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
	typedef typename itk::ImageRegionIterator< ImageType > RegionIterator;

	// Maximum is taken between the two images
	static ImageTypePointer fregl_util_fuse_images(ImageTypePointer image1, ImageTypePointer image2);

	static ImageTypePointer fregl_util_read_image( std::string const & file_name, bool channel_set = false, int channel = 0, bool denoise = false);

	static ImageType2DPointer fregl_util_max_projection(ImageTypePointer image, float sigma = 0);

	static FloatImageType2DPointer fregl_util_min_projection(FloatImageTypePointer image, float sigma = 0);

	static ImageType2DPointer fregl_util_fast_max_projection(ImageTypePointer image);

	//: Computer the ammount of overlap between two images given the transformation
	static double fregl_util_overlap(TransformType::Pointer transform, itk::Size<3> size_from, itk::Size<3> size_to);

	static ColorImageType2D::Pointer fregl_util_max_projection_color(ColorImageType::Pointer image);

	//: High frequency image (LoG) is removed from the original image
	static void fregl_util_reduce_noise(ImageTypePointer image);

	//: convert vil_image_view to itk::Image 
	//
	//  The code is hacked from ril_itk_convert.h in rpi package
	static ImageTypePointer fregl_util_convert_vil_to_itk( vil3d_image_view<InputPixelType> img );

	//: extract single channel from lsm image
	//ImageTypePointer fregl_util_lsm_one_channel(std::string filename, int channel);
};
#endif
