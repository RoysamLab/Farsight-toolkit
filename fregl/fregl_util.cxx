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

#include <fregl/fregl_util.h>
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

typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
typedef itk::ImageRegionIterator< ImageType > RegionIterator;

ImageType::Pointer
fregl_util_fuse_images(ImageType::Pointer image1, ImageType::Pointer image2)
{
	ImageType::Pointer image_out = ImageType::New();
	image_out->SetRegions(image1->GetRequestedRegion());
	image_out->Allocate();

	//Set the iterator
	RegionConstIterator inputIt1( image1, image1->GetRequestedRegion() );
	RegionConstIterator inputIt2( image2, image2->GetRequestedRegion() );
	RegionIterator outputIt( image_out, image_out->GetRequestedRegion() );

	for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); 
		!inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++outputIt)
	{
		outputIt.Set( vnl_math_max(inputIt1.Get(), inputIt2.Get()) );
	}

	return image_out;
}

ImageType::Pointer
fregl_util_read_image( std::string const & file_name, bool channel_set, int channel, bool denoise )
{
	// Read in the image using ftkImage, which serves as a reader. It is
	// capable of reading images of different format and dealing with
	// multiple channels, whether data or color channels.

	ftk::Image::Pointer imageReader = ftk::Image::New();
	imageReader->LoadFile( file_name );
	ImageType::Pointer final_image = ImageType::New();

	// get the extension
	std::string ext;
	const std::string dot = ".";
	std::string::size_type po = file_name.find_last_of(dot);
	ext = file_name.substr(po+1,file_name.length()-1);
	const ftk::Image::Info* info = imageReader->GetImageInfo();
	int numChannels = info->numChannels;

	// Have to give special care to lsm images, since channels are
	// extracted and fused for registration if channel_set is false
	if ( !channel_set && (ext=="lsm" || ext=="LSM") ) {  
		ImageType::Pointer image = imageReader->GetItkPtr<InputPixelType>(0, 0);
		final_image->SetRegions( image->GetBufferedRegion() );
		final_image->Allocate();
		final_image->FillBuffer(0);
		for (int counter=0; counter < numChannels; counter++) {
			image = imageReader->GetItkPtr<InputPixelType>(0, counter, ftk::Image::RELEASE_CONTROL);
			final_image = fregl_util_fuse_images(final_image, image);
		}
	}
	else { // return one single channel
		std::cout<<"Number of channels = "<<numChannels<<std::endl;
		std::cout<<"Channel "<<channel<<" extracted"<<std::endl;
		std::cout<<"bytePerPix = "<<int(info->bytesPerPix)<<std::endl;
		ImageType::Pointer image = imageReader->GetItkPtr<InputPixelType>(0, channel,ftk::Image::RELEASE_CONTROL);
		final_image->SetRegions( image->GetBufferedRegion() );
		final_image->Allocate();
		final_image->FillBuffer(0);
		// What it does here is simply to copy the content of image to
		// final_image
		final_image = fregl_util_fuse_images(final_image, image);
	}

	// Smooth the image if needed
	if (denoise) fregl_util_reduce_noise( final_image );

	return final_image;
}

ImageType2D::Pointer
fregl_util_max_projection(ImageType::Pointer image, float sigma)
{
	typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
	typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
	typedef itk::DiscreteGaussianImageFilter< ImageType2D,FloatImageType2D > SmoothingFilterType;
	typedef itk::CastImageFilter< FloatImageType2D, ImageType2D > CastFilterType;

	ImageType2D::RegionType region;
	ImageType2D::RegionType::SizeType size;
	ImageType2D::RegionType::IndexType index;
	ImageType::RegionType requestedRegion = image->GetRequestedRegion();
	index[ 0 ] = 0;
	index[ 1 ] = 0;
	size[ 0 ] = requestedRegion.GetSize()[ 0 ];
	size[ 1 ] = requestedRegion.GetSize()[ 1 ];
	region.SetSize( size );
	region.SetIndex( index );
	ImageType2D::Pointer image2D = ImageType2D::New();
	image2D->SetRegions( region );
	image2D->Allocate();

	//Set the iterator
	SliceIteratorType output3DIt( image, image->GetRequestedRegion() );
	LinearIteratorType output2DIt( image2D, image2D->GetRequestedRegion() );

	unsigned int direction[2];
	direction[0] = 0;
	direction[1] = 1;

	output3DIt.SetFirstDirection( direction[1] );
	output3DIt.SetSecondDirection( direction[0] );
	output2DIt.SetDirection( 1 - direction[0] );

	// Initialized the 2D image
	output2DIt.GoToBegin();
	while ( ! output2DIt.IsAtEnd() ) {
		while ( ! output2DIt.IsAtEndOfLine() ) {
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
	SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
	CastFilterType::Pointer caster = CastFilterType::New();
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

//: Determine the approximate amount of overlap between two images
//
//  Whether two images overlap is roughly determined by the
//  translation in x and y, since we assume not much distortion from
//  other parameters.
double
fregl_util_overlap(TransformType::Pointer transform, itk::Size<3> size_from, itk::Size<3> size_to)
{
	bool is_overlap = true;
	double overlap = 0.0;

	if ( !transform ) return 0;

	// A quick check to see if the two volumes overlap
	TransformType::ParametersType params = transform->GetParameters();
	double tx = params[9];
	double ty = params[10];

	if (tx > size_to[0] || ty > size_to[1]) is_overlap = false;
	if (-tx> size_from[0] || -ty> size_from[1]) is_overlap = false;

	if (!is_overlap) 
		return overlap;

	//Work out the appoximate amount of overlapping using only the x-y
	//translation since the z-dimension changes.
	double total_area = size_to[0]*size_to[1];
	double top_x = tx>0?tx:0;
	double top_y = ty>0?ty:0;
	double bot_x = size_to[0]>size_from[0]+tx?size_from[0]-1+tx:size_to[0]-1;
	double bot_y = size_to[1]>size_from[1]+ty?size_from[1]-1+ty:size_to[1]-1;
	overlap = (bot_x-top_x)*(bot_y-top_y)/total_area;

	return overlap;
}

ColorImageType2D::Pointer
fregl_util_max_projection_color(ColorImageType::Pointer image)
{
	typedef itk::ImageLinearIteratorWithIndex< ColorImageType2D > LinearIteratorType;
	typedef itk::ImageSliceIteratorWithIndex< ColorImageType > SliceIteratorType;
	//typedef itk::DiscreteGaussianImageFilter< ColorImageType2D,FloatImageType2D > SmoothingFilterType;
	// typedef itk::CastImageFilter< FloatImageType2D, ImageType2D > CastFilterType;

	ColorImageType2D::RegionType region;
	ColorImageType2D::RegionType::SizeType size;
	ColorImageType2D::RegionType::IndexType index;
	ColorImageType::RegionType requestedRegion = image->GetRequestedRegion();
	index[ 0 ] = 0;
	index[ 1 ] = 0;
	size[ 0 ] = requestedRegion.GetSize()[ 0 ];
	size[ 1 ] = requestedRegion.GetSize()[ 1 ];
	region.SetSize( size );
	region.SetIndex( index );
	ColorImageType2D::Pointer image2D = ColorImageType2D::New();
	image2D->SetRegions( region );
	image2D->Allocate();
	image2D->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));

	//Set the iterator
	SliceIteratorType output3DIt( image, image->GetRequestedRegion() );
	LinearIteratorType output2DIt( image2D, image2D->GetRequestedRegion() );

	unsigned int direction[2];
	direction[0] = 0;
	direction[1] = 1;

	output3DIt.SetFirstDirection( direction[1] );
	output3DIt.SetSecondDirection( direction[0] );
	output2DIt.SetDirection( 1 - direction[0] );
	/*
	// Initialized the 2D image
	output2DIt.GoToBegin();
	while ( ! output2DIt.IsAtEnd() ) {
	while ( ! output2DIt.IsAtEndOfLine() ) {
	output2DIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
	++output2DIt;
	}
	output2DIt.NextLine();
	}
	*/

	// Now do the max projection, 
	output3DIt.GoToBegin();
	output2DIt.GoToBegin();

	while( !output3DIt.IsAtEnd() ) {
		while ( !output3DIt.IsAtEndOfSlice() ) {
			while ( !output3DIt.IsAtEndOfLine() ) {
				//output2DIt.Set( vnl_math_max( output3DIt.Get(), output2DIt.Get() ));
				ColorPixelType pix = output2DIt.Get();
				pix.SetRed(vnl_math_max( output3DIt.Get().GetRed(), output2DIt.Get().GetRed()));
				pix.SetGreen(vnl_math_max( output3DIt.Get().GetGreen(), output2DIt.Get().GetGreen()));
				pix.SetBlue(vnl_math_max( output3DIt.Get().GetBlue(), output2DIt.Get().GetBlue()));
				output2DIt.Set( pix );
				++output3DIt;
				++output2DIt;
			}
			output2DIt.NextLine();
			output3DIt.NextLine();
		}
		output2DIt.GoToBegin();
		output3DIt.NextSlice();
	}
	return image2D;

	/*
	if (sigma <= 0) return image2D;

	// Perform Gaussian smoothing if the 
	SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
	CastFilterType::Pointer caster = CastFilterType::New();
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
	*/
}

void
fregl_util_reduce_noise(ImageType::Pointer image)
{
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType2D, FloatImageType2D > LoGFilterType2D;
	typedef itk::SubtractImageFilter< ImageType2D, FloatImageType2D, FloatImageType2D> SubtractFilterType;
	typedef itk::RescaleIntensityImageFilter< FloatImageType2D, ImageType2D> RescaleFilterType;
	typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
	typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;

	ImageType2D::Pointer image_2D = fregl_util_max_projection(image,0);

	// run LoG on the 2D image to define the max intensity for each pixel.
	//
	LoGFilterType2D::Pointer laplacian = LoGFilterType2D::New();
	SubtractFilterType::Pointer subtracter = SubtractFilterType::New();
	RescaleFilterType::Pointer rescale = RescaleFilterType::New();
	laplacian->SetInput( image_2D );
	laplacian->SetSigma( 2 );
	try {
		laplacian->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cout << "ExceptionObject: "<< err << std::endl;
		return;
	}
	FloatImageType2D::Pointer image_high=laplacian->GetOutput();
	image_high->DisconnectPipeline();

	subtracter->SetInput1( image_2D );
	subtracter->SetInput2( image_high );
	rescale->SetInput( subtracter->GetOutput() );
	rescale->SetOutputMinimum(   0 );
	rescale->SetOutputMaximum( 255 );
	try {
		rescale->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cout << "ExceptionObject: "<< err << std::endl;
		return;
	}
	ImageType2D::Pointer image_low_2D = rescale->GetOutput();

	/* It is not good to apply any sort of smoothing.
	// Perform anisotropic diffusion image filtering to smooth the image
	typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType2D, 
	FloatImageType2D > ADFilterType;
	ADFilterType::Pointer adfilter = ADFilterType::New();
	adfilter->SetInput( rescale->GetOutput() );       
	adfilter->SetNumberOfIterations( 5 );
	adfilter->SetTimeStep( 0.125 );
	adfilter->SetConductanceParameter( 3.0 );
	rescale->SetInput(adfilter->GetOutput());
	try {
	rescale->Update();
	}
	catch( itk::ExceptionObject & err ) {
	std::cout << "ExceptionObject: "<< err << std::endl;
	return;
	}               
	*/

	// scale down each pixel using the value in image as using the ratio
	// of low_max/org_max where low_max is the max in along the
	// z-direction in image_low, and org_max for the original image.
	//
	SliceIteratorType output3DIt( image, image->GetRequestedRegion() );
	LinearIteratorType output2DIt( image_2D, image_2D->GetRequestedRegion() );
	LinearIteratorType output2DIt_low( image_low_2D, image_low_2D->GetRequestedRegion() );

	unsigned int direction[2];
	direction[0] = 0;
	direction[1] = 1;

	output3DIt.SetFirstDirection( direction[1] );
	output3DIt.SetSecondDirection( direction[0] );
	output2DIt.SetDirection( 1 - direction[0] );
	output2DIt_low.SetDirection( 1 - direction[0] );
	output3DIt.GoToBegin();
	output2DIt.GoToBegin();
	output2DIt_low.GoToBegin();
	while( !output3DIt.IsAtEnd() ) {
		while ( !output3DIt.IsAtEndOfSlice() ) {
			while ( !output3DIt.IsAtEndOfLine() ) {
				if (output2DIt.Get() > 0)
					output3DIt.Set(int(output2DIt_low.Get()/(float)output2DIt.Get() * output3DIt.Get() ));
				++output3DIt;
				++output2DIt;
				++output2DIt_low;
			}
			output2DIt.NextLine();
			output2DIt_low.NextLine();
			output3DIt.NextLine();
		}
		output2DIt.GoToBegin();
		output2DIt_low.GoToBegin();
		output3DIt.NextSlice();
	}
}

ImageType::Pointer
fregl_util_convert_vil_to_itk( vil3d_image_view<InputPixelType> img )
{
	// Then the image object can be created
	ImageType::Pointer image = ImageType::New();

	if( !img.is_contiguous() ) {
		vcl_cerr <<"Cannot convert a non-contiguous vil image"<< vcl_endl;
		std::cout<<"ROI dimension"<<img.ni()<<","<<img.nj()<<","<<img.nk()<<std::endl;
		return image;
	}

	ImageType::IndexType start;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::SizeType  size;
	size[0]  = img.ni();  // size along X
	size[1]  = img.nj();  // size along Y
	size[2]  = img.nk(); 

	ImageType::RegionType region;
	region.SetSize( size );  
	region.SetIndex( start );

	image->SetRegions( region );
	// image->SetBufferedRegion( region );

	image->GetPixelContainer()->SetImportPointer( img.origin_ptr(),
		img.size_bytes(),
		false );
	return image;
}

/*
ImageType::Pointer
fregl_util_lsm_one_channel(std::string filename, int channel)
{
vtkLSMReader * lsmR = vtkLSMReader::New();
lsmR->SetFileName(file_name.c_str());
lsmR->OpenFile();
int numChannels = lsmR->GetNumberOfChannels();
int numTimes = lsmR->GetNumberOfTimePoints();
std::cout<<"numChannels = "<<numChannels<<std::endl;
std::cout<<"numTimes = "<<numTimes<<std::endl;
if (numChannels == 0 || numTimes > 1){
std::cerr <<"Cannot handle the current format."<<std::endl;
return 0;
}
vtkImageData * vimdata;
int imsize, size_x, size_y, size_z;
int extent[6]; //recording dimensions of the images. extent[0,2,4]
//is the min location, and extent[1,3,5] is the max
//location in 3D
int time_stamp = 0;
vimdata = lsmR->GetTimePointOutput(time_stamp, channel);
lsmR->Update();
vimdata->GetExtent(extent);
size_x = extent[1]-extent[0]+1;
size_y = extent[3]-extent[2]+1;
size_z = extent[5]-extent[4]+1;
imsize = (size_z)*(size_y)*(size_x);

unsigned char *channel_data = (unsigned char*) vimdata->GetScalarPointer();
ImageType::Pointer image = ImageType::New();
ImageType::IndexType start;
start[0] = 0; // first index on X
start[1] = 0; // first index on Y
start[2] = 0; // first index on Z

ImageType::SizeType size;
size[0] = size_x; // size along X
size[1] = size_y; // size along Y
size[2] = size_z; // size along Z

ImageType::RegionType region;
region.SetSize( size );
region.SetIndex( start );

image->SetRegions( region );
image->Allocate();

image->GetPixelContainer()->SetImportPointer( channel_data, sizeof(unsigned char), false );

return image;
}
*/
