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

//: Executable program to blend 2 3D images and generate 1 2D-projected output
//
#include <vul/vul_arg.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>

#include "itkTIFFImageIO.h"
#include "itkImageFileReader.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkImageFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

/*
ImageType::Pointer
read_image( std::string const & file_name, int channel )
{
std::cout<<"Reading the image "<<file_name<<std::endl;

ImageType::Pointer image;

// Get pixel information
itk::TIFFImageIO::Pointer io = itk::TIFFImageIO::New();
io->SetFileName(file_name);
io->ReadImageInformation();
int pixel_type = (int)io->GetPixelType();
std::cout<<"Pixel Type = "<<pixel_type<<std::endl; //1 - grayscale, 2-RGB, 3-RGBA, etc.,

if (pixel_type == 3) { //RGBA pixel type
typedef fsc_channel_accessor<itk::RGBAPixel<unsigned char>,3 > ChannelAccType;
ChannelAccType channel_accessor(file_name);
image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
}
else if (pixel_type == 2) { //RGA pixel type
typedef fsc_channel_accessor<itk::RGBPixel<unsigned char>,3 > ChannelAccType;
ChannelAccType channel_accessor(file_name);
image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
}
else {// Gray image
typedef itk::ImageFileReader< ImageType > ReaderType;
ReaderType::Pointer reader = ReaderType::New();
reader->SetFileName( file_name );
try {
reader->Update();
}
catch(itk::ExceptionObject& e) {
vcl_cout << e << vcl_endl;
}
image =  reader->GetOutput();
}
fregl_util_reduce_noise( image );

return image ;
}
*/

int
main(  int argc, char* argv[] )
{
	typedef unsigned short																InputPixelType;
	typedef itk::Image< InputPixelType, 3 >												ImageType;
	typedef itk::Image< InputPixelType, 2 >												ImageType2D;	
	typedef itk::Image< float, 3 >														FloatImageType;
	typedef itk::Image< float, 2 >												        FloatImageType2D;
	typedef itk::ImageRegionConstIterator< ImageType >									RegionConstIterator;
	typedef itk::ImageRegionIterator< ImageType2D >										RegionIterator2D;
	typedef fregl_space_transformer< InputPixelType >::TransformType					TransformType;
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType2D, FloatImageType2D >	LoGFilterType2D;
	//typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType, FloatImageType > LoGFilterType;
	typedef itk::SubtractImageFilter< ImageType2D, FloatImageType2D, FloatImageType2D>	SubtractFilterType;
	typedef itk::RescaleIntensityImageFilter< FloatImageType2D, ImageType2D>			RescaleFilterType;
	typedef itk::CastImageFilter< FloatImageType, ImageType >							CastFilterType;
	typedef itk::CastImageFilter< ImageType, FloatImageType >							CastFilterType2;

	vul_arg< vcl_string > arg_xml_file  ( 0, "xml file for transforms" );
	vul_arg< vcl_string > arg_file_from ( 0, "From_image name" );
	vul_arg< vcl_string > arg_file_to   ( 0, "To_image name" );
	vul_arg< int >        arg_channel   ("-channel", "The color channel (0-red, 1-green, 2-blue).",0); 
	vul_arg< vcl_string > arg_img_path  ("-path","The path of the image files.",".");
	vul_arg< vcl_string > arg_old_str   ("-old_str","The old substr in the image names to be replaced");
	vul_arg< vcl_string > arg_new_str   ("-new_str","The new substr in the image names");
	vul_arg< vcl_string > arg_outfile   ("-output","The name of the otuput image", "fused_image.tif");

	vul_arg< bool > arg_in_anchor       ("-in_anchor","The output image is the size of the anchor image", false);
	vul_arg< float > arg_alpha          ("-alpha","The value to power of the photobleaching",4);    

	vul_arg_parse( argc, argv );

	// Cosntruct the graph of joint registration
	fregl_joint_register< InputPixelType >::Pointer joint_register = new fregl_joint_register< InputPixelType >( arg_xml_file() );
	if (arg_old_str.set() && arg_new_str.set()) {
		std::cout<<"Replace the name substr"<<std::endl;
		joint_register->replace_image_name_substr(arg_old_str(), arg_new_str());
	}

	std::string from_image_name = arg_img_path()+std::string("/")+arg_file_from();
	std::string to_image_name = arg_img_path()+std::string("/")+arg_file_to();
	ImageType::Pointer from_image, to_image;
	from_image = fregl_util< InputPixelType >::fregl_util_read_image( from_image_name, arg_channel.set(), arg_channel() );
	to_image = fregl_util< InputPixelType >::fregl_util_read_image( to_image_name, arg_channel.set(), arg_channel() );

	// Set the space transformer
	//
	fregl_space_transformer< InputPixelType > space_transformer(joint_register);

	bool overlap_only = true;
	std::cout<<"Transform "<<arg_file_from()<<" to "<<arg_file_to()<<std::endl;
	space_transformer.set_anchor( arg_file_to(), arg_in_anchor(), overlap_only );

	// Locate the image indices of from_image and to_image
	//
	std::vector<std::string> image_names = space_transformer.image_names();
	int to_index = -1;
	int from_index = -1;
	for (unsigned i = 0; i<image_names.size(); i++){
		if (arg_file_from() == image_names[i]) 
			from_index = i;
		if (arg_file_to() == image_names[i]) 
			to_index = i;
	}
	if (to_index<0) {
		std::cerr<<"To_image is not found!!!"<<std::endl; 
		return 0;
	}
	if (from_index<0) {
		std::cerr<<"From_image is not found!!!"<<std::endl; 
		return 0;
	}

	// Transform the images in 3D
	//
	//ImageType::Pointer xformed_from, xformed_to;
	//xformed_from = space_transformer.transform_image(from_image, from_index);
	//xformed_to = space_transformer.transform_image(to_image, to_index);

	// compute the 2D weight map for photobleaching
	//
	space_transformer.set_individual_weight_map(from_index, from_image, arg_alpha());
	space_transformer.set_individual_weight_map(to_index, to_image, arg_alpha());
	space_transformer.normalize_individual_weight_maps();

	ImageType2D::Pointer from_image_2D = fregl_util< InputPixelType >::fregl_util_max_projection(from_image,0);

	typedef itk::ImageFileWriter< ImageType2D >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "from_proj.tif" );
	writer->SetInput( from_image_2D );
	writer->Update();

	ImageType2D::Pointer to_image_2D = fregl_util< InputPixelType >::fregl_util_max_projection(to_image,0);
	writer->SetFileName( "to_proj.tif" );
	writer->SetInput( to_image_2D );
	writer->Update();


	/*
	typedef itk::ImageFileWriter< ColorImageType >  WriterType3D;
	WriterType3D::Pointer writer3D = WriterType3D::New();
	writer3D->SetFileName( arg_outfile() );
	writer3D->SetInput( out_image );
	writer3D->Update();
	*/

	return 0;
}

