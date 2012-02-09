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

//: Executable program to mosaic an image pair with given transformations
//
//  The input is an xml file containing the registration result. The
//  images can be gray, rgb color or rgba color. The input images are
//  assumed TIF images. The output montage is a color image with
//  images montaged to different color channels. The usage is
//  
//  mosaic_images xml_file from_image to_image
//
//  where
//    xml_file      Name of the xml_file containing the xforms
//    from_image    Name of the from_image
//    to_image      Name of the to_image
//
//  Optional arguments"
//    -path         The path of the image files.
//    -old_str      The old substr in the image names to be replaced
//    -new_str      The replacement of the old substr
//    -output       The output image name.
//    -channel      Needed if the input image is rgb or rgba

#include <iostream>
using std::cerr;
using std::endl;

#include <vul/vul_arg.h>
#include <vul/vul_file.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>

#include "itkTIFFImageIO.h"
#include "itkImageFileReader.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"



int
main(  int argc, char* argv[] )
{
	typedef unsigned short								InputPixelType;
	typedef itk::Image< InputPixelType, 3>				ImageType;
	typedef itk::RGBPixel< unsigned char >				OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 >			ColorImageType;
	typedef itk::Image< OutputPixelType, 2 >			ColorImageType2D;
	typedef itk::ImageRegionConstIterator< ImageType >	RegionConstIterator;
	typedef itk::ImageRegionIterator< ColorImageType >	RegionIterator;

	vul_arg< vcl_string > arg_xml_file  ( 0, "xml file for transforms" );
	vul_arg< vcl_string > arg_file_from ( 0, "From_image name" );
	vul_arg< vcl_string > arg_file_to   ( 0, "To_image name" );
	vul_arg< int >        arg_channel   ("-channel", "The color channel (0-red, 1-green, 2-blue), or the image channel if the original image is lsm image.",0); 
	vul_arg< vcl_string > arg_img_path  ("-path","The path of the image files.",".");
	vul_arg< vcl_string > arg_old_str   ("-old_str","The old substr in the image names to be replaced");
	vul_arg< vcl_string > arg_new_str   ("-new_str","The new substr in the image names");
	vul_arg< vcl_string > arg_outfile   ("-output","The name of the otuput image");

	vul_arg< bool > arg_in_anchor       ("-in_anchor","The output image is the size of the anchor image", false);

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

	// bool overlap_only = true;
	bool overlap_only = false;
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

	// Transform the images
	//
	ImageType::Pointer xformed_from, xformed_to;
	xformed_from = space_transformer.transform_image(from_image, from_index);
	from_image = 0;
	xformed_to = space_transformer.transform_image(to_image, to_index);
	to_image = 0;

	/*
	typedef itk::ImageFileWriter< ImageType >  WriterType3D;

	WriterType3D::Pointer writer3D = WriterType3D::New();
	writer3D->SetFileName( arg_outfile() );
	writer3D->SetInput( xformed_from );
	writer3D->Update();
	*/

	// Fuse the images
	ColorImageType::Pointer out_image = ColorImageType::New();
	out_image->SetRegions( xformed_to->GetRequestedRegion() );
	out_image->Allocate();
	out_image->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));

	RegionConstIterator inputIt1( xformed_from, xformed_from->GetRequestedRegion() );
	RegionConstIterator inputIt2( xformed_to, xformed_to->GetRequestedRegion() );
	RegionIterator outputIt( out_image, out_image->GetRequestedRegion() );

	for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); 
		!inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++outputIt) {
			OutputPixelType pix = outputIt.Get();
			pix.SetRed( inputIt1.Get() );
			pix.SetGreen( inputIt2.Get() );
			outputIt.Set( pix );
	}

	// dump the 3d images as 2d slices in a directory
	std::string name_prefix = std::string("pairwise_montage_")+vul_file::strip_extension( arg_file_to() );
	if ( arg_outfile.set() ) {
		name_prefix = arg_outfile();
	}

	/*
	std::string command = std::string("mkdir ")+name_prefix;
	if(vcl_system(command.c_str()) != 0)
	{
	cerr << "mkdir returned nonzero" << endl;
	}

	typedef itk::NumericSeriesFileNames NameGeneratorType;
	typedef itk::ImageSeriesWriter< ColorImageType, ColorImageType2D > SeriesWriterType;

	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	seriesWriter->SetInput( out_image );

	int last=out_image->GetLargestPossibleRegion().GetSize()[2] - 1;
	nameGenerator->SetStartIndex( 0 );
	nameGenerator->SetEndIndex( last );
	nameGenerator->SetIncrementIndex( 1 );

	#if defined(VCL_WIN32) && !defined(__CYGWIN__)
	std::string name_pattern = name_prefix+std::string("\\slice")+std::string("_%03d.png");
	#else
	std::string name_pattern = name_prefix+std::string("/slice")+std::string("_%03d.png");
	#endif
	nameGenerator->SetSeriesFormat( name_pattern );
	seriesWriter->SetFileNames( nameGenerator->GetFileNames() );
	seriesWriter->Update();*/

	// doing the 2d maximum projection and dump it out
	ColorImageType2D::Pointer image_2d = fregl_util< InputPixelType >::fregl_util_max_projection_color(out_image);
	typedef itk::ImageFileWriter< ColorImageType2D >  WriterType2D;
	WriterType2D::Pointer writer2D = WriterType2D::New();
	std::string name_2d = name_prefix + std::string("_2d_proj.png");
	writer2D->SetFileName( name_2d );
	writer2D->SetInput( image_2d );
	writer2D->Update();

	// dump out 3D the montage
	typedef itk::ImageFileWriter< ColorImageType >  WriterType3D;
	WriterType3D::Pointer writer3D = WriterType3D::New();
	writer3D->SetFileName( name_prefix+".tif" );
	writer3D->SetInput( out_image );
	writer3D->Update();


	// dump the output to xml
	//std::string name_no_ext = vul_file::strip_extension( arg_outfile() );
	//std::string xml_name = name_no_ext+std::string(".xml");
	//space_transformer.write_xml( xml_name, arg_outfile(), std::string(""), true, arg_in_anchor(), arg_channel());
	std::string xml_name = name_prefix+std::string(".xml");
	space_transformer.write_xml( xml_name, name_prefix, name_2d, true, arg_in_anchor(), arg_channel());
	return 0;
}
