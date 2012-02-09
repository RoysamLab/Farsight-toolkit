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

//: Executable to max_project top/bottom n slices to generate 2D images. 
//
// The interface images are needed for alignment of inter-layer
// regisration for registration of the entire brain. If [ext] is
// provided, the images are read from a directory instead. All files
// in the directory will be read.
// 
// For example:
// extract_stack_interfaces directory_name n png  

#include "itkImage.h"
#include <fregl/fregl_util.h>
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include <vul/vul_file_iterator.h>
#include <vul/vul_arg.h>

std::vector<std::string>
gen_names(const std::string dir, const std::string ext)
{
	std::vector<std::string> output;
	std::string pattern = dir + std::string("/*.")+ext;
	for (vul_file_iterator fn=pattern; fn; ++fn) {
		std::cout<<fn()<<std::endl;
		output.push_back(fn());
	}
	return output;
}

int
main(  int argc, char* argv[] )
{
	typedef unsigned short InputPixelType;
	typedef itk::Image< InputPixelType, 3 > ImageType;
	typedef itk::Image< InputPixelType, 2 > ImageType2D;	
	typedef itk::ImageFileWriter< ImageType2D >  WriterType;

	vul_arg< vcl_string > image_arg   ( 0, "The image" );
	vul_arg< int >        slices_arg  ( 0, "Number of slices" );
	vul_arg< int >        channel_arg ("-channel", "The color channel (0-red, 1-green, 2-blue) if not grey-scale.",0);
	/*
	if (argc < 3) {
	std::cerr<<"Syntex: "<<argv[0]<<" tif_image n (slices) [ext]"<<std::endl;
	return EXIT_FAILURE; 
	}
	*/
	vul_arg_parse( argc, argv );

	// Get the base name of the tif image
	std::string image_name = image_arg();
	std::string base_name;
	const std::string slash = "\\/";
	const std::string dot = ".";
	std::string::size_type po = image_name.find_last_of(slash);
	if (po != std::string::npos) image_name = image_name.substr(po+1,image_name.size()-po-1);
	po = image_name.find_last_of(dot);
	base_name = image_name.substr (0, po);
	std::string top_name = base_name+std::string("_top_")+std::string(argv[2])+std::string(".tif");
	std::string bottom_name = base_name+std::string("_bottom_")+std::string(argv[2])+std::string(".tif");

	int num_slices = slices_arg();
	//std::stringstream( argv[2] ) >> num_slices;
	ImageType::Pointer image = fregl_util< InputPixelType >::fregl_util_read_image( image_arg(), channel_arg(), true ); 

	/*
	ImageType::Pointer image; 
	if (argc == 3)
	image = fregl_util_read_image( argv[1], 0, true );
	else { //read in the image from a directory
	ReaderType::Pointer reader = ReaderType::New();
	std::vector<std::string> filenames = gen_names(argv[1], argv[3]);
	reader->SetFileNames(filenames);
	reader->Update();
	image = reader->GetOutput();
	}
	*/

	// Generate the top image
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();
	size[2] = num_slices;
	region.SetSize( size );
	image->SetRequestedRegion(region);
	ImageType2D::Pointer top_image = fregl_util< InputPixelType >::fregl_util_max_projection( image );
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( top_name );
	writer->SetInput( top_image );
	writer->Update();

	// Generate the bottom image
	region = image->GetLargestPossibleRegion();
	size = region.GetSize();
	ImageType::IndexType start = region.GetIndex();
	start[2] = size[2]-num_slices;
	size[2] = num_slices;
	region.SetSize( size );
	region.SetIndex( start );
	image->SetRequestedRegion(region);
	ImageType2D::Pointer bottom_image = fregl_util< InputPixelType >::fregl_util_max_projection( image );
	writer->SetFileName( bottom_name );
	writer->SetInput( bottom_image );
	writer->Update();

	return 0;
}
