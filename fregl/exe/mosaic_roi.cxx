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

//: Executable program to mosaic a set of images with given transformations
//
//  The input is an xml file containing either one pairwise
//  registration result or a set of transformations from joint
//  registration. The images can be gray, rgb color or rgba color. The
//  input images are assumed TIF images. The output montage is a gray
//  image. The usage is
//
//  mosaic_roi xml_file anchor_image -path path_to_images
//
//  where
//    xml_file      Name of the xml_file containing the xforms
//    anchor_image  Name of the anchor image
//
//  Optional arguments"
//    -path         The path of the image files.
//    -old_str      The old substr in the image names to be replaced
//    -new_str      The replacement of the old substr
//    -output       The output image name.


#include <iostream>
#include <string>
using std::cerr;
using std::endl;

#include <vul/vul_arg.h>
#include <vul/vul_file.h>

#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>
#include <fregl/fregl_image_manager.h>
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"

#include "fregl_roi.h"



int
main(int argc, char* argv[]) {
	typedef unsigned short InputPixelType;
	typedef itk::Image< InputPixelType, 3 > ImageType;
	typedef itk::Image< InputPixelType, 2 > ImageType2D;
	typedef ImageType::PointType PointType; //physical space
	typedef ImageType::IndexType IndexType; //physical space
	typedef ImageType::SizeType SizeType;

	vul_arg< vcl_string > arg_xml_file(0, "A xml file containing transformations");
	vul_arg< vcl_string > arg_anchor(0, "Anchor image name");
	vul_arg< double > arg_roi_originx(0, "The x of the roi origin");
	vul_arg< double > arg_roi_originy(0, "The y of the roi origin");
	vul_arg< double > arg_roi_originz(0, "The z of the roi origin");
	vul_arg< unsigned int > arg_roi_sizex(0, "The x-size of the roi origin");
	vul_arg< unsigned int > arg_roi_sizey(0, "The y-size of the roi origin");
	vul_arg< unsigned int > arg_roi_sizez(0, "The z-size of the roi origin");
	vul_arg< vcl_string > arg_outfile(0, "The name of the output file displaying the ROI");
	vul_arg< int > arg_channel("-channel", "The color channel (0-red, 1-green, 2-blue), or the image channel if the original image is a lsm image.", 0);
	vul_arg< vcl_string > arg_img_path("-path", "The path of the image files.", ".");
	vul_arg< vcl_string > arg_old_str("-old_str", "The old substr in the image names to be replaced");
	vul_arg< vcl_string > arg_new_str("-new_str", "The new substr in the image names");
	vul_arg< bool > arg_3d("-3d", "Generate a 3D image as well", false);
	vul_arg< bool > arg_in_anchor("-in_anchor", "The final space is set to the anchor image", false);
	vul_arg< bool > arg_overlap("-overlap_only", "Only consider images that overlap the anchor image", false);
	vul_arg< bool > arg_nn("-nn", "Use Nearest-Neighbor interpolation", false);
	vul_arg< int > arg_blending("-blending", "0: max (default), 1: even weighted, 2: photopleaching weighted (the fanciest).", 0);
	vul_arg< bool > arg_denoise("-denoise", "Making an attempt to remove noise of high frequencies", false);
	vul_arg< bool > arg_use_roi("-use_roi", "Work in roi", false);

	vul_arg_parse(argc, argv);

	//Set up the roi origin and size from the args
	PointType roi_origin;
	SizeType roi_size;
	std::string post_string;
	roi_origin[0] = arg_roi_originx();
	roi_origin[1] = arg_roi_originy();
	roi_origin[2] = arg_roi_originz();
	roi_size[0] = arg_roi_sizex();
	roi_size[1] = arg_roi_sizey();
	roi_size[2] = arg_roi_sizez();	

	
	ImageType::Pointer final_image;
	fregl_roi< InputPixelType > *roi_filter = new fregl_roi< InputPixelType >(arg_xml_file(), arg_img_path(), arg_anchor(), arg_nn());
	roi_filter->setROI(roi_origin, roi_size);
	final_image = roi_filter->getROI();


	// Costruct the montage
	/*fregl_image_manager::Pointer region_montage = new fregl_image_manager(arg_xml_file(), arg_img_path(), arg_anchor(),
			arg_nn());
	
	region_montage->set_regionofinterest(roi_origin, roi_size);*/
	//Taking the maximum
	std::string image_name = arg_img_path() + std::string("/") + arg_anchor();
	//region_montage->Update();
	//final_image = region_montage->GetOutput();
	
	std::string output_name = arg_outfile();

	/*
	//dump the 3d image as series of 2d images
	std::string command = std::string("mkdir ")+output_name;
	  if(vcl_system(command.c_str()) != 0)
	  {
	  cerr << "mkdir returned nonzero" << endl;
	  }

	  typedef itk::NumericSeriesFileNames NameGeneratorType;
	  typedef itk::ImageSeriesWriter< ImageType, ImageType2D > SeriesWriterType;

	  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	  seriesWriter->SetInput( final_image );

	  int last=final_image->GetLargestPossibleRegion().GetSize()[2] - 1;
	  nameGenerator->SetStartIndex( 0 );
	  nameGenerator->SetEndIndex( last );
	  nameGenerator->SetIncrementIndex( 1 );

#if defined(VCL_WIN32) && !defined(__CYGWIN__)
	std::string name_pattern = output_name+std::string("\\slice")+std::string("_%03d.png");
#else
	std::string name_pattern = output_name+std::string("/slice")+std::string("_%03d.png");
#endif
	nameGenerator->SetSeriesFormat( name_pattern );
	seriesWriter->SetFileNames( nameGenerator->GetFileNames() );
	useriesWriter->Update();*/

	if (arg_3d()) {
		std::string name_3d = output_name;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(name_3d);
		writer->SetInput(final_image);
		writer->Update();
	}

	// doing the 2d maximum projection and dump it out
	ImageType2D::Pointer image_2d = fregl_util< InputPixelType >::fregl_util_max_projection(final_image);
	typedef itk::ImageFileWriter< ImageType2D > WriterType2D;
	WriterType2D::Pointer writer2D = WriterType2D::New();
	std::string name_2d = vul_file::strip_extension(output_name) + std::string("_2d_proj.png");
	writer2D->SetFileName(name_2d);
	writer2D->SetInput(image_2d);
	writer2D->Update();

	/*// dump the output to xml
	std::string xml_name = output_name + std::string(".xml");
	fregl_space_transformer::Pointer space_transformer = region_montage->get_space_transformer();
	space_transformer->write_xml(xml_name, output_name, name_2d, arg_overlap(), arg_in_anchor(), arg_channel(), arg_blending(), arg_nn(), arg_denoise());
	*/

	return 0;
}
