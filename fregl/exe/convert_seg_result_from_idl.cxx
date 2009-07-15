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

//: Executable program to convert the IDL result to the new format
//
//  The usage:
//  
//    convert_seg_result_from_idl _seg_final.dat size_x size_y size_z
//
//  Options:
//    -image    The image_name field for the xml file
//    -xml      Name for the xml file.
//    -class    Name of the _class.[dat/txt] file. Each line contains either "ID  class" or "class" where the line number indicate the ID.
//

#include <sstream>

#include <maciej_seg/nuclear_IDL_util.h>
#include <maciej_seg/xml_util.h>
#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/farsight_maciejseg.h>

#include "itkImageFileWriter.h"

#include <vul/vul_arg.h>
#include <vul/vul_file.h>

typedef itk::Image< unsigned char, 3 > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType >  WriterType;

int
main( int argc, char* argv[] )
{
  /*
  if (argc<5) {
    std::cerr << "Usage: " << argv[0] 
              << " xml_path xml_file image_path image_name"<<std::endl;
    return EXIT_FAILURE;
  }
  
  const char* xml_path = argv[1];
  const char* xml_file = argv[2];
  const char* image_path = argv[3];
  const char* image_name = argv[4];
  */
  
  vul_arg<vcl_string> arg_final_seg (0,"Name of the _seg_final.dat file containing the label information");
  vul_arg<int> arg_sizex("-x","x-dimension");
  vul_arg<int> arg_sizey("-y","y-dimension");
  vul_arg<int> arg_sizez("-z","z-dimension");
  vul_arg<vcl_string> arg_image_name ("-image","The image_name field for the xml file");
  vul_arg<vcl_string> arg_xml_name ("-xml","Name for the xml file. If not specified, it will be image_name.xml");
  vul_arg<vcl_string> arg_class("-class","Name of the _class.[dat/txt] file");
  vul_arg<bool> arg_boundary("-boundary","Showing the cell boundaries", false);
  vul_arg<vcl_string> arg_image_size("-size","The name of an image which provide the size information of the original image");
  
  vul_arg_parse( argc, argv );

  // Get the image size
  if (!arg_image_size.set() && 
      !(arg_sizez.set()&& arg_sizey.set() && arg_sizex.set())) {
    std::cout<<"Either x-y-z dimensions or an image should be provided"<<std::endl;
    return 0;
  }
  int size_z, size_y, size_x;
  if (!arg_image_size.set()) {
    size_x = arg_sizex();
    size_y = arg_sizey();
    size_z = arg_sizez();
  }
  else {
    ImageType::Pointer image;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( arg_image_size() );
    try {
      reader->Update();
    }
    catch(itk::ExceptionObject& e) {
      vcl_cout << e << vcl_endl;
    }
    image =  reader->GetOutput();
    size_x = image->GetRequestedRegion().GetSize()[0];
    size_y = image->GetRequestedRegion().GetSize()[1];
    size_z = image->GetRequestedRegion().GetSize()[2];
  }

  vcl_string image_name, xml_name;
  if (!arg_image_name.set()) { //get the original lsm filename
    std::string image_no_dir = vul_file::strip_directory(arg_final_seg());
    std::string image_no_ext = vul_file::strip_extension(image_no_dir);
    std::string::size_type size = image_no_ext.size();
    image_name = image_no_ext.substr(0,size-10);//_seg_final is length 10
  }
  else image_name = arg_image_name();

  if (!arg_xml_name.set()) { //
    xml_name = vul_file::strip_extension(image_name)+vcl_string(".xml");
  }
  else xml_name = arg_xml_name();

  // Read in the result
  maciejSegmentation marciejSeg;
  if (arg_class.set())
    nuclear_IDL_read_class( arg_class(), arg_final_seg(), size_x, size_y, size_z, marciejSeg);
  else nuclear_IDL_read_seg( arg_final_seg(), size_x, size_y, size_z, marciejSeg);
  
   std::string path = "./";
  // Write out again in the format for maciej segmentation
  xml_util_write( path, xml_name, path, image_name, marciejSeg); 

  // To see the output file
  if (arg_boundary()) {
    farsight_maciejseg ms_wrapper;
    ms_wrapper.read_xml(path, xml_name);
    WriterType::Pointer writer = WriterType::New();
    std::string b_name =  vul_file::strip_extension(image_name)+vcl_string("_boundary_mask.tif");
    writer->SetFileName(b_name);
    writer->SetInput( ms_wrapper.display() );
    writer->Update();
  }

  return 0;
}
