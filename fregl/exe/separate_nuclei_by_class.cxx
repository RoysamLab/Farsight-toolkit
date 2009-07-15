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

//: Executable program to separate the mask image by class type.
//
//  The usage:
//   
//      separate_nuclei_by_class xml_file -class class_type -image
//
//  Either or both options should be chosen to generate the output. 

#include <maciej_seg/xml_util.h>
#include <maciej_seg/farsight_maciejseg.h>

#include <string>

#include <vul/vul_file.h>
#include <vul/vul_arg.h>
#include <vnl/vnl_vector_fixed.h>

#include "itkImage.h"
#include "itkImageFileWriter.h"

/*
void separate_dirname_from_filename( const std::string& input,
                                     std::string& dirname,
                                     std::string& filename)
{
  std::string base_name;
  const std::string slash = "\\/";
  std::string::size_type po = input.find_last_of(slash);
  if (po != std::string::npos) {
    filename = input.substr(po+1,filename.size()-po-1);
    dirname = input.substr(0,po+1);
  }
  else {
    filename = input;
    dirname = "";
  }
}
*/

typedef itk::Image< unsigned char, 3 > ImageType;

int 
main( int argc, char* argv[] )
{
  /*
  if (argc < 2) {
    std::cerr<<"Usage: "<<argv[0]<<" xml_file"<<std::endl;
    return 0;
  }
  */
  
  vul_arg<vcl_string> arg_xml_file (0,".xml file containing nuclear segmentation.");
  vul_arg<int> arg_class("-class","Only separate out the specific class type"); 
  vul_arg<bool> arg_image("-image","Generate ID mask images of individual classes", false);
  vul_arg<bool> arg_text("-text","Generate the text files, one for each class, containing the locations of the boundary points only", false);
  
  vul_arg_parse( argc, argv );

  if (!arg_image.set() && !arg_text.set()) {
    std::cerr<<"Either -image or -text should be specified"<<std::endl;
    return 0;
  }

  
  std::string path = vul_file::dirname( arg_xml_file() ) + std::string("/");
  std::string xml_name = vul_file::strip_directory( arg_xml_file() );
  //std::string path, xml_name, name_no_xml;
  //separate_dirname_from_filename(argv[1], path, xml_name);
  std::string name_no_xml = vul_file::strip_extension( xml_name );
  maciejSegmentation::Pointer seg_result = new maciejSegmentation();
  std::string image_path, image_name;
  xml_util_read( path, xml_name, image_path, image_name, *seg_result );

  // count the number of classes
  int class_count = 0;
  std::vector<rich_cell::Pointer> const& cells = seg_result->all_cells();
  for (unsigned i = 0; i<cells.size(); i++) {
    if (class_count < cells[i]->class_type_) class_count = cells[i]->class_type_;
  }
  if (arg_class()) {
    class_count = vnl_math_max( class_count, arg_class());
  }

  // Get the image dimension
  vnl_vector_fixed<int,3> const&  image_size = seg_result->image_size(); 
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  
  ImageType::SizeType size;
  size[0] = image_size[0]; // size along X
  size[1] = image_size[1]; // size along Y
  size[2] = image_size[2]; // size along Z
  
  ImageType::RegionType region;
  region.SetSize( size );

  // For each class generate one mask image and one text file if specified
  //
  for (int ci = 1; ci<= class_count; ci++) {
    if (arg_class.set() && ci != arg_class()) continue;

    char class_num[5] ;
    std::sprintf( class_num, "%d", ci);
  
    // For mask image
    if ( arg_image.set() ) {
      std::cout<<"Constructing the boundary image for class "<<ci<<std::endl;
      std::string file_name = name_no_xml+std::string("_class")+ class_num+std::string(".tif");
      ImageType::Pointer class_label_image = ImageType::New();
      class_label_image->SetRegions( region );
      class_label_image->Allocate();
      class_label_image->FillBuffer(0);

      for (unsigned int i = 0; i<cells.size(); i++) {
        if (cells[i]->class_type_ != ci || !cells[i]->valid_ || cells[i]->dup_) continue;
        
        const std::vector< vnl_vector_fixed< float, 3 > >& points = cells[i]->boundary_points_;
        for (unsigned int pi = 0; pi<points.size(); pi++) {
          vnl_vector_fixed< float, 3 > const & pt = points[pi];
          ImageType::IndexType pos;
          pos[0] = pt[0];
          pos[1] = pt[1];
          pos[2] = pt[2];
          class_label_image->SetPixel(pos, 255);
        }
      }

      std::cout<<"Generate boundary image "<<file_name<<std::endl;
      typedef itk::ImageFileWriter< ImageType >  WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(file_name);
      writer->SetInput( class_label_image );
      writer->Update();
    }

    // For text file
    if (arg_text.set()) {
      std::string file_name = name_no_xml+std::string("_class")+ class_num+std::string(".txt");
      std::ofstream outfile(file_name.c_str());
      for (int i = 0; i<cells.size(); i++) {
        if (cells[i]->class_type_ != ci || !cells[i]->valid_ || cells[i]->dup_) continue;
        
        const std::vector< vnl_vector_fixed< float, 3 > >& points = cells[i]->boundary_points_;
        for (int pi = 0; pi<points.size(); pi++) {
          // z,y,x
          outfile << points[pi][2]<<"\t"<<points[pi][1]<<"\t"<<points[pi][0]<<"\t0\t0\n";
        }
      
      }
      outfile.close();
    }

  }

  return 0;
}
