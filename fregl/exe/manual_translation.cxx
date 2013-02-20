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

//: Executable program to generate an xml file containing only translations given by the user
//
//  The input is a file containing the image names, image sizes, and the translation values
//  
//  manual_translation input_file anchor sx sy sz
//
//  where
//   input_file  A text file containing image names and the given translations
//   anchor      Name of the anchor image
//   sx          x-dimension of the anchor image
//   sy          y-dimension of the anchor image
//   sz          z-dimension of the anchor image
//
//  Optional:
//    -output    The xml file containing the transformations in the format needed for the mosaicking executables
//
//   The input file is in the format of:
//     from_image_1 sx1 sy1 sz1 tx1 ty1 tz1
//     from_image_2 sx2 sy2 sz2 tx2 ty2 tz2
//     ...
//

#include <vul/vul_arg.h>

#include <iostream>
#include <fstream>

#include <fregl/fregl_reg_record.h>
#include <fregl/fregl_util.h>
#include <fregl/fregl_joint_register.h>

#include "itkAffineTransform.h"
#include "itkSize.h"

typedef unsigned short						InputPixelType;
typedef itk::AffineTransform< double, 3>	TransformType;
typedef itk::Size<3>						SizeType;

int 
main( int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_input  ( 0, "Name of the input file containing the image names, image sizes, and the translation values" );
 
  vul_arg< vcl_string > arg_anchor  ( 0, "The anchor image" );
  vul_arg< int > arg_sx  ( 0, "x-dimension of the anchor image" );
  vul_arg< int > arg_sy  ( 0, "y-dimension of the anchor image" );
  vul_arg< int > arg_sz  ( 0, "z-dimension of the anchor image" );
  vul_arg< vcl_string > arg_output ( "-output", "Name of the xml file containing the transformations. If not set, it is the the anchor imgae name appended with _joint.xml" );
  vul_arg_parse( argc, argv );

  // Set the anchor image if given
  std::string anchor_image = arg_anchor();
  SizeType anchor_size;
  anchor_size[0]=arg_sx();
  anchor_size[1]=arg_sy();
  anchor_size[2]=arg_sz();
  
  // Read in the input file
  std::ifstream inList;
  inList.open(arg_input().c_str());
  if ( !inList ){
    std::cerr<<"Couldn't open "<<arg_input()<<std::endl;
    exit( 0 );
  }

  std::string line_str;
  std::string from_image;
  SizeType from_size;
  TransformType::ParametersType parameters(12);
  for (int i = 0; i<12; i++) parameters[i] = 0;
  parameters[0] = 1;
  parameters[4] = 1;
  parameters[8] = 1;
  
  std::vector<fregl_reg_record::Pointer> reg_records;
  int image_count = 1;
  while ( inList ) {
    std::getline( inList, line_str );
    if (line_str.length() == 0) continue;

    std::istringstream line_stream(line_str);
    line_stream>> from_image >> from_size[0] >>from_size[1] >> from_size[2]>> parameters[9] >> parameters[10] >> parameters[11];
    fregl_reg_record::Pointer reg_record = new fregl_reg_record(from_image, anchor_image, from_size, anchor_size);
    TransformType::Pointer transform = TransformType::New();
    transform->SetParameters( parameters );
    reg_record->set_transform(transform);
    reg_record->set_obj_value( 0 );
	double vol_overlap = fregl_util< InputPixelType >::fregl_util_overlap(transform, from_size, anchor_size);
    reg_record->set_overlap(vol_overlap);
    reg_records.push_back(reg_record);
    
    image_count++;
  }
  
  // Run joint-registration without mutual consistency
  fregl_joint_register< InputPixelType >::Pointer joint_register = new fregl_joint_register< InputPixelType >( reg_records, 0, 1 );
  joint_register->infer_graph();
    
  // Output to xml file
  std::string output_file;
  if (!arg_output.set()) output_file = anchor_image+std::string("_joint.xml");
  else output_file = arg_output();
  joint_register->write_xml(output_file, false, false);
}
