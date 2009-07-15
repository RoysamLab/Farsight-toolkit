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

#include <maciej_seg/xml_util.h>
#include <maciej_seg/farsight_maciejseg.h>

#include <string>
#include <vector>

#include "itkImageFileWriter.h"

int main( int argc , char * argv[] )
{
  /*
  // To test if the result stored in xml can be read back to memory
  // properly
  // The input list contains: 1. filepath, 2. xml_filename
  if (argc != 3) {
    std::cerr <<"command: nucleus_seg xml_filepath xml_filename" << std::endl;
    return 0;
  }

  std::string xml_path = argv[1];
  std::string xml_filename = argv[2];
  farsight_maciejseg ms_wrapper;
  ms_wrapper.read_xml(xml_path, xml_filename);
  */

  // The input list contains: 1. image path, 2. image_name,
  // 3. param_xml, 4. output_xml path, 5. output_xml_file
  if (argc != 6) {
    std::cerr <<"command: nucleus_seg image_path image_name param_xml_file output_xml_path output_xml" << std::endl;
    return 0;
  }
  
  farsight_maciejseg ms_wrapper;
  std::string xml_param_filename = argv[3];
  ms_wrapper.set_parameters( xml_param_filename );
  std::string image_path = argv[1];
  std::string image_name = argv[2];
  ms_wrapper.set_image(image_path, image_name);
  ms_wrapper.run();
  std::string xml_path = argv[4];
  std::string xml_filename = argv[5];
  ms_wrapper.write_xml(xml_path, xml_filename);


  // To see the output file
  /*
  typedef itk::Image< unsigned char, 3 > ImageType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "seg_result2.tif" );

  writer->SetInput( ms_wrapper.display() );
  writer->Update();
  */
  /*
  // Write out the default parameters to the xml file
  std::vector<std::string> parameters;
  std::string param = "1"; //channel
  parameters.push_back(param);
  param = "500"; //threshold
  parameters.push_back(param);
  param = "4"; // grid
  parameters.push_back(param);
  param = "1"; //rM
  parameters.push_back(param);
  param = "1"; //rG
  parameters.push_back(param);
  param = "1"; //mode
  parameters.push_back(param);
  param = ""; //glia
  parameters.push_back(param);
  param = ""; //neuron
  parameters.push_back(param);

  std::string xml_filename = "nucleus_params.xml";
  xml_util_write_param( xml_filename, parameters);
  */

  return 0;
}
