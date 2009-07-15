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

//: Executable program to store the segmentation xml file names with the image ID for registration
//
//  The useage is:
//  
//       update_result_sets result_sets_xml_file reg_image_id
//
//  optional arguments: 
//  -nuclear   
//  -vessel

#include <fregl/fregl_result_record.h>
#include <vul/vul_arg.h>

int
main(  int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_xml_file  ( 0, "xml file" );
  vul_arg< vcl_string > arg_reg_image_id (0,"Registration image ID");
  vul_arg< vcl_string > arg_nuclear ("-nuclear","Nuclear segmentation xml");
  vul_arg< vcl_string > arg_vessel ("-vessel","Vessel .npts file");
  vul_arg< vcl_string > arg_microglia ("-microglia","microglia trace point file (*TracePoints.txt)");
  vul_arg< vcl_string > arg_astrocyte ("-astrocyte","astrocyte trace point file (*TracePoints.txt)");
  vul_arg_parse( argc, argv );

  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(arg_xml_file(), result_records);

  // Locate the result record with the same reg_image_id. If not
  // found, create a new record.
  fregl_result_record::Pointer record;
  for (unsigned i = 0; i<result_records.size(); i++) {
    if (result_records[i]->registration_image() == arg_reg_image_id()) {
      record = result_records[i];
    }
  }
  
  if (!record) {
    record = new fregl_result_record(arg_reg_image_id());
    result_records.push_back( record );
  }

  if (arg_nuclear.set()) record->set_nuclear_xml( arg_nuclear() );
  if (arg_vessel.set()) record->set_vessel_xml( arg_vessel() );
  if (arg_microglia.set()) record->set_microglia_xml( arg_microglia() );
  if (arg_astrocyte.set()) record->set_astrocyte_xml( arg_astrocyte() );

  result_record_write_xml(arg_xml_file(), result_records);

  return 0;
}
