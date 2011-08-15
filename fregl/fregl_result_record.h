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

#ifndef fregl_result_record_h_
#define fregl_result_record_h_

//:
// \file
// \brief The class which stores record to relate image ID for segmentaion to registration result
// \author Charlene Tsai
// \date 01 Jan 2008
//

#include <string>
#include <vector>

#include <libxml/parser.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>

class fregl_result_record: public vbl_ref_count
{
public:
  /*
  typedef  fregl_reg_record                  Self;
  typedef  itk::SmartPointer< Self >         Pointer;
  typedef  itk::SmartPointer< const Self >   ConstPointer;

  itkNewMacro(Self) 
  */
  typedef vbl_smart_ptr< fregl_result_record >  Pointer;

  //: Default Constructor
  fregl_result_record(){};
  
  //: Constructor
  fregl_result_record( std::string const & image_id);

  //: Destructor
  ~fregl_result_record(){};

  //: Set the registration image
  void set_registration_image(std::string const & id) {reg_image_ = id;}

  //: Set the nuclear file id
  void set_nuclear_xml(std::string const & id) {nuclear_xml_ = id;}

    //: set the vessel file id
  void set_vessel_xml(std::string const & id) {vessel_xml_ = id;}

  //: set the microglia file id
  void set_microglia_xml(std::string const & id) {microglia_xml_ = id;}

  //: set the astrocyte file id
  void set_astrocyte_xml(std::string const & id) {astrocyte_xml_ = id;}

  //: get registration image id
  std::string const & registration_image() {return reg_image_;}

  //: get nuclear file id
  std::string const & nuclear_xml() {return nuclear_xml_;}

  //: get vessel file id
  std::string const & vessel_xml() {return vessel_xml_;}

  //: get microglia file id
  std::string const & microglia_xml() {return microglia_xml_;}

  //: get astrocyte file id
  std::string const & astrocyte_xml() {return astrocyte_xml_;}

  // IO with xml file

  void write_xml_node(xmlNodePtr parent_node);
  void read_xml_node(xmlNodePtr parent_node);
 
private:
  std::string ToString(double val);

private:
  std::string reg_image_;
  std::string nuclear_xml_;
  std::string vessel_xml_;
  std::string microglia_xml_;
  std::string astrocyte_xml_;
};

bool result_record_read_xml(std::string const & filename,
                            std::vector<fregl_result_record::Pointer>& result_records);
void result_record_write_xml(std::string const & filename,
                             std::vector<fregl_result_record::Pointer> const& result_records);

#endif
