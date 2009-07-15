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

#ifndef fregl_reg_record_h_
#define fregl_reg_record_h_

//:
// \file
// \brief The class which stores 3D registration result
// \author Charlene Tsai
// \date 21 Nov 2007
//

#include <string>

#include "itkAffineTransform.h"
#include "itkSize.h"

//#include "itkLightObject.h"

#include <tinyxml/tinyxml.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>

class fregl_reg_record: public vbl_ref_count
{
public:
  /*
  typedef  fregl_reg_record                  Self;
  typedef  itk::SmartPointer< Self >         Pointer;
  typedef  itk::SmartPointer< const Self >   ConstPointer;

  itkNewMacro(Self) 
  */
  typedef vbl_smart_ptr< fregl_reg_record >  Pointer;
  typedef itk::AffineTransform< double, 3>   TransformType;
  typedef itk::Size<3>                       SizeType;

  //: Default Constructor
  fregl_reg_record();
  
  //: Constructor
  fregl_reg_record( std::string const & from_image_id, 
                    std::string const & to_image_id, 
                    SizeType const & from_image_size, 
                    SizeType const & to_image_size);

  //: Destructor
  ~fregl_reg_record(){};

  //: Set the from_image ID
  void set_from_image_info(std::string const & id,
                           SizeType const & image_size);

  //: Set the to_image ID
  void set_to_image_info(std::string const & id,
                         SizeType const & image_size);

  //: set the transform
  void set_transform(TransformType::Pointer transform) {transform_ = transform;}

  //: set the overlap %
  void set_overlap(double percentage) {overlap_ = percentage;}

  //: set the objective function value
  void set_obj_value(double obj) {obj_value_ = obj;}

  //: get from_image id
  std::string const & from_image() {return from_image_id_;}

  //: get from_image size
  SizeType const& from_image_size() {return from_image_size_;}
  
  //: get to_image_id
  std::string const & to_image() {return to_image_id_;}

  //: get to_image size
  SizeType const& to_image_size() {return to_image_size_;}

  //: get overlap
  double overlap() {return overlap_;}

  //: get objective value
  double obj() {return obj_value_;}

  //: get transform
  TransformType::Pointer transform() {return transform_;}

  // IO with xml file
  void read_xml(std::string const & filename);
  void write_xml(std::string const & filename);
  void write_xml_node(TiXmlElement* parent_node);
  void read_xml_node(TiXmlElement* parent_node);
  
private:
  std::string ToString(double val);

private:
  std::string from_image_id_;
  std::string to_image_id_;
  SizeType from_image_size_;
  SizeType to_image_size_;
  TransformType::Pointer transform_;
  double overlap_;
  double obj_value_;
};

#endif
