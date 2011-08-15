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

#include "fregl_reg_record.h"

#include <sstream>

fregl_reg_record::
fregl_reg_record()
{
  transform_ = NULL;
  overlap_ = 0;
  obj_value_ = 1;
}

fregl_reg_record::
fregl_reg_record( std::string const & from_image_id, 
                  std::string const & to_image_id,
                  SizeType const & from_image_size, 
                  SizeType const & to_image_size)
{
  from_image_id_ = from_image_id;
  to_image_id_ = to_image_id;
  from_image_size_ = from_image_size;
  to_image_size_ = to_image_size;
  transform_ = NULL;
  overlap_ = 0;
  obj_value_ = 1;
}

void 
fregl_reg_record::
set_from_image_info(std::string const & id,
                    SizeType const & image_size)
{
  from_image_id_ = id;
  from_image_size_ = image_size;
}

void 
fregl_reg_record::
set_to_image_info(std::string const & id,
                  SizeType const & image_size)
{
  to_image_id_ = id;
  to_image_size_ = image_size;
}

std::string 
fregl_reg_record::
ToString(double val)
{
    std::ostringstream strm;
    strm<< val<<std::endl;
    return strm.str();
}

void
fregl_reg_record::
read_xml(std::string const & filename)
{
  TiXmlDocument doc;
    
  //Parse the resource
  if ( !doc.LoadFile( filename.c_str() ) ) {
    vcl_cerr<<"Unable to load XML File"<<vcl_endl;
    return ;
  }

  /*Get the root element node */
  TiXmlElement* root_element = doc.RootElement();
  const char* docname = root_element->Value();
  if ( strcmp( docname, "Pairwise_Registration" ) != 0 &&
       strcmp( docname, "Joint_Registration" ) != 0 ) {
    vcl_cerr<<"Incorrect XML root Element "<<vcl_endl;
    return;
  }
  
  //Read the parameters
  TiXmlElement* node = root_element->FirstChildElement();
  read_xml_node( node );
  
}

void 
fregl_reg_record::
write_xml(std::string const & filename)
{
  TiXmlDocument doc; 
  TiXmlElement* root_node = NULL; /* node pointers */

  /* 
   * Creates a new document, a node and set it as a root node
   */

  root_node = new TiXmlElement( "Pairwise_Registration" );
  doc.LinkEndChild( root_node );
  
  write_xml_node(root_node);
 
  /* 
   * Dumping document to a file
   */
  doc.SaveFile(filename.c_str());
  
}

void 
fregl_reg_record::
write_xml_node(TiXmlElement* root_node) 
{
  std::string str;

  TiXmlElement* parent_node = new TiXmlElement("Transform");
  root_node->LinkEndChild(parent_node);

  // from image_id
  TiXmlElement* node = new TiXmlElement("from_image_ID");
  node->LinkEndChild(new TiXmlText(from_image_id_.c_str()));
  str = ToString( from_image_size_[0] );
  node->SetAttribute("size_x", str.c_str());
  str = ToString( from_image_size_[1] );
  node->SetAttribute("size_y", str.c_str());
  str = ToString( from_image_size_[2] );
  node->SetAttribute("size_z", str.c_str());
  parent_node->LinkEndChild(node);

  //to image_id
  node = new TiXmlElement("to_image_ID");
  node->LinkEndChild(new TiXmlText( to_image_id_.c_str() ));
  str = ToString( to_image_size_[0] );
  node->SetAttribute("size_x", str.c_str());
  str = ToString( to_image_size_[1] );
  node->SetAttribute("size_y", str.c_str());
  str = ToString( to_image_size_[2] );
  node->SetAttribute("size_z", str.c_str());
  parent_node->LinkEndChild(node);

  // overlap
  str = ToString( overlap_ );
  node = new TiXmlElement("overlap_percentage");
  node->LinkEndChild(new TiXmlText( str.c_str() ));
  parent_node->LinkEndChild(node);

  // objective value
  str = ToString( obj_value_ );
  node = new TiXmlElement("obj_value");
  node->LinkEndChild(new TiXmlText( str.c_str() ));
  parent_node->LinkEndChild(node);

  // Transformation
  if (transform_) {
    node = new TiXmlElement("parameters");
    parent_node->LinkEndChild(node);
    
    TransformType::ParametersType params = transform_->GetParameters();
    str = ToString( params[0] );
    node->SetAttribute("a00", str.c_str());
    str = ToString( params[1] );
    node->SetAttribute("a01", str.c_str());
    str = ToString( params[2] );
    node->SetAttribute("a02", str.c_str());
    str = ToString( params[3] );
    node->SetAttribute("a10", str.c_str());
    str = ToString( params[4] );
    node->SetAttribute("a11", str.c_str());
    str = ToString( params[5] );
    node->SetAttribute("a12", str.c_str());
    str = ToString( params[6] );
    node->SetAttribute("a20", str.c_str());
    str = ToString( params[7] );
    node->SetAttribute("a21", str.c_str());
    str = ToString( params[8] );
    node->SetAttribute("a22", str.c_str());
    str = ToString( params[9] );
    node->SetAttribute("tx", str.c_str());
    str = ToString( params[10] );
    node->SetAttribute("ty", str.c_str());
    str = ToString( params[11] );
    node->SetAttribute("tz", str.c_str());

	//std::cout << params[11] << std::endl;
  }
}

void 
fregl_reg_record::
read_xml_node(TiXmlElement* parent_node)
{
  std::string str;

  TiXmlElement* cur_node =  parent_node->FirstChildElement();
  for ( ; cur_node; cur_node = cur_node->NextSiblingElement()) {

    const char * value = cur_node->Value();
    
    // from_image_id
    if (strcmp( value, "from_image_ID") == 0 ) {
      from_image_id_ = cur_node->GetText();
      
      from_image_size_[0] = atoi(cur_node->Attribute("size_x"));
      from_image_size_[1] = atoi(cur_node->Attribute("size_y"));
      from_image_size_[2] = atoi(cur_node->Attribute("size_z"));
      
      continue;
    }
    
    // to_image_id
    if (strcmp( value, "to_image_ID") == 0) {
      to_image_id_ = cur_node->GetText();
      to_image_size_[0] = atoi(cur_node->Attribute("size_x"));
      to_image_size_[1] = atoi(cur_node->Attribute("size_y"));
      to_image_size_[2] = atoi(cur_node->Attribute("size_z"));
      
      continue;
    }

    // overlapping
    if (strcmp( value, "overlap_percentage") == 0 ) {
      std::stringstream( cur_node->GetText() ) >> overlap_;
      continue;
    }

    // obj_value
    if ( strcmp( value, "obj_value") == 0 ) {
      std::stringstream( cur_node->GetText() ) >> obj_value_;
      continue;
    }

    // parameters
    if ( strcmp( value, "parameters") == 0 ) {
      TransformType::ParametersType params(12);
      std::stringstream( cur_node->Attribute("a00") ) >> params[0];
      std::stringstream( cur_node->Attribute("a01") ) >> params[1];
      std::stringstream( cur_node->Attribute("a02") ) >> params[2];
      std::stringstream( cur_node->Attribute("a10") ) >> params[3];
      std::stringstream( cur_node->Attribute("a11") ) >> params[4];
      std::stringstream( cur_node->Attribute("a12") ) >> params[5];
      std::stringstream( cur_node->Attribute("a20") ) >> params[6];
      std::stringstream( cur_node->Attribute("a21") ) >> params[7];
      std::stringstream( cur_node->Attribute("a22") ) >> params[8];
      std::stringstream( cur_node->Attribute("tx") ) >> params[9];
      std::stringstream( cur_node->Attribute("ty") ) >> params[10];
      std::stringstream( cur_node->Attribute("tz") ) >> params[11];
      
      transform_ = TransformType::New();
      transform_->SetParameters(params);
      continue;
    }
    
  }
}

