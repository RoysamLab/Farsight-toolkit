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

#include "fregl_result_record.h"

#include <sstream>
#include <iostream>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>

#define MY_ENCODING "ISO-8859-1"

fregl_result_record::
fregl_result_record( std::string const & image_id)
  : reg_image_(image_id)
{}

#if defined(LIBXML_TREE_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

void
fregl_result_record::
write_xml_node(xmlNodePtr root_node)
{
  xmlNodePtr parent_node, node;
  std::string str;
  parent_node = xmlNewChild(root_node, NULL, BAD_CAST "result_set", NULL);
  node = xmlNewChild(parent_node, NULL, BAD_CAST "registration_image", 
                     BAD_CAST reg_image_.c_str());
  if ( nuclear_xml_.length() > 0) {
    node = xmlNewChild(parent_node, NULL, BAD_CAST "nuclei", 
                       BAD_CAST nuclear_xml_.c_str());
  }
  if ( vessel_xml_.length() > 0) {
    node = xmlNewChild(parent_node, NULL, BAD_CAST "vessel", 
                       BAD_CAST vessel_xml_.c_str());
  }
  if ( microglia_xml_.length() > 0) {
    node = xmlNewChild(parent_node, NULL, BAD_CAST "microglia", 
                       BAD_CAST microglia_xml_.c_str());
  }
  if ( astrocyte_xml_.length() > 0) {
    node = xmlNewChild(parent_node, NULL, BAD_CAST "astrocyte", 
                       BAD_CAST astrocyte_xml_.c_str());
  }
}

void
fregl_result_record::
read_xml_node(xmlNodePtr parent_node)
{
  char *contents;
  xmlNodePtr cur_node = NULL;
  std::string str;

  
  cur_node =  parent_node->children;
  for ( ; cur_node; cur_node = cur_node->next) {

    if (cur_node->type != XML_ELEMENT_NODE) continue;

    if (!xmlStrcmp(cur_node->name, BAD_CAST "registration_image") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      reg_image_ = contents;
      xmlFree( contents );
      continue;
    }

    if (!xmlStrcmp(cur_node->name, BAD_CAST "nuclei") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      nuclear_xml_ = contents;
      xmlFree( contents );
      continue;
    }

    if (!xmlStrcmp(cur_node->name, BAD_CAST "vessel") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      vessel_xml_ = contents;
      xmlFree( contents );
      continue;
    }

    if (!xmlStrcmp(cur_node->name, BAD_CAST "microglia") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      microglia_xml_ = contents;
      xmlFree( contents );
      continue;
    }

    if (!xmlStrcmp(cur_node->name, BAD_CAST "astrocyte") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      astrocyte_xml_ = contents;
      xmlFree( contents );
      continue;
    }
  }
}

bool
result_record_read_xml(std::string const & filename, 
                       std::vector<fregl_result_record::Pointer>& result_records)
{
  xmlDocPtr doc; /* the resulting document tree */
  xmlNodePtr root_element, cur_node;

  LIBXML_TEST_VERSION;
    
  //Parse the resource 
  doc = xmlReadFile(filename.c_str(), NULL, 0);
  if (doc == NULL) {
    std::cout<<filename<<" does not yet exist."<<std::endl;
    return false;
  }

  /*Get the root element node */
  root_element = xmlDocGetRootElement(doc);
  
  //Read the parameters
  cur_node =  root_element->children;
  for (; cur_node; cur_node = cur_node->next) {
    if (cur_node->type !=  XML_ELEMENT_NODE) continue;

    fregl_result_record::Pointer result_rec = new fregl_result_record();
    result_rec->read_xml_node(cur_node);
    result_records.push_back( result_rec );
  }

  /*
   * Cleanup function for the XML library.
   */
  xmlCleanupParser();
  /*
   * this is to debug memory for regression tests
   */
  xmlMemoryDump();
  
  return true;
}

void
result_record_write_xml(std::string const & filename, 
                        std::vector<fregl_result_record::Pointer> const& result_records)
{
  xmlDocPtr doc = NULL;       /* document pointer */
  xmlNodePtr root_node = NULL, node = NULL;/* node pointers */
  xmlDtdPtr dtd = NULL;       /* DTD pointer */
  std::string str;

  LIBXML_TEST_VERSION;

  /* 
   * Creates a new document, a node and set it as a root node
   */
  doc = xmlNewDoc(BAD_CAST "1.0");

  root_node = xmlNewNode(NULL, BAD_CAST "Result_sets_for_registration");
  xmlDocSetRootElement(doc, root_node);
  
  for (unsigned i=0; i< result_records.size(); i++) {
    result_records[i]->write_xml_node(root_node);
  }
  
  /* 
   * Dumping document to a file
   */
  xmlSaveFormatFileEnc(filename.c_str(), doc, MY_ENCODING, 1);

  /*free the document */
  xmlFreeDoc(doc);
  
  /*
   *Free the global variables that may
   *have been allocated by the parser.
   */
  xmlCleanupParser();
}

#else 
void
result_record_read_xml(std::string const & filename)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}

void 
result_record_write_xml(std::string const & filename)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}

void
fregl_result_record::
write_xml_node(xmlNodePtr parent_node)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}

void
fregl_result_record::
read_xml_node(xmlNodePtr parent_node)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}
#endif
