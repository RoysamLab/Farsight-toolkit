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

#include "xml_util.h"

//#include <libxml/encoding.h>
//#include <libxml/xmlwriter.h>

#include <iostream>
#include <sstream>
#include <string>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>

#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#define MY_ENCODING "ISO-8859-1"


#if defined(LIBXML_TREE_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

static rich_cell::Pointer read_one_nucleus( xmlNodePtr node );

/***************** xml_util_write ***********************************/

bool xml_util_write( std::string const& file_path,
                     std::string const& xml_filename, 
                     std::string const& image_path,
                     std::string const& image_name, 
                     maciejSegmentation const& marciejSeg )
{
  xmlDocPtr doc = NULL;       /* document pointer */
  xmlNodePtr root_node = NULL, node = NULL, node1 = NULL, node2 = NULL, node3 = NULL;/* node pointers */
  xmlDtdPtr dtd = NULL;       /* DTD pointer */

  LIBXML_TEST_VERSION;

  /* 
   * Creates a new document, a node and set it as a root node
   */
  std::vector<std::string> parameters = marciejSeg.parameters();

  doc = xmlNewDoc(BAD_CAST "1.0");
  root_node = xmlNewNode(NULL, BAD_CAST "nuclei");
  xmlNewProp(root_node, BAD_CAST "program", BAD_CAST "Maciej Segmentation");
  xmlNewProp(root_node, BAD_CAST "version", BAD_CAST "1.0");
  xmlNewProp(root_node, BAD_CAST "image_path", BAD_CAST image_path.c_str());
  xmlNewProp(root_node, BAD_CAST "image", BAD_CAST image_name.c_str());

  //xmlNewProp(root_node, BAD_CAST "type", BAD_CAST parameters[2].c_str());
  //xmlNewProp(root_node, BAD_CAST "sstack", BAD_CAST parameters[3].c_str());
  //xmlNewProp(root_node, BAD_CAST "estack", BAD_CAST parameters[4].c_str());
  //xmlNewProp(root_node, BAD_CAST "swidth", BAD_CAST parameters[5].c_str());

  xmlNewProp(root_node, BAD_CAST "channel", BAD_CAST parameters[0].c_str());
  xmlNewProp(root_node, BAD_CAST "threshold", BAD_CAST parameters[1].c_str());
  xmlNewProp(root_node, BAD_CAST "grid", BAD_CAST parameters[2].c_str());
  xmlNewProp(root_node, BAD_CAST "median_filter_radius", BAD_CAST parameters[3].c_str());
  xmlNewProp(root_node, BAD_CAST "morph_opt_radius", BAD_CAST parameters[4].c_str());
  xmlNewProp(root_node, BAD_CAST "mode", BAD_CAST parameters[5].c_str());
  xmlNewProp(root_node, BAD_CAST "gliadat", BAD_CAST parameters[6].c_str());
  xmlNewProp(root_node, BAD_CAST "neurdat", BAD_CAST parameters[7].c_str());

  std::string label_image_name = image_name + std::string("_label_img.tif");
  xmlNewProp(root_node, BAD_CAST "label_image", BAD_CAST label_image_name.c_str());

  xmlDocSetRootElement(doc, root_node);
  
  /* 
   * Iterate through the cells and attach each child node
   * to root_node node. 
   */

  std::vector<rich_cell::Pointer> const& cells = marciejSeg.all_cells();
  for (unsigned i = 0; i<cells.size(); i++) {
    std::ostringstream sstr;
    std::string str;

    //if (cells[i] == NULL) continue;

    node = xmlNewChild(root_node, NULL, BAD_CAST "nucleus", NULL);

    str = ToString((double)cells[i]->valid_);    
    xmlNewChild(node, NULL, BAD_CAST "validity", BAD_CAST str.c_str());

    str = ToString((double)cells[i]->dup_);    
    xmlNewChild(node, NULL, BAD_CAST "duplicated", BAD_CAST str.c_str());

    str = ToString((double)cells[i]->label_);
    xmlNewChild(node, NULL, BAD_CAST "label", BAD_CAST str.c_str());
 
    str = ToString(cells[i]->volume_);
    xmlNewChild(node, NULL, BAD_CAST "volume", BAD_CAST str.c_str());
    
    str = ToString(cells[i]->ave_intensity_);
    xmlNewChild(node, NULL, BAD_CAST "intensity", BAD_CAST str.c_str());
   
    str = ToString(cells[i]->texture_);
    xmlNewChild(node, NULL, BAD_CAST "texture", BAD_CAST str.c_str());

    str = ToString(cells[i]->eccentricity_);
    xmlNewChild(node, NULL, BAD_CAST "eccentricity", BAD_CAST str.c_str());

    str = ToString(cells[i]->average_radius_);
    xmlNewChild(node, NULL, BAD_CAST "average_radius", BAD_CAST str.c_str());

    str = ToString(cells[i]->neuronal_signal_);
    xmlNewChild(node, NULL, BAD_CAST "neuronal_signal", BAD_CAST str.c_str());

    str = ToString(cells[i]->nearest_nbr_dist_);
    xmlNewChild(node, NULL, BAD_CAST "nearest_nbr_dist", BAD_CAST str.c_str());

    str = ToString(cells[i]->nearest_nbr_label_);
    xmlNewChild(node, NULL, BAD_CAST "nearest_nbr_label", BAD_CAST str.c_str());

    str = ToString(cells[i]->score_);
    xmlNewChild(node, NULL, BAD_CAST "score", BAD_CAST str.c_str());

    str = ToString(cells[i]->class_type_);
    xmlNewChild(node, NULL, BAD_CAST "class", BAD_CAST str.c_str());

    /*
    str = ToString(cells[i]->convexity_);
    xmlNewChild(node, NULL, BAD_CAST "convexity", BAD_CAST str.c_str());

    str = ToString(cells[i]->shape_factor_);
    xmlNewChild(node, NULL, BAD_CAST "shape_factor", BAD_CAST str.c_str());

    str = ToString(cells[i]->bending_energy_);
    xmlNewChild(node, NULL, BAD_CAST "bending_energy", BAD_CAST str.c_str());

    str = ToString(cells[i]->ave_bound_gradient_);
    xmlNewChild(node, NULL, BAD_CAST "boundary_gradient", BAD_CAST str.c_str());

    str = ToString(cells[i]->vol_grad_);
    xmlNewChild(node, NULL, BAD_CAST "volume_gradient", BAD_CAST str.c_str());

    str = ToString(cells[i]->radius_variance_);
    xmlNewChild(node, NULL, BAD_CAST "radius_variance", BAD_CAST str.c_str());

    str = ToString(cells[i]->bound_ints_ratio_);
    xmlNewChild(node, NULL, BAD_CAST "intensity_ratio", BAD_CAST str.c_str());

    str = ToString(cells[i]->percent_nbr_);
    xmlNewChild(node, NULL, BAD_CAST "shared_boundary", BAD_CAST str.c_str());


    str = ToString(cells[i]->eng_intensity_dist_);
    xmlNewChild(node, NULL, BAD_CAST "engergy_intensity_distribution", BAD_CAST str.c_str());

    str = ToString(cells[i]->entropy_intensity_dist_);
    xmlNewChild(node, NULL, BAD_CAST "entropy_intensity_distribution", BAD_CAST str.c_str());

    str = ToString(cells[i]->intensity_variation_);
    xmlNewChild(node, NULL, BAD_CAST "intensity_variation", BAD_CAST str.c_str());

    str = ToString(cells[i]->num_poles_on_surface_);
    xmlNewChild(node, NULL, BAD_CAST "num_of_poles_on_surface", BAD_CAST str.c_str());

    str = ToString(cells[i]->skew_intensity_dist_);
    xmlNewChild(node, NULL, BAD_CAST "skew_intensity_distribution", BAD_CAST str.c_str());

    str = ToString(cells[i]->surface_area_);
    xmlNewChild(node, NULL, BAD_CAST "surface_area", BAD_CAST str.c_str());

    str = ToString(cells[i]->orientation_);
    xmlNewChild(node, NULL, BAD_CAST "orientation", BAD_CAST str.c_str());

    */
    str = ToString(cells[i]->ave_nnbr_dist_);
    xmlNewChild(node, NULL, BAD_CAST "ave_nnbr_distance", BAD_CAST str.c_str());

    // center node
    node1 = xmlNewChild(node, NULL, BAD_CAST "center", NULL);
    str = ToString(cells[i]->center_[0]);
    xmlNewChild(node1, NULL, BAD_CAST "x", BAD_CAST str.c_str());
    str = ToString(cells[i]->center_[1]);
    xmlNewChild(node1, NULL, BAD_CAST "y", BAD_CAST str.c_str());
    str = ToString(cells[i]->center_[2]);
    xmlNewChild(node1, NULL, BAD_CAST "z", BAD_CAST str.c_str());
    
    // bounding box
    rich_cell::RegionType::IndexType start = cells[i]->bounding_box_.GetIndex();
    rich_cell::RegionType::SizeType size = cells[i]->bounding_box_.GetSize();
    node1 = xmlNewChild(node, NULL, BAD_CAST "bounding_box", NULL);
    str = ToString(start[0]);
    xmlNewChild(node1, NULL, BAD_CAST "min_x", BAD_CAST str.c_str());
    str = ToString(start[1]);
    xmlNewChild(node1, NULL, BAD_CAST "min_y", BAD_CAST str.c_str());
    str = ToString(start[2]);
    xmlNewChild(node1, NULL, BAD_CAST "min_z", BAD_CAST str.c_str());
    str = ToString(size[0]);
    xmlNewChild(node1, NULL, BAD_CAST "size_x", BAD_CAST str.c_str());
    str = ToString(size[1]);
    xmlNewChild(node1, NULL, BAD_CAST "size_y", BAD_CAST str.c_str());
    str = ToString(size[2]);
    xmlNewChild(node1, NULL, BAD_CAST "size_z", BAD_CAST str.c_str());

    /*
    // neighbors
    if (cells[i]->neighbors_.size() != 0) {
      node1 = xmlNewChild(node, NULL, BAD_CAST "neighbors", NULL);
      for (unsigned int j=0; j<cells[i]->neighbors_.size(); j++) {
        str = ToString((double)cells[i]->neighbors_[j]);
        xmlNewChild(node1, NULL, BAD_CAST "label", BAD_CAST str.c_str());
      }
    }
    */

    // edit history
    if (cells[i]->edit_history_.size() != 0) {
      node1 = xmlNewChild(node, NULL, BAD_CAST "edit_history", NULL);
      unsigned int num_rec = cells[i]->edit_history_.size();
      for (unsigned int j = 0; j<num_rec; j++) {
        edit_record const& rec =  cells[i]->edit_history_[j];
        node2 = xmlNewChild(node1, NULL, BAD_CAST "record", NULL);
        xmlNewProp(node2, BAD_CAST "date", BAD_CAST rec.date_.c_str() );
        xmlNewProp(node2, BAD_CAST "name", BAD_CAST rec.name_.c_str() );

        switch (rec.status_) {
        case edit_record::ADDED:
          xmlNewProp(node2, BAD_CAST "status", BAD_CAST "ADDED");
          break;
        case edit_record::DELETED:
          xmlNewProp(node2, BAD_CAST "status", BAD_CAST "DELETED");
          break;
        case edit_record::MERGED:
          xmlNewProp(node2, BAD_CAST "status", BAD_CAST "MERGED");
          break; 
        case edit_record::SPLIT:
          xmlNewProp(node2, BAD_CAST "status", BAD_CAST "SPLIT");
          break;
        case edit_record::DUP:
          xmlNewProp(node2, BAD_CAST "status", BAD_CAST "DUP");
          break; 
        default:
          {
            std::cerr <<"No such editing status"<< std::endl;
            exit(0);
          }
        } //switch

        for (unsigned int k = 0; k<rec.replacements_.size(); k++) {
          str = ToString((double)rec.replacements_[k]);
          xmlNewChild(node2, NULL, BAD_CAST "replaced_by", BAD_CAST str.c_str());
        }
      }
    }
  }

  /* 
   * Dumping document to a file
   */
  std::string full_xml_name = file_path+xml_filename;
  xmlSaveFormatFileEnc(full_xml_name.c_str(), doc, MY_ENCODING, 1);

  /*free the document */
  xmlFreeDoc(doc);
  
  /*
   *Free the global variables that may
   *have been allocated by the parser.
   */
  xmlCleanupParser();

  std::string full_label_image_name = file_path+label_image_name;
  marciejSeg.output_cell_pixels( full_label_image_name );
  //generate_label_image(label_image_name, marciejSeg.image_size(), cells);

  return 0;
}

/********************** xml_util_read ***********************************/

void xml_util_read( std::string const& file_path,
                    std::string const& xml_filename, 
                    std::string & image_path,
                    std::string & image_name,
                    maciejSegmentation& marciejSeg )
{
  xmlDocPtr doc; /* the resulting document tree */
  xmlNodePtr root_element, cur_node = NULL;

  LIBXML_TEST_VERSION;
    
  //Parse the resource 
  std::string  xml_full_filename = file_path+xml_filename;
  doc = xmlReadFile(xml_full_filename.c_str(), NULL, 0);
  if (doc == NULL) {
    std::cerr<<"Failed to parse "<<xml_full_filename<<std::endl;
    return;
  }

  /*Get the root element node */
  root_element = xmlDocGetRootElement(doc);

  //Read the parameters
  xmlChar* attribute;
  std::vector<std::string>        parameters;
  std::string label_image_name;

  attribute = xmlGetProp(root_element, BAD_CAST "program");
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "version");
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "image_path");
  image_path = (char*)attribute;
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "image");
  image_name = (char*)attribute;
  xmlFree( attribute );

  /*
  attribute = xmlGetProp(root_element, BAD_CAST "type");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "sstack");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "estack");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "swidth");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  */

  parameters.clear();
  attribute = xmlGetProp(root_element, BAD_CAST "channel");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "threshold");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "grid");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "median_filter_radius");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "morph_opt_radius");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "mode");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "gliadat");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "neurdat");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );

  marciejSeg.set_parameters( parameters );

  attribute = xmlGetProp(root_element, BAD_CAST "label_image");
  label_image_name = (char*)attribute;
  xmlFree( attribute );

  //This may not be the most elegant way of parsing the tree, but it
  //works.
  cur_node =  root_element->children;
  for ( ; cur_node; cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE && 
        !xmlStrcmp(cur_node->name, BAD_CAST "nucleus") ) {
      rich_cell::Pointer cell = read_one_nucleus( cur_node );
      marciejSeg.add_cell( cell );
    }
  }
 
  /*
   * Cleanup function for the XML library.
   */
  xmlCleanupParser();

  /*
   * this is to debug memory for regression tests
   */
  xmlMemoryDump();
  
  // Read in the label image and update the cell pixels
  typedef itk::RGBPixel<unsigned char>      pixelType;
  typedef itk::Image< pixelType,  3 >       imageType;
  typedef itk::ImageFileReader< imageType > readerType;

  std::string label_image_filename = file_path+label_image_name;

  readerType::Pointer reader = readerType::New();
  reader->SetFileName(label_image_filename);
  reader->Update();
  imageType::Pointer label_image = reader->GetOutput();
  marciejSeg.update_cell_pixels(label_image);
}

/***************** xml_util_write_param ********************************/

void xml_util_write_param(std::string const& xml_filename, 
                          std::vector<std::string> const & parameters)
{
  xmlDocPtr doc = NULL;       /* document pointer */
  xmlNodePtr root_node = NULL, node = NULL, node1 = NULL, node2 = NULL, node3 = NULL;/* node pointers */
  xmlDtdPtr dtd = NULL;       /* DTD pointer */

  LIBXML_TEST_VERSION;

  /* 
   * Creates a new document, a node and set it as a root node
   */
  //std::vector<std::string> parameters = marciejSeg.parameters();

  doc = xmlNewDoc(BAD_CAST "1.0");
  root_node = xmlNewNode(NULL, BAD_CAST "Parameters");
  xmlNewProp(root_node, BAD_CAST "program", BAD_CAST "Maciej Segmentation");
  xmlNewProp(root_node, BAD_CAST "version", BAD_CAST "1.0");

  xmlDocSetRootElement(doc, root_node);

  xmlNewChild(root_node, NULL, BAD_CAST "channel", BAD_CAST parameters[0].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "threshold", BAD_CAST parameters[1].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "grid", BAD_CAST parameters[2].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "median_filter_radius", BAD_CAST parameters[3].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "morph_opt_radius", BAD_CAST parameters[4].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "mode", BAD_CAST parameters[5].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "gliadat", BAD_CAST parameters[6].c_str());
  xmlNewChild(root_node, NULL, BAD_CAST "neurdat", BAD_CAST parameters[7].c_str());
  

  //xmlNewProp(root_node, BAD_CAST "image_path", BAD_CAST parameters[0].c_str());
  //xmlNewProp(root_node, BAD_CAST "image", BAD_CAST parameters[1].c_str());
  //xmlNewProp(root_node, BAD_CAST "type", BAD_CAST parameters[2].c_str());
  //xmlNewProp(root_node, BAD_CAST "sstack", BAD_CAST parameters[3].c_str());
  //xmlNewProp(root_node, BAD_CAST "estack", BAD_CAST parameters[4].c_str());
  //xmlNewProp(root_node, BAD_CAST "swidth", BAD_CAST parameters[5].c_str());

  /* 
   * Dumping document to a file
   */
  xmlSaveFormatFileEnc(xml_filename.c_str(), doc, MY_ENCODING, 1);

  /*free the document */
  xmlFreeDoc(doc);
  
  /*
   *Free the global variables that may
   *have been allocated by the parser.
   */
  xmlCleanupParser();
}

/***************** xml_util_read_param ********************************/

void xml_util_read_param(std::string const& xml_filename,
                         std::vector<std::string> & parameters)
{
  xmlDocPtr doc; /* the resulting document tree */
  xmlNodePtr root_node, cur_node = NULL;

  LIBXML_TEST_VERSION;
    
  //Parse the resource 
  doc = xmlReadFile(xml_filename.c_str(), NULL, 0);
  if (doc == NULL) {
    std::cerr<<"Failed to parse "<<xml_filename<<std::endl;
    return;
  }

  /*Get the root element node */
  root_node = xmlDocGetRootElement(doc);

  //Read the parameters
  char *contents;

  /*
  attribute = xmlGetProp(root_element, BAD_CAST "image_path");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "image");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "type");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "sstack");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "estack");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "swidth");
  parameters.push_back( (char*)attribute );
  xmlFree( attribute );
  */

  parameters.resize(8);

  cur_node =  root_node->children;
  for ( ; cur_node; cur_node = cur_node->next) {

    if (cur_node->type != XML_ELEMENT_NODE) continue; 

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "channel") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[0];
	  parameters[0] = contents;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "threshold") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[1];
	  parameters[1] = contents;
      xmlFree( contents );
      continue;
    }
    
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "grid") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[2];
      parameters[2] = contents;
	  xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "median_filter_radius") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[3];
      parameters[3] = contents;
	  xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "morph_opt_radius") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[4];
      parameters[4] = contents;
	  xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "mode") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[5];
      parameters[5] = contents;
	  xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "gliadat") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[6];
      parameters[6] = contents;
	  xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "neurdat") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      //std::stringstream( contents ) >> parameters[7];
      parameters[7] = contents;
	  xmlFree( contents );
      continue;
    }
  }

  /*
   * Cleanup function for the XML library.
   */
  xmlCleanupParser();
  /*
   * this is to debug memory for regression tests
   */
  xmlMemoryDump();
  
}

/******************* xml_util_read_idl *********************************/
/*
void xml_util_read_idl(std::string const& file_path,
                       std::string const& xml_filename, 
                       maciejSegmentation& marciejSeg)
{
  xmlDocPtr doc; //the resulting document tree 
  xmlNodePtr root_element, cur_node = NULL;

  LIBXML_TEST_VERSION;
    
  //Parse the resource 
  std::string  xml_full_filename = file_path+xml_filename;
  doc = xmlReadFile(xml_full_filename.c_str(), NULL, 0);
  if (doc == NULL) {
    std::cerr<<"Failed to parse "<<xml_full_filename<<std::endl;
    return;
  }

  //Get the root element node
  root_element = xmlDocGetRootElement(doc);

  //Read the parameters
  xmlChar* attribute;
  std::vector<std::string>        parameters;
  std::string label_image_name;

  attribute = xmlGetProp(root_element, BAD_CAST "size_X");
  int size_x = atoi((char*)attribute);
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "size_Y");
  int size_y = atoi((char*)attribute);
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "size_Z");
  int size_z = atoi((char*)attribute);
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "Number_Of_Nuclei");
  xmlFree( attribute );
  attribute = xmlGetProp(root_element, BAD_CAST "Segmentation_Output_File_Name");
  label_image_name = file_path + std::string((char*)attribute);
  xmlFree( attribute );
  
  parameters.clear();
  parameters.push_back( "1" ); //channel
  parameters.push_back("500" ); //threshold
  parameters.push_back("4"); //grid
  parameters.push_back("1"); // "median_filter_radius"
  parameters.push_back("1"); //morph_opt_radius
  parameters.push_back("1"); // "mode"
  parameters.push_back("");
  parameters.push_back("");


  marciejSeg.set_parameters( parameters );
  cur_node =  root_element->children;
  for ( ; cur_node; cur_node = cur_node->next) {
    if (!cur_node->type == XML_ELEMENT_NODE) continue;
    
    rich_cell::Pointer cell = rich_cell::New();

    attribute = xmlGetProp(root_element, BAD_CAST "Class_Membership");
    std::stringstream((char*)attribute ) >> cell->class_type_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "ID");
    std::stringstream( (char*)attribute ) >> cell->label_;
    xmlFree( attribute );

    // Center
    attribute = xmlGetProp(root_element, BAD_CAST "X");
    std::stringstream( (char*)attribute ) >> cell->center_[0];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "Y");
    std::stringstream( (char*)attribute ) >> cell->center_[1];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "Z");
    std::stringstream( (char*)attribute ) >> cell->center_[2];
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "bending_energy");
    std::stringstream( (char*)attribute ) >> cell->bending_energy_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "boundary_sharing");
    std::stringstream( (char*)attribute ) >> cell->percent_nbr_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "convexity");
    std::stringstream( (char*)attribute ) >> cell->convexity_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "depth");
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "eccentricity");
    std::stringstream( (char*)attribute ) >> cell->eccentricity_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "eccentricity");
    std::stringstream( (char*)attribute ) >> cell->eccentricity_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "energy_of_intensity_ditribution");
    std::stringstream( (char*)attribute ) >> cell->eng_intensity_dist_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "entropy_of_intensity_ditribution");
    std::stringstream( (char*)attribute ) >> cell->entropy_intensity_dist_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "intensity");
    std::stringstream( (char*)attribute ) >> cell->ave_intensity_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "intensity_ratio");
    std::stringstream( (char*)attribute ) >> cell->bound_ints_ratio_;
    xmlFree( attribute );
    
    attribute = xmlGetProp(root_element, BAD_CAST "intensity_variation");
    std::stringstream( (char*)attribute ) >> cell->intensity_variation_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "interior_gradient");
    std::stringstream( (char*)attribute ) >> cell->vol_grad_;
    xmlFree( attribute );

    // bounding box
    rich_cell::RegionType::IndexType start, end;
    rich_cell::RegionType::SizeType size;
    attribute = xmlGetProp(root_element, BAD_CAST "max_x");
    std::stringstream( (char*)attribute ) >> end[0];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "max_y");
    std::stringstream( (char*)attribute ) >> end[1];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "max_z");
    std::stringstream( (char*)attribute ) >> end[2];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "min_x");
    std::stringstream( (char*)attribute ) >> start[0];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "min_y");
    std::stringstream( (char*)attribute ) >> start[1];
    xmlFree( attribute );
    attribute = xmlGetProp(root_element, BAD_CAST "min_z");
    std::stringstream( (char*)attribute ) >> start[2];
    xmlFree( attribute );
    size[0] = end[0]-start[0]+1;
    size[1] = end[1]-start[1]+1;
    size[2] = end[2]-start[2]+1;
    cell->bounding_box_.SetSize( size );
    cell->bounding_box_.SetIndex( start );

    attribute = xmlGetProp(root_element, BAD_CAST "number_of_poles_on_surface");
    std::stringstream( (char*)attribute ) >> cell->num_poles_on_surface_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "orientation");
    std::stringstream( (char*)attribute ) >> cell->orientation_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "radius_variation");
    std::stringstream( (char*)attribute ) >> cell->radius_variance_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "radius_variation");
    std::stringstream( (char*)attribute ) >> cell->radius_variance_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "shape");
    std::stringstream( (char*)attribute ) >> cell->shape_factor_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "skew_of_intensity_ditribution");
    std::stringstream( (char*)attribute ) >> cell->skew_intensity_dist_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "surface_area");
    std::stringstream( (char*)attribute ) >> cell->surface_area_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "surface_gradient");
    std::stringstream( (char*)attribute ) >> cell->ave_bound_gradient_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "texture");
    std::stringstream( (char*)attribute ) >> cell->texture_;
    xmlFree( attribute );

    attribute = xmlGetProp(root_element, BAD_CAST "volume");
    std::stringstream( (char*)attribute ) >> cell->volume_;
    xmlFree( attribute );

    marciejSeg.add_cell( cell );
  }

  // Read in the label image and update the cell pixels
  marciejSeg.update_cell_pixels(label_image_name, size_x, size_y, size_z);
}
*/

/******************* Local functions ***********************************/

static rich_cell::Pointer read_one_nucleus( xmlNodePtr node )
{
  char *contents;
  xmlChar* attribute;

  rich_cell::Pointer cell = rich_cell::New();
  
  xmlNodePtr cur_node =  node->children;
  for ( ; cur_node; cur_node = cur_node->next) {

    if (cur_node->type != XML_ELEMENT_NODE) continue; 

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "validity") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->valid_;
      xmlFree( contents );
      continue;
    }
    
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "duplicated") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->dup_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "label") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->label_;
      xmlFree( contents );
      continue;
    } 
  
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "volume") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->volume_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "intensity") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->ave_intensity_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "texture") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->texture_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "eccentricity") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->eccentricity_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "average_radius") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->average_radius_;
      xmlFree( contents );
      continue;
    }
    
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "neuronal_signal") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->neuronal_signal_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "class") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->class_type_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "nearest_nbr_dist") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->nearest_nbr_dist_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "nearest_nbr_label") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->nearest_nbr_label_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "score") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->score_;
      xmlFree( contents );
      continue;
    }
    
    /*
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "convexity") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->convexity_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "shape_factor") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->shape_factor_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "bending_energy") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->bending_energy_;
      xmlFree( contents );
      continue;
    }
    
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "boundary_gradient") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->ave_bound_gradient_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "volume_gradient") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->vol_grad_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "radius_variance") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->radius_variance_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "intensity_ratio") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->bound_ints_ratio_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "shared_boundary") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->percent_nbr_;
      xmlFree( contents );
      continue;
    }

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "engergy_intensity_distribution") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->eng_intensity_dist_;
      xmlFree( contents );
      continue;
    }

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "entropy_intensity_distribution") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->entropy_intensity_dist_;
      xmlFree( contents );
      continue;
    }

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "intensity_variation") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->intensity_variation_;
      xmlFree( contents );
      continue;
    }

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "num_of_poles_on_surface") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->num_poles_on_surface_;
      xmlFree( contents );
      continue;
    }

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "skew_intensity_distribution") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->skew_intensity_dist_;
      xmlFree( contents );
      continue;
    }

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "surface_area") ) {
      contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->surface_area_;
      xmlFree( contents );
      continue;
     }
     
     if ( !xmlStrcmp(cur_node->name, BAD_CAST "orientation") ) {
       contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->orientation_;
      xmlFree( contents );
      continue;
     }
    */

     if ( !xmlStrcmp(cur_node->name, BAD_CAST "ave_nnbr_distance") ) {
       contents = (char*)xmlNodeGetContent(cur_node);
      std::stringstream( contents ) >> cell->ave_nnbr_dist_;
      xmlFree( contents );
      continue;
     }

    //center
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "center") ) {
      xmlNodePtr child_node =  cur_node->children;
      for (; child_node; child_node = child_node->next) {
        if (child_node->type == XML_ELEMENT_NODE && 
            !xmlStrcmp(child_node->name, BAD_CAST "x")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> cell->center_[0];
          xmlFree( contents );
        }
        if (child_node->type == XML_ELEMENT_NODE && 
            !xmlStrcmp(child_node->name, BAD_CAST "y")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> cell->center_[1];
          xmlFree( contents );
        }  
        if (child_node->type == XML_ELEMENT_NODE && 
            !xmlStrcmp(child_node->name, BAD_CAST "z")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> cell->center_[2];
          xmlFree( contents );
        }  
      }
      continue;
    }

    //bounding_box
    rich_cell::RegionType::IndexType start;
    rich_cell::RegionType::SizeType size;

    if ( !xmlStrcmp(cur_node->name, BAD_CAST "bounding_box") ) {
      xmlNodePtr child_node =  cur_node->children;
      for (; child_node; child_node = child_node->next) {
        if (child_node->type != XML_ELEMENT_NODE) continue;

        if (!xmlStrcmp(child_node->name, BAD_CAST "min_x")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> start[0];
          xmlFree( contents );
        }
        if (!xmlStrcmp(child_node->name, BAD_CAST "min_y")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> start[1];
          xmlFree( contents );
        }  
        if (!xmlStrcmp(child_node->name, BAD_CAST "min_z")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> start[2];
          xmlFree( contents );
        }  
        if (!xmlStrcmp(child_node->name, BAD_CAST "size_x")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> size[0];
          xmlFree( contents );
        }
        if (!xmlStrcmp(child_node->name, BAD_CAST "size_y")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> size[1];
          xmlFree( contents );
        }  
        if (!xmlStrcmp(child_node->name, BAD_CAST "size_z")) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> size[2];
          xmlFree( contents );
        } 
      }
      cell->bounding_box_.SetSize( size );
      cell->bounding_box_.SetIndex( start );

      continue;
    }

    /*
    //neighbor
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "neighbors") ) {
      xmlNodePtr child_node =  cur_node->children;
      int nb_label;
      for (; child_node; child_node = child_node->next) {
        if (child_node->type == XML_ELEMENT_NODE) {
          contents = (char*)xmlNodeGetContent(child_node);
          std::stringstream( contents ) >> nb_label;
          cell->neighbors_.push_back(nb_label);
          xmlFree( contents );
        }
      }
      continue;
    }
    */

    //Edit history
    if ( !xmlStrcmp(cur_node->name, BAD_CAST "edit_history") ) {
      xmlNodePtr child_node =  cur_node->children;
      for (; child_node; child_node = child_node->next) {
        if (child_node->type == XML_ELEMENT_NODE) {
          edit_record rec;
          attribute = xmlGetProp(child_node, BAD_CAST "date");
          rec.date_ = (char*)attribute;
          xmlFree( attribute );
          attribute = xmlGetProp(child_node, BAD_CAST "name");
          rec.name_ = (char*)attribute;
          xmlFree( attribute );

          attribute = xmlGetProp(child_node, BAD_CAST "status");
          if (!xmlStrcmp(attribute, BAD_CAST "ADDED"))
            rec.status_ = edit_record::ADDED;
          else if (!xmlStrcmp(attribute, BAD_CAST "DELETED"))
            rec.status_ = edit_record::DELETED;
          else if (!xmlStrcmp(attribute, BAD_CAST "MERGED"))
            rec.status_ = edit_record::MERGED;
          else if (!xmlStrcmp(attribute, BAD_CAST "SPLIT"))
            rec.status_ = edit_record::SPLIT;
          else //(!xmlStrcmp(attribute, BAD_CAST "DUP"))
            rec.status_ = edit_record::DUP;
          xmlFree( attribute );

          xmlNodePtr grandchild_node =  child_node->children;
          int replaced_by;
          for (; grandchild_node; grandchild_node = grandchild_node->next) {
            if (grandchild_node->type == XML_ELEMENT_NODE) {
              contents = (char*)xmlNodeGetContent(grandchild_node);
              std::stringstream( contents ) >> replaced_by;
              rec.replacements_.push_back(replaced_by);
              xmlFree( contents );
            } //if
          } //for grandchild
          cell->edit_history_.push_back(rec);
        } //if
      }
      continue;
    }

  }

  return cell;
}
 
#else
int xml_util_write(std::string filename, 
                   std::vector<rich_cell::Pointer> const& cell) 
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
  return 0;
}

void xml_util_read(std::string filename,  
                   std::vector<rich_cell::Pointer>& cells)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}

void xml_util_read_param(std::string xml_filename)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}

void xml_util_write_param(std::string xml_filename)
{
  std::cerr<<"tree support not compiled in\n"<<std::endl;
}
#endif


