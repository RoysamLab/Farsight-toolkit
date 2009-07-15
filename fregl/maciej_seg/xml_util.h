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

/** @file xml_util.h
*   @brief utility functions for xml
*
*   @author Chia-Ling Tsai
*/

#ifndef __XML_UTIL_H_
#define __XML_UTIL_H_

#include <string>
#include <iostream>
#include <sstream>

#include "rich_cell.h"
#include "maciejSegmentation.h"

//: Write the cells to an xml file specified by the filename
bool xml_util_write(std::string const& file_path,
                    std::string const& xml_filename, 
                    std::string const& image_path,
                    std::string const& image_name,
                    maciejSegmentation const& marciejSeg);

//: Read in the xml file of the segmentation result back to maciejSegmentation. 
void xml_util_read( std::string const& file_path,
                    std::string const& xml_filename, 
                    std::string & image_path,
                    std::string & image_name,
                    maciejSegmentation& marciejSeg );

//: Read the parameter xml file (path included in the filename)
void xml_util_read_param(std::string const& xml_filename, 
                         std::vector<std::string> & parameters);

//: Write the parameter xml file (path included in the filename)
void xml_util_write_param(std::string const& xml_filename,
                          std::vector<std::string> const & parameters);

//: Read the result writen by the IDL farsight 
void xml_util_read_idl(std::string const& file_path,
                       std::string const& xml_filename, 
                       maciejSegmentation& marciejSeg);

inline std::string ToString(double val)
{
    std::ostringstream strm;
    strm<< val<<std::endl;
    return strm.str();
}
#endif
