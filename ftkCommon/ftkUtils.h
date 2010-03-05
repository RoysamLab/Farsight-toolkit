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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkUtils_h
#define __ftkUtils_h

#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>
#include <vtkDoubleArray.h>

#include "ftkImage/ftkImage.h"
#include <tinyxml/tinyxml.h>

#include <string>
#include <ctime>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace ftk
{

bool FileExists(std::string filename);
bool AppendTextFile(std::string filename, std::string text);			//Add new line to the file with the given text
bool SaveTable(std::string filename, vtkSmartPointer<vtkTable> table);
vtkSmartPointer<vtkTable> LoadTable(std::string filename);
bool SaveXMLImage(std::string filename, ftk::Image::Pointer image);
ftk::Image::Pointer LoadXMLImage(std::string filename);
std::string NumToString(double d);
std::string NumToString(int i);
std::string NumToString(double d, int p);
std::string TimeStamp();
std::string GetExtension(std::string filename);
std::string SetExtension(std::string filename, std::string ext); // if ext = "", remove extension (and .)

}  // end namespace ftk

#endif	// end __ftkUtils_h
