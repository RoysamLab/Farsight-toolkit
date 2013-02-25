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
#include <ftkFeatures/ftkObjectAssociation.h>

#include <string>
#include <ctime>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace ftk
{

typedef struct { int number; std::string name; std::string type; } Channel;
std::vector<Channel> ReadChannels(TiXmlElement * inputElement);
bool FileExists(std::string filename);
bool AppendTextFile(std::string filename, std::string text);			//Add new line to the file with the given text
bool SaveTable(std::string filename, vtkSmartPointer<vtkTable> table);
bool SaveTableAppend(std::string filename, vtkSmartPointer<vtkTable> table, int id);
bool SaveTableSeries(std::string filename,std::vector< vtkSmartPointer<vtkTable> >  table4DImage,std::string path);
bool SaveImageSeries(std::string seriesfilename, ftk::Image::Pointer image,std::string path);
bool SaveLabelSeries(std::string seriesfilename, ftk::Image::Pointer image,std::string path);

vtkSmartPointer<vtkTable> LoadTable(std::string filename);
vtkSmartPointer<vtkTable> LoadRotatedTable(std::string filename);
vtkSmartPointer<vtkTable> LoadXYZTable(std::string filename);
vtkSmartPointer<vtkTable> AppendLoadTable(std::string filename, vtkSmartPointer<vtkTable> table , double tx, double ty, double tz);
bool SaveXMLImage(std::string filename, ftk::Image::Pointer image);
ftk::Image::Pointer LoadXMLImage(std::string filename);
ftk::Image::Pointer LoadImageSeries(std::string filename);
std::vector<std::string> GetSeriesPaths(std::string filename);
ftk::Image::Pointer LoadImageSeriesLabels(std::string filename);
vtkSmartPointer<vtkTable> AppendTables(vtkSmartPointer<vtkTable> table_initial,vtkSmartPointer<vtkTable> table_new );
std::vector< vtkSmartPointer<vtkTable> > LoadTableSeries(std::string filename);
std::string NumToString(double d);
std::string NumToString(int i);
std::string NumToString(double d, int p);
std::string TimeStamp();
std::string GetExtension(std::string filename);
std::string SetExtension(std::string filename, std::string ext); // if ext = "", remove extension (and .)
std::string GetFilenameFromFullPath(std::string f);
std::string GetFilePath(std::string f);
std::vector<std::string> GetColumsWithString( std::string colName, vtkSmartPointer<vtkTable> table );
std::string GetStringInCaps( std::string in_srting );
bool Load(std::string filename);
vtkSmartPointer<vtkTable> CopyTable(vtkSmartPointer<vtkTable> featureTable );

bool SaveTableSeriesActive(std::string filename,std::vector< vtkSmartPointer<vtkTable> >  table4DImage);
double GetMean(std::vector<double> data);
double GetStd(std::vector<double> data);

//std::vector<ftk::AssociationRule> ReadAssociationRules(TiXmlElement * inputElement);

typedef struct { std::string regionChannelName; std::string targetChannelName; int mode;
				 std::string outputFilename; int radius; int erodeRadius; } PixelAnalysisDefinitions;
}  // end namespace ftk

#endif	// end __ftkUtils_h
