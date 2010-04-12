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
#ifndef PROJECTMANAGER_H_
#define PROJECTMANAGER_H_

#include "tinyxml/tinyxml.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

namespace ftk
{

#define MY_ENCODING "ISO-8859-1"
//valid file types are
//	image:	tiff pic mhd
//	soma:	tiff pic 
//	trace:	swc vtk rpi.xml

struct FileInfoManager
{
	std::string fileName;
	std::string fileType;
	double tx, ty, tz;
};

class ProjectManager
{
public:
	ProjectManager(const char * filename);
	ProjectManager();
	void addFile(std::string fileName, std::string fileType, double x, double y, double z);
	void addOutputTraceFile(unsigned int i, std::string fileName);
	void readProject(const char * filename);
	bool writeProject(const char* filename);
	unsigned int size();
	std::string GetFileName(int i);
	std::string GetFileType(int i);
	double GetTranslationX(int i);
	double GetTranslationY(int i);
	double GetTranslationZ(int i);

	void ReplaceTranslations(std::string fileName, bool zOnly=false);

private:
	std::vector<FileInfoManager> fileInfo;
};

} // end namespace ftk

#endif