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
#ifndef TRACEPROJECTMANAGER_H_
#define TRACEPROJECTMANAGER_H_

#include "tinyxml/tinyxml.h"
#include <string>
#include <vector>

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
	ProjectManager(char * filename);
	bool writeProject(char* filename);
private:
	std::vector<FileInfoManager> fileInfo;
};
#endif;