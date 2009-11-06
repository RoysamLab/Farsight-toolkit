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

#include "ftkNuclearSegmentation.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout<<"Usage: classify_nuclei <xmlResultFile> <ClassFile>\n";
		return 0;
	}

	std::string xmlFullName = argv[1];

	size_t found = xmlFullName.find_last_of(".");
	std::string xmlBaseName = xmlFullName.substr(0,found);

	std::string classFile = argv[2];

	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();
	segmentation->RestoreFromXML(xmlFullName);
	//segmentation->LoadClassInfoFromFile(classFile);
	//segmentation->WriteToXML(xmlFullName);

	//Now load up the label image and split it into a separate image for each class
	segmentation->SaveLabelByClass();
	
	delete segmentation;

	return 1;
}
