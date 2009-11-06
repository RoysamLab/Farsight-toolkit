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
	if(argc < 4)
	{
		std::cout<<"Usage: compute_nuclei_features <InputImageFileName> <InputLabelFileName> <SegmentationResultsFileName>\n";
		return 0;
	}
    
	std::string imageName = argv[1];
	std::string labelName = argv[2];
	std::string resultsName = argv[3];
	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();	

	segmentation->LoadFromImages(imageName,labelName);
	
	if(argc == 5)
	{
		segmentation->LoadAssociationsFromFile(argv[4]);
	}
	if(argc == 6)
	{
		//segmentation->LoadClassInfoFromFile(argv[5]);
	}

	//segmentation->WriteToXML(resultsName);	
	
	delete segmentation;

	return 1;
}
