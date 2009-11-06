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

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout<<"Usage: convert_dat_to_label <InputImageFileName> <InputdatlFileName>\n";
		return 0;
	}
    
	std::string imageName = argv[1];
	std::string datName = argv[2];
	
	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();	

	std::cout<<" Reading from dat segmentation output...";
	segmentation->LoadFromDAT(imageName,imageName);
	std::cout<<"done"<<std::endl;
	
	//segmentation->SaveLabel();
	
	delete segmentation;

	return 1;
}
