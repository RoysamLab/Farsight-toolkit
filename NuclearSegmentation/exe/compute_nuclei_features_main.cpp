#include "ftkNuclearSegmentation.h"

#include <iostream>

using namespace std;

int main(unsigned int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout<<"Usage: compute_nuclei_features <ProjectPath> <InputImageFileName> <SegmentationResultsFileName>\n";
		return 0;
	}	
    
	std::string projectPath = argv[1];
	std::string imageName = argv[2];
	std::string resultsName = argv[3];
	std::string projectName = imageName.substr(0,imageName.find_first_of("."));
	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation(projectPath,projectName);	
	segmentation->LoadFromResult(imageName.c_str(), resultsName.c_str());	
	segmentation->SaveAll();	
	
	delete segmentation;

	return 1;
}
