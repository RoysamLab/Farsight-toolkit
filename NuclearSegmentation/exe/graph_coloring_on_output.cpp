#include "ftkNuclearSegmentation.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 4)
	{
		std::cout<<"Usage: graph_coloring_on_output <ProjectPath> <SegmentationFileName> <OutputfileName>\n";
		return 0;
	}
    
	std::string projPath = argv[1];
	std::string segmentationFileName = argv[2];
	std::string outputfileName = argv[3];

	size_t found = segmentationFileName.find_last_of(".");
	std::string projName = segmentationFileName.substr(0,found);
	

	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation(projPath,projName);	
	segmentation->AddResultFile(segmentationFileName);

	//run graph coloring
	segmentation->RunGraphColoring(outputfileName.c_str());

		
	delete segmentation;

	return 1;
}
