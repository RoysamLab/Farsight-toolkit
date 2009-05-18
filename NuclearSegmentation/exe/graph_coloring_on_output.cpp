#include "ftkNuclearSegmentation.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout<<"Usage: graph_coloring_on_output <LabelFileName> <OutputfileName>\n";
		return 0;
	}
    
	std::string segmentationFileName = argv[1];
	std::string outputfileName = argv[2];

	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();	
	segmentation->RunGraphColoring(segmentationFileName, outputfileName);

	delete segmentation;

	return 1;
}
