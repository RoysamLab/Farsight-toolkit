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
		segmentation->LoadClassInfoFromFile(argv[5]);
	}

	segmentation->WriteToXML(resultsName);	
	
	delete segmentation;

	return 1;
}
