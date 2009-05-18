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
	segmentation->LoadClassInfoFromFile(classFile);
	segmentation->WriteToXML(xmlFullName);

	//Now load up the label image and split it into a separate image for each class
	segmentation->SaveLabelByClass();
	
	delete segmentation;

	return 1;
}
