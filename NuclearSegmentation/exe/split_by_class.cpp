#include "ftkNuclearSegmentation.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 4)
	{
		std::cout<<"Usage: classify_nuclei <ProjectPath> <xmlResultFile> <ClassFile>\n";
		return 0;
	}
    
	std::string projPath = argv[1];

	std::string xmlFullName = argv[2];

	size_t found = xmlFullName.find_last_of(".");
	std::string xmlBaseName = xmlFullName.substr(0,found);

	std::string classFile = argv[3];

	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation(projPath,xmlBaseName);
	segmentation->RestoreFromXML();
	segmentation->LoadClassInfoFromFile(classFile);
	segmentation->WriteToXML();
	
	delete segmentation;

	return 1;
}
