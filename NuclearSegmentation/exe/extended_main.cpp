#include "nuclear_segmentation/yousef_core/yousef_seg.h"
#include "nuclear_segmentation/gui/nuclei.h"

#include <iostream>

using namespace std;

int main(unsigned int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout<<"Usage: nucseg2 <InputImageFileName> <ParametersFileName>\n";
		return 0;
	}
	
	FTKAbstractSegmentation *segmentation = new nuclei();
	segmentation->setup(argv[1], argv[2]);
	segmentation->executeModule(0);	
	segmentation->executeModule(1);
	segmentation->executeModule(2);
	segmentation->save();
	segmentation->generateXML();

	segmentation->initMetaNeural();

	delete segmentation;

	return 1;
}
