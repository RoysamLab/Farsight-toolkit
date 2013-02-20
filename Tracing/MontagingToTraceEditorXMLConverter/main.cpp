#include "MontagingToTraceEditorXMLConverter.h"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "No input file specified" << std::endl;
		return -1;
	}

	MontagingToTraceEditorXMLConverter *converter = new MontagingToTraceEditorXMLConverter(argv[1]);
	converter->readXMLFile();
	converter->writeXMLFile();

	return 0;
}