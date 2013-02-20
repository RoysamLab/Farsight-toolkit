#include "ImageTransformData.h"
#include <iostream>
#include <tinyxml/tinyxml.h>
#include <vector>



class MontagingToTraceEditorXMLConverter
{
public:
	MontagingToTraceEditorXMLConverter(std::string fileName);
	~MontagingToTraceEditorXMLConverter();

	void readXMLFile();
	void writeXMLFile();
private:
	std::vector<ImageTransformData*> image_transforms_data;
	std::string fileName;
	
};