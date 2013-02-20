#include "MontagingToTraceEditorXMLConverter.h"

MontagingToTraceEditorXMLConverter::MontagingToTraceEditorXMLConverter(std::string fileName)
{
	this->fileName = fileName;
}

void MontagingToTraceEditorXMLConverter::readXMLFile()
{
	TiXmlDocument doc;	
	
	if (!doc.LoadFile(fileName.c_str()))
		std::cerr << "Unable to load XML File" << std::endl;
	
	TiXmlElement* rootElement = doc.RootElement();
	TiXmlNode* transform_node = rootElement->FirstChild()->NextSibling()->NextSibling(); 

	while (transform_node != NULL)
	{
		TiXmlNode* to_image_ID_node = transform_node->FirstChild()->NextSibling();	//2nd child is to_image_ID
		
		TiXmlNode* parameters_node = transform_node->FirstChild()->NextSibling()->NextSibling()->NextSibling()->NextSibling(); //5th child is parameters_node
		TiXmlElement* parameters_element = parameters_node->ToElement();

		std::string fileName = to_image_ID_node->FirstChild()->Value();

		std::cout << fileName << " ";

		double tX, tY, tZ;

		parameters_element->QueryDoubleAttribute("tx", &tX);
		parameters_element->QueryDoubleAttribute("ty", &tY);
		parameters_element->QueryDoubleAttribute("tz", &tZ);

		std::cout << tX << " " << tY << " " << tZ << std::endl;

		image_transforms_data.push_back(new ImageTransformData(fileName, tX, tY, tZ));

		transform_node = transform_node->NextSibling();
	}
}

void MontagingToTraceEditorXMLConverter::writeXMLFile()
{
	TiXmlDocument doc;
	TiXmlDeclaration* decl = new TiXmlDeclaration("1.0", "", "");
	doc.LinkEndChild(decl);

	TiXmlElement* root = new TiXmlElement("Source");
	doc.LinkEndChild(root);

	std::vector<ImageTransformData*>::iterator image_transforms_data_iter;

	int counter = 0;
	for (image_transforms_data_iter = image_transforms_data.begin(); image_transforms_data_iter != image_transforms_data.end(); image_transforms_data_iter++)
	{
		TiXmlElement *File = new TiXmlElement("File");
		root->LinkEndChild(File);

		ImageTransformData *image_transform = *image_transforms_data_iter;

		File->SetDoubleAttribute("tX", -1.0 * (*image_transforms_data_iter)->gettX());
		File->SetDoubleAttribute("tY", -1.0 * (*image_transforms_data_iter)->gettY());
		File->SetDoubleAttribute("tZ", -1.0 * (*image_transforms_data_iter)->gettZ());
		File->SetAttribute("FileName", (*image_transforms_data_iter)->getFileName());
		File->SetAttribute("Type", "Image");
	}

	doc.SaveFile("project.xml");
}