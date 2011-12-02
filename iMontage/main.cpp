#include "TraceObject.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
#include <string>
#include <limits>
#include <tinyxml/tinyxml.h>

int main(int argc, char* argv[])
{
	TraceObject* tobj = new TraceObject();

	if (argc < 2)
	{
		std::cout << "Must run program with correct number of arguments" << std::endl;
		std::cout << "Usage: " << argv[0] << " project_xml_file" << std::endl;
		return -1;
	}
	
	char* xml_file = argv[1];
	TiXmlDocument project_doc(argv[1]);
	
	//Handles if the file does not exist or is not a valid XML file
	if (!project_doc.LoadFile())
	{
		std::cerr << "Could not load file: " << argv[1] << std::endl;
		return -1;
	}

	TiXmlElement* current_element = project_doc.FirstChildElement();

	if (!current_element)
	{
		std::cerr << "No root element! Is this a valid XML file?" << std::endl;
		return -2;
	}

	TiXmlHandle root = TiXmlHandle(current_element);
	current_element = root.FirstChild().Element();
	std::cout << current_element->Value() << std::endl;


	for (; current_element->NextSibling() != 0; current_element = current_element->NextSiblingElement())
	{
		//File lines
		std::string fileName, type;
		int tX, tY, tZ;
		
		type = current_element->Attribute("Type");
	
		if (type == "Image")
		{
			fileName = current_element->Attribute("FileName");
			current_element->QueryIntAttribute("tX", &tX);
			current_element->QueryIntAttribute("tY", &tY);
			current_element->QueryIntAttribute("tZ", &tZ);
			std::cout << fileName << " " << tX << " " << tY << " " << tZ << std::endl;
		}
	}

	
/*

	//char* swc_file = argv[1];
	std::cout << "Reading SWC file: " << swc_file << std::endl;
	
	tobj->ReadFromSWCFile(swc_file);
	CellTraceModel* cells = new CellTraceModel(tobj->CalculateCellFeatures());

	double maxX = 0;
	double minX = std::numeric_limits<double>::max();
	double maxY = 0;
	double minY = std::numeric_limits<double>::max();
	double maxZ = 0;
	double minZ = std::numeric_limits<double>::max();


	for (int k = 0; k < cells->getCellCount(); k++)
	{
		CellTrace* cell = cells->GetCellAt(k);
		
		double bounds[6];
		cell->getCellBounds(bounds);

		std::cout << "Cell " << k << ": Min X: " << bounds[0] << " Max X: " << bounds[1] << " Min Y: " << bounds[2] << " Max Y: " << bounds[3] << " Min Z: " << bounds[4] << " Max Z: " << bounds[5] << std::endl;
	
		if (bounds[1] > maxX)
			maxX = bounds[1];
		if (minX > bounds[0])
			minX = bounds[0];

		if (bounds[3] > maxY)
			maxY = bounds[3];
		if (minY > bounds[2])
			minY = bounds[2];

		if (bounds[5] > maxZ)
			maxZ = bounds[5];
		if (minZ > bounds[4])
			minZ = bounds[4];
	}

	std::cout << "The extreme bounds of the cells in the image are: (" << minX << ", " << minY << ", " << minZ << ") to (" << maxX << ", " << maxY << ", " << maxZ << ")" << std::endl; 
	
	int image_x_size = maxX - minX + 1;
	int image_y_size = maxY - minY + 1;
	int image_z_size = maxZ - minZ + 1;

	for (int k = 0; k < cells->getCellCount(); k++)
	{		
		CellTrace* cell = cells->GetCellAt(k);
		
		double bounds[6];
		cell->getCellBounds(bounds);
		
		if (bounds[0] < image_x_size * 0.1)	//if the left side of the traces falls inside the 10% on the left of the image, it is most likely in the left overlap region
		{
			//std::cout << "Cell " << k << " rejected for overlap on the left side of the image" << std::endl;
			continue;				//for now do nothing
		}
		
		if (bounds[1] > image_x_size * 0.9)	//if the right side of the traces falls inside the 10% on the right of the image, it is most likely in the right overlap region
		{
			//std::cout << "Cell " << k << " rejected for overlap on the right side of the image" << std::endl;
			continue;				//for now do nothing
		}

		if (bounds[2] < image_y_size * 0.1)	//if the top of the traces falls inside the 10% on the top of the image, it is most likely in the top overlap region
		{
			//std::cout << "Cell " << k << " rejected for overlap on the top of the image" << std::endl;
			continue;
		}

		if (bounds[3] > image_y_size * 0.9)	//if the bottom of the traces falls inside the 10% on the bottom of the image, it is most likely in the bottom overlap region
		{
			//std::cout << "Cell " << k << " rejected for overlap on the bottom of the image" << std::endl;
			continue;
		}

		//Note that we do not have Z-overlap since we do not register images in that dimension. However, if we do, it is a simple matter of extending the above if statements to include those cases in the z-dimension 		
		std::cout << "Cell " << k << " found inside non-overlap region" << std::endl;
		
		std::ostringstream fileNameStream;
		fileNameStream << k << ".swc";

		tobj->WriteToSWCFile(cell->getSegments(), fileNameStream.str().c_str());	
	}*/	
}
