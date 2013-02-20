#include "TraceObject.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
#include <string>
#include <limits>
#include <tinyxml/tinyxml.h>
#include "vul/vul_file.h"
#include "ftkCommon/ftkUtils.h"
#include "VolumeOfInterest.h"

int main(int argc, char* argv[])
{

	if (argc < 3)
	{
		std::cout << "Must run program with correct number of arguments" << std::endl;
		std::cout << "Usage: " << argv[0] << " project_xml_file.xml vtkpolydata.vtp" << std::endl;
		return -1;
	}

	char* xml_file = argv[1];
	TiXmlDocument project_doc(argv[1]);

	vcl_string dir_path = vul_file::dirname(xml_file);

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

	TraceObject *tobj = new TraceObject();
		
	//For each line in the XML file
	while (current_element != 0)
	{
		//File lines
		std::string fileName, type;
		int tX, tY, tZ;

		
		type = current_element->Attribute("Type");
		//std::cout << "Type: " << type << std::endl;

		if (type == "Trace")
		{
			fileName = dir_path + "/" + vul_file::strip_directory(current_element->Attribute("FileName"));
			current_element->QueryIntAttribute("tX", &tX);
			current_element->QueryIntAttribute("tY", &tY);
			current_element->QueryIntAttribute("tZ", &tZ);
			//std::cout << fileName << " " << tX << " " << tY << " " << tZ << std::endl;			
					
			char *swc_file = (char *)fileName.c_str();	//converting const char * to char *
			std::cout << "Reading SWC file: " << swc_file << std::endl;
			
			//Read into the TraceObject from swc_file
			tobj->SetTraceOffset(tX, tY, tZ);
			tobj->ReadFromSWCFile(swc_file);		//DANGEROUS, ReadFromSWCFile will change fileName
			
			//reset fileName and swc_file
			fileName = dir_path + "/" + vul_file::strip_directory(current_element->Attribute("FileName"));
			swc_file = (char *)fileName.c_str();
			//std::cout << swc_file << std::endl;
			
		}	
		current_element = current_element->NextSiblingElement();
	}
	CellTraceModel* cells = new CellTraceModel(tobj->CalculateCellFeatures());
	
        VolumeOfInterest *voi = new VolumeOfInterest();
        voi->ReadVTPVOI(argv[2]);
        voi->CalculateCellDistanceToVOI(cells);

	cells->setCells(tobj->CalculateCellFeatures());

	ftk::SaveTable(vul_file::strip_extension(argv[1]) + ".txt", cells->getDataTable());
}
