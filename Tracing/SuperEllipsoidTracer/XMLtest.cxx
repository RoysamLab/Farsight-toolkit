// vector::begin
#include <iostream>
#include <vector>
#include "tinyxml.h"
//using namespace std;

int main ()
{
	//std::vector<int> myvector;
 // for (int i=1; i<=5; i++) myvector.push_back(i);

 // std::vector<int>::iterator it;

 // std::cout << "myvector contains:";
 // for ( it=myvector.begin() ; it < myvector.end(); it++ )
 //   std::cout << " " << *it;

 // std::cout << std::endl;

 // return 0;

	TiXmlDocument doc("test.xml");
	if (!doc.LoadFile()) {
		return false;
	}

	TiXmlElement* root = doc.FirstChildElement( "SuperElliposoidTracing" );
	while ( root )		{
		TiXmlAttribute* rAttrib = root->FirstAttribute();
		while (rAttrib)	{
			if (!strcmp(rAttrib->Name(),"InputFileName"))	{
				std::cout << "Input file name: " << rAttrib->Value() <<std::endl;
			}
			else if (!strcmp(rAttrib->Name(),"OutputFileName"))	{
				std::cout << "Output file name: " << rAttrib->Value() <<std::endl;
			}
			rAttrib=rAttrib->Next();
		}

		TiXmlElement* param = root->FirstChildElement( "Parameter" );
		while ( param )	{
			TiXmlAttribute* pAttrib = param->FirstAttribute();
			if (pAttrib){
				if (!strcmp(pAttrib->Name(),"GridSpacing"))	{
					std::cout << "GridSpacing:  " << pAttrib->Value() << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"StepRatio"))	{
					std::cout << "StepRatio:  " << pAttrib->Value() << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"AspectRatio"))	{
					std::cout << "AspectRatio:  " << pAttrib->Value() << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"THRESHOLD"))	{
					std::cout << "THRESHOLD:  " << pAttrib->Value() << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"minContrast"))	{
					std::cout << "minContrast:  " << pAttrib->Value() << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"MaximumVesselWidth"))	{
					std::cout << "MaximumVesselWidth:  " << pAttrib->Value() << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"MinimumVesselWidth"))	{
					std::cout << "MinimumVesselWidth:  " << pAttrib->Value() << std::endl;
				}

				else if(!strcmp(pAttrib->Name(),"MinimumVesselLength"))	{
					std::cout << "MinimumVesselLength:  " << pAttrib->Value() << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"StartTHRESH"))	{
					std::cout << "StartTHRESH:  " << pAttrib->Value() << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"StartIterations"))	{
					int temp = -1;
					if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
						n->ID = temp;
					std::cout << "StartIterations:  " << pAttrib->Value() << std::endl;
				}
				else {
					std::cout << "UNRECOGNIZED TAG:  " << pAttrib->Value() << std::endl;
				}
			}
			param = param->NextSiblingElement();
		}
		root = root->NextSiblingElement();
	}
}