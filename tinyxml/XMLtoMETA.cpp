/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "tinyxml.h"
#include <fstream>	
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

bool header_dumped;

int dump_header(TiXmlElement* pElement, const char* outFilename) 
{
	if ( !pElement ) return 0;

	TiXmlAttribute* pAttrib=pElement->FirstAttribute();

	int i=0;

	int s = strlen(outFilename);
	char* headFilename = (char *)malloc(s+9);
	strncpy(headFilename,outFilename,s-4);
	headFilename[s-4] = '\0';
	strcat(headFilename,"_sel_lbls.txt");
	
	ofstream outFile;
	outFile.open(headFilename, ios::out | ios::trunc );
	if (!outFile.is_open() ) return 0;
	
	while (pAttrib)
	{
		if ( strstr(pAttrib->Name(),"ID") )
		{
		}
		else if ( strstr(pAttrib->Name(),"Class_Membership") )
		{
		}
		else	//Regular attribute just print it
		{
			outFile << pAttrib->Name() << endl;
			i++;
		}
		pAttrib=pAttrib->Next();
	}

	//Now print y and ID
	outFile << "RESPONSE" << endl;
	i++;
	outFile << "ID" << endl;
	i++;

	outFile.close();
	header_dumped = true;

	return i;
}

int dump_attribs_to_file(TiXmlElement* pElement, const char* outFilename)
{
	if ( !pElement ) return 0;

	TiXmlAttribute* pAttrib=pElement->FirstAttribute();

	int i=0;
	int ival;
	double dval;

	int ID; //The ID of the object
	int y;  //This is the class/label of the object

	ofstream outFile;
	outFile.open(outFilename, ios::out | ios::app );

	if (!outFile.is_open() ) return 0;

	outFile << " ";
	long startPos = outFile.tellp(); //position at beginning of row
	long pos;				//for holding current position later on
	long neededSpaces;		

	while (pAttrib)
	{
		//*outFilePtr << pIndent << pAttrib->Name() << ": value = [" << pAttrib->Value() << "]";
		if ( strstr(pAttrib->Name(),"ID") )
		{
			if( pAttrib->QueryIntValue(&ival) == TIXML_SUCCESS ) 
				ID = ival;
		}
		else if ( strstr(pAttrib->Name(),"Class_Membership") )
		{
			if ( pAttrib->QueryIntValue(&ival) == TIXML_SUCCESS ) //Check for int
			{
				y = ival;
			}
			else //must be string, so check for names is this field
			{
				const char *nm = pAttrib->Value();
				if ( strstr(nm, "Neurons") )
				{
					y = 1;
				}
				else if ( strstr(nm, "Astrocytes") )
				{
					y = 2;
				}
				else if ( strstr(nm, "Microglia") )
				{
					y = 3;
				}
				else if ( strstr(nm, "Endothelials") )
				{
					y = 4;
				}
			}
		}
		else	//Regular attribute just print it
		{
			pos = outFile.tellp();
			neededSpaces = (startPos + (i*15)) - pos;
			if (neededSpaces)
			{
				for (int i=0; i<neededSpaces; ++i)
				{
					outFile << " ";
				}
			}

			if (pAttrib->QueryDoubleValue(&dval)==TIXML_SUCCESS) 
			{
				outFile << setprecision(10) << dval;
			}
			else if (pAttrib->QueryIntValue(&ival)==TIXML_SUCCESS)
			{
				outFile << ival;
			}
			else
			{
				//outFile << pAttrib->Value();
				outFile << "-1";
			}
			i++;
		}

		pAttrib=pAttrib->Next();
	}

	//Now print y and ID
	pos = outFile.tellp();
	neededSpaces = (startPos + (i*15)) - pos;
	if (neededSpaces)
	{
		for (int i=0; i<neededSpaces; ++i)
		{
			outFile << " ";
		}
	}
	outFile << y;
	i++;

	pos = outFile.tellp();
	neededSpaces = (startPos + (i*15)) - pos;
	if (neededSpaces)
	{
		for (int i=0; i<neededSpaces; ++i)
		{
			outFile << " ";
		}
	}
	outFile << ID;
	i++;

	outFile << endl;
	outFile.close();

	return i;
}

void dump_to_file ( TiXmlNode* pParent, const char* outFilename )
{
	if ( !pParent ) return;

	TiXmlNode* pChild;
	int t = pParent->Type();

	switch ( t )
	{
		case TiXmlNode::ELEMENT:
			if ( strstr(pParent->Value(), "Nuclear_Features") )
			{
				dump_attribs_to_file(pParent->ToElement(), outFilename);
				if (!header_dumped)
				{
					dump_header(pParent->ToElement(), outFilename);
				}
			}
		break;

		default:
		break;
	}
	
	for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) 
	{
		dump_to_file( pChild, outFilename);
	}
}

//Attempt to open both documents
void dump_to_file(const char* xmlFilename)
{
	TiXmlDocument doc(xmlFilename);
	bool loadOkay = doc.LoadFile();

	char * outFilename = (char *)xmlFilename;
	char * pch = strstr (outFilename,".xml");
	strncpy(pch,".txt",4);
	
	ofstream outFile;
	outFile.open(outFilename, ios::out | ios::trunc );
	if (!outFile.is_open() )
	{
		loadOkay = false;
	}
	outFile.close();
	
	if ( loadOkay)
	{
		std::cerr << "Loaded Document: " << xmlFilename << std::endl;
		header_dumped = false;
		dump_to_file( &doc, outFilename );
	}
	else
	{
		std::cerr << "Failed to Load Document: " << xmlFilename << std::endl;
	}
}

bool isXML(const char* filename)
{
	if ( !strstr(filename, ".xml") )
	{
		return false;
	}
	else
	{
		return true;
	}
}

// ----------------------------------------------------------------------
// main() for printing file named on the command line
// ----------------------------------------------------------------------
int main(int argc, char* argv[])
{
	for (int i=1; i<argc; i++)
	{
		if ( isXML(argv[i]) )
			dump_to_file(argv[i]);
	}
	return 0;
}
