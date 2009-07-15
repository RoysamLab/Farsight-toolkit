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

// tutorial demo program
//#include <stdafx.h>
#include "tinyxml.h"
#include <fstream>		
using namespace std;

// ----------------------------------------------------------------------
// STDOUT dump and indenting utility functions
// ----------------------------------------------------------------------
const unsigned int NUM_INDENTS_PER_SPACE=2;

const char * getIndent( unsigned int numIndents )
{
	static const char * pINDENT="                                      + ";
	static const unsigned int LENGTH=strlen( pINDENT );
	unsigned int n=numIndents*NUM_INDENTS_PER_SPACE;
	if ( n > LENGTH ) n = LENGTH;

	return &pINDENT[ LENGTH-n ];
}

// same as getIndent but no "+" at the end
const char * getIndentAlt( unsigned int numIndents )
{
	static const char * pINDENT="                                        ";
	static const unsigned int LENGTH=strlen( pINDENT );
	unsigned int n=numIndents*NUM_INDENTS_PER_SPACE;
	if ( n > LENGTH ) n = LENGTH;

	return &pINDENT[ LENGTH-n ];
}

int dump_attribs_to_file(TiXmlElement* pElement, ofstream* outFilePtr, unsigned int indent)
{
	if ( !pElement ) return 0;

	TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	int i=0;
	int ival;
	double dval;
	const char* pIndent=getIndent(indent);
	*outFilePtr << endl;
	while (pAttrib)
	{
		*outFilePtr << pIndent << pAttrib->Name() << ": value = [" << pAttrib->Value() << "]";
		//printf( "%s%s: value=[%s]", pIndent, pAttrib->Name(), pAttrib->Value());

		if (pAttrib->QueryIntValue(&ival)==TIXML_SUCCESS)    *outFilePtr << "int=" << ival;
		if (pAttrib->QueryDoubleValue(&dval)==TIXML_SUCCESS) *outFilePtr << " d=" << dval;
		*outFilePtr << endl;
		i++;
		pAttrib=pAttrib->Next();
	}
	return i;
}

int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent)
{
	if ( !pElement ) return 0;

	TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	int i=0;
	int ival;
	double dval;
	const char* pIndent=getIndent(indent);
	printf("\n");
	while (pAttrib)
	{
		printf( "%s%s: value=[%s]", pIndent, pAttrib->Name(), pAttrib->Value());

		if (pAttrib->QueryIntValue(&ival)==TIXML_SUCCESS)    printf( " int=%d", ival);
		if (pAttrib->QueryDoubleValue(&dval)==TIXML_SUCCESS) printf( " d=%1.1f", dval);
		printf( "\n" );
		i++;
		pAttrib=pAttrib->Next();
	}
	return i;	
}

void dump_to_file ( TiXmlNode* pParent, ofstream* outFilePtr, unsigned int indent = 0)
{
	if ( !pParent ) return;

	TiXmlNode* pChild;
	TiXmlText* pText;
	int t = pParent->Type();
	*outFilePtr << getIndent(indent);
	int num;

	switch ( t )
	{
		case TiXmlNode::DOCUMENT:
			*outFilePtr << "Document";
		break;

		case TiXmlNode::ELEMENT:
			*outFilePtr << "Element [" << pParent->Value() << "]";
			num=dump_attribs_to_file(pParent->ToElement(), outFilePtr, indent+1);
			switch(num)
			{
				//case 0:  printf( " (No attributes)"); break;
				case 0:  *outFilePtr << " (No attributs)"; break;
				//case 1:  printf( "%s1 attribute", getIndentAlt(indent)); break;
				case 1:  *outFilePtr << getIndentAlt(indent) << "1 attribute"; break;
				//default: printf( "%s%d attributes", getIndentAlt(indent), num); break;
				default: *outFilePtr << getIndentAlt(indent) << num << " attributes"; break;
			}
		break;

		case TiXmlNode::COMMENT:
			//printf( "Comment: [%s]", pParent->Value());
			*outFilePtr << "Comment: [" << pParent->Value() << "]";
		break;

		case TiXmlNode::UNKNOWN:
			//printf( "Unknown" );
			*outFilePtr << "Unknown";
		break;

		case TiXmlNode::TEXT:
			pText = pParent->ToText();
			//printf( "Text: [%s]", pText->Value() );
			*outFilePtr << "Text: [" << pText->Value() << "]";
		break;

		case TiXmlNode::DECLARATION:
			//printf( "Declaration" );
			*outFilePtr << "Declaration";
		break;

		default:
		break;
	}

	//printf( "\n" );
	*outFilePtr << endl;
	
	for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) 
	{
		dump_to_file( pChild, outFilePtr, indent+1 );
	}


}
void dump_to_stdout( TiXmlNode* pParent, unsigned int indent = 0 )
{
	if ( !pParent ) return;

	TiXmlNode* pChild;
	TiXmlText* pText;
	int t = pParent->Type();
	printf( "%s", getIndent(indent));
	int num;

	switch ( t )
	{
	case TiXmlNode::DOCUMENT:
		printf( "Document" );
		break;

	case TiXmlNode::ELEMENT:
		printf( "Element [%s]", pParent->Value() );
		num=dump_attribs_to_stdout(pParent->ToElement(), indent+1);
		switch(num)
		{
			case 0:  printf( " (No attributes)"); break;
			case 1:  printf( "%s1 attribute", getIndentAlt(indent)); break;
			default: printf( "%s%d attributes", getIndentAlt(indent), num); break;
		}
		break;

	case TiXmlNode::COMMENT:
		printf( "Comment: [%s]", pParent->Value());
		break;

	case TiXmlNode::UNKNOWN:
		printf( "Unknown" );
		break;

	case TiXmlNode::TEXT:
		pText = pParent->ToText();
		printf( "Text: [%s]", pText->Value() );
		break;

	case TiXmlNode::DECLARATION:
		printf( "Declaration" );
		break;
	default:
		break;
	}
	printf( "\n" );
	for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) 
	{
		dump_to_stdout( pChild, indent+1 );
	}
}

// load the named file and dump its structure to STDOUT
void dump_to_stdout(const char* pFilename)
{
	TiXmlDocument doc(pFilename);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		printf("\n%s:\n", pFilename);
		dump_to_stdout( &doc ); // defined later in the tutorial
	}
	else
	{
		printf("Failed to load file \"%s\"\n", pFilename);
	}
}

void dump_to_file(const char* xmlFilename, const char* txtFilename)
{
	TiXmlDocument doc(xmlFilename);
	bool loadOkay = doc.LoadFile();
	ofstream outFile; 
	outFile.open(txtFilename, ios::out | ios::trunc );
	if ( loadOkay && outFile.is_open() )
	{
		std::cerr << "Loaded Document: " << xmlFilename << std::endl;
		dump_to_file( &doc, &outFile );
	}
	else
	{
		std::cerr << "Failed to Load Document: " << xmlFilename << std::endl;
	}
}

// ----------------------------------------------------------------------
// main() for printing file named on the command line
// ----------------------------------------------------------------------
int main(int argc, char* argv[])
{
	if (argc == 2)
		dump_to_stdout(argv[1]);
	else if (argc == 3)
		dump_to_file(argv[1],argv[2]);	

	return 0;
}
