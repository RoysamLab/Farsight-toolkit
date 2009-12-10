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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include "ftkProjectFiles.h"

namespace ftk
{

ProjectFiles::ProjectFiles()
{
	name = "";

	//Filenames:
	input = "";
	output = "";
	log = "";
	definition = "";
	table = "";

	inputSaved = true;
	outputSaved = true;
	logSaved = true;
	definitionSaved = true;
	tableSaved = true;
}

//************************************************************************
//************************************************************************
//************************************************************************
// READ TOOLS
//************************************************************************
//************************************************************************
bool ProjectFiles::Read(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ProjectFiles" ) != 0 )
		return false;

	name = rootElement->Attribute("name");

	TiXmlElement * parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "input" ) == 0 )
		{
			input = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "output" ) == 0 )
		{
			output = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "log" ) == 0 )
		{
			 log = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "definition" ) == 0 )
		{
			 definition = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "table" ) == 0 )
		{
			table = parentElement->Attribute("file");
		}

		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	return true;
}

bool ProjectFiles::Write(std::string filename)
{
	TiXmlDocument doc;   
	TiXmlElement * root = new TiXmlElement( "ProjectFiles" );
	root->SetAttribute("name", name.c_str());
	doc.LinkEndChild( root );  
 
	TiXmlElement * file;

	file = new TiXmlElement("input");
	file->SetAttribute("file", input.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("output");
	file->SetAttribute("file", output.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("log");
	file->SetAttribute("file", log.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("definition");
	file->SetAttribute("file", definition.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("table");
	file->SetAttribute("file", table.c_str());
	root->LinkEndChild(file);

	if(doc.SaveFile( filename.c_str() ))
		return true;
	else
		return false;
}
//************************************************************************
//************************************************************************
//************************************************************************


}  // end namespace ftk
