/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

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
	this->ClearAll();
}

void ProjectFiles::ClearAll()
{
	path = "";
	name = "";
	type = "";

	//Filenames:
	input = "";
	output = "";
	log = "";
	definition = "";
	table = "";

	inputSaved = true;
	outputSaved = true;
	//logSaved = true;
	definitionSaved = true;
	tableSaved = true;
	adjTablesSaved = true;
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

	path = rootElement->Attribute("path");
	name = rootElement->Attribute("name");
	if(rootElement->Attribute("type"))		
		type = rootElement->Attribute("type"); 

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
			 nucSegValidated = parentElement->Attribute("validated");
		}
		else if( strcmp( parent, "definition" ) == 0 )
		{
			 definition = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "table" ) == 0 )
		{
			table = parentElement->Attribute("file");
		}
		else if( strcmp( parent, "adjTables" ) == 0 )
		{
			adjTables = parentElement->Attribute("file");
		}

		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	return true;
}

bool ProjectFiles::Write()
{
	TiXmlDocument doc;   
	TiXmlElement * root = new TiXmlElement( "ProjectFiles" );
	root->SetAttribute("path", path.c_str());
	root->SetAttribute("name", name.c_str());
	if(strcmp(  type.c_str(), "multi" ) == 0 )
		root->SetAttribute("type", type.c_str());

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
	file->SetAttribute("validated", ftk::NumToString(nucSegValidated));
	root->LinkEndChild(file);

	file = new TiXmlElement("definition");
	file->SetAttribute("file", definition.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("table");
	file->SetAttribute("file", table.c_str());
	root->LinkEndChild(file);

	file = new TiXmlElement("adjTables");
	file->SetAttribute("file", adjTables.c_str());
	root->LinkEndChild(file);

	if(doc.SaveFile( this->GetFullName().c_str() ))
		return true;
	else
		return false;
}

std::string ProjectFiles::GetFullName()
{
	std::string full = path + name + ".xml";
	return full;
}

std::string ProjectFiles::GetFullInput()
{
	std::string full = path + input;
	return full;
}

std::string ProjectFiles::GetFullOutput()
{
	std::string full = path + output;
	return full;
}

std::string ProjectFiles::GetFullLog()
{
	std::string full = path + log;
	return full;
}

std::string ProjectFiles::GetFullDef()
{
	std::string full = path + definition;
	return full;
}

std::string ProjectFiles::GetFullTable()
{
	std::string full = path + table;
	return full;
}

std::string ProjectFiles::GetFullAdjTables()
{
	std::string full = path + adjTables;
	return full;
}
//************************************************************************
//************************************************************************
//************************************************************************


}  // end namespace ftk
