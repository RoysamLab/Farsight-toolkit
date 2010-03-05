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

#include "ftkParameters.h"

namespace ftk
{

void Parameters::AddParameter(std::string name, Type type, std::string value)
{
	Parameter p;
	p.name = name;
	p.type = type;
	p.value = value;
	m_Parameters.push_back(p);
}


int Parameters::QueryParameter(std::string parameterName)
{
	int retVal = -1;

	for(int i=0; i<(int)m_Parameters.size(); ++i)
	{
		if(m_Parameters.at(i).name == parameterName)
		{
			retVal = i;
			break;
		}
	}
	return retVal;
}

Parameters::Type Parameters::GetType(int idx)
{
	if(idx < (int)m_Parameters.size())
		return m_Parameters.at(idx).type;
	else
		return STRING;
}
std::string Parameters::GetValue(int idx)
{
	if(idx < (int)m_Parameters.size())
		return m_Parameters.at(idx).value;
	else
		return std::string("");
}
int Parameters::GetValueAsInt(int idx)
{
	return atoi( GetValue(idx).c_str() );
}
double Parameters::GetValueAsDouble(int idx)
{
	return atof( GetValue(idx).c_str() );
}
bool Parameters::GetValueAsBool(int idx)
{
	std::string val = GetValue(idx);
	if( val == "true" || val == "TRUE")
		return true;
	else
		return false;
}

bool Parameters::LoadFromFile(std::string filename)
{
	if( GetExtension(filename) != this->extension )
		return false;

	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	docname = rootElement->Value();

	this->ReadFromTinyXML(rootElement);
	return true;
}

void Parameters::ReadFromTinyXML(TiXmlElement * inputElement)
{
	parent = inputElement->Value();

	m_Parameters.clear();	//Clear parameters list

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * ch = parameterElement->Value();
		if ( strcmp( ch, "parameter" ) == 0 )
		{
			const char * name = parameterElement->Attribute("name");
			const char * type = parameterElement->Attribute("type");
			const char * valu = parameterElement->Attribute("value");

			Parameter p;
			if(name != NULL)
				p.name = name;
			if(type != NULL)
			{
				if(strcmp(type, "STRING") == 0)
					p.type = STRING;
				else if(strcmp(type, "BOOL") == 0)
					p.type = BOOL;
				else if(strcmp(type, "INT") == 0)
					p.type = INT;
				else if(strcmp(type, "DOUBLE") == 0)
					p.type = DOUBLE;
			}
			if(valu != NULL)
				p.value = valu;

			m_Parameters.push_back(p);
		}
		parameterElement = parameterElement->NextSiblingElement();
	} // end while(parentElement)
}

bool Parameters::WriteToFile(std::string filename)
{
	if( GetExtension(filename) != this->extension )
		return false;

	TiXmlDocument doc;   
	TiXmlElement * root = new TiXmlElement( parent.c_str() );
	
	doc.LinkEndChild( root );

	this->AddToTinyXML(root);

	return doc.SaveFile( filename.c_str() );

}

void Parameters::AddToTinyXML(TiXmlElement * rootElement)
{
	if((int)m_Parameters.size() > 0)
	{
		for(int i=0; i<(int)m_Parameters.size(); ++i)
		{
			Parameter p = m_Parameters.at(i);
			TiXmlElement * pElement = new TiXmlElement("parameter");
			pElement->SetAttribute("name", p.name);
			switch(p.type)
			{
			case STRING:
				pElement->SetAttribute("type", "STRING");
				break;
			case INT:
				pElement->SetAttribute("type", "INT");
				break;
			case DOUBLE:
				pElement->SetAttribute("type", "DOUBLE");
				break;
			case BOOL:
				pElement->SetAttribute("type", "BOOL");
				break;
			}
			pElement->SetAttribute("value", p.value);
			rootElement->LinkEndChild( pElement );
		}
	}
}

}  // end namespace ftk
