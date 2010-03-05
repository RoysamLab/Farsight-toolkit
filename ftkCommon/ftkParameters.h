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
#ifndef __ftkParameters_h
#define __ftkParameters_h

#include <tinyxml/tinyxml.h>
#include <ftkUtils.h>
#include <string>
#include <vector>

namespace ftk
{

class Parameters
{
	enum Type { STRING, BOOL, INT, DOUBLE };
	typedef struct { std::string name; Type type; std::string value; } Parameter;

public:
	Parameters(){ extension = "prm"; };
	bool LoadFromFile(std::string filename);
	void ReadFromTinyXML(TiXmlElement * inputElement);
	bool WriteToFile(std::string filename);
	void AddToTinyXML(TiXmlElement * rootElement);

	void AddParameter(std::string name, Type type, std::string value);

	int GetNumberOfParameters(){ return (int)m_Parameters.size(); };
	int QueryParameter(std::string parameterName);

	Type GetType(int idx);
	std::string GetValue(int idx);
	int GetValueAsInt(int idx);
	double GetValueAsDouble(int idx);
	bool GetValueAsBool(int idx);

	std::string GetParent(){ return parent; };
	void SetParent(std::string prnt){ parent = prnt; };

private:
	std::vector<Parameter> m_Parameters;
	std::string parent;			//The value of the parent element
	std::string docname;
	
	std::string extension;
};

}  // end namespace ftk

#endif	// end __ftkParameters_h
