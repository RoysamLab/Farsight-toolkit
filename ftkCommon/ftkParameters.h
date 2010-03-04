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
	void LoadFromFile(std::string filename){};
	void WriteToFile(std::string filename){};

	void AddParameter(std::string name, Type type, std::string value);

private:
	std::vector<Parameter> m_Parameters;
	
	std::string extension;
};

}  // end namespace ftk

#endif	// end __ftkParameters_h
