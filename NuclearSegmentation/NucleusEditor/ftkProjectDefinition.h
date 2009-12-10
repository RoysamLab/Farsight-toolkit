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
#ifndef __ftkProjectDefinition_h
#define __ftkProjectDefinition_h

#include <tinyxml/tinyxml.h>
#include <ftkCommon/ftkUtils.h>

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace ftk
{

class ProjectDefinition
{
public:
	//ENUMS & STRUCTS:
	enum TaskType { NUCLEAR_SEGMENTATION, CYTOPLASM_SEGMENTATION, RAW_ASSOCIATIONS, CLASSIFY, ANALYTE_MEASUREMENTS };
	typedef struct { int number; std::string name; std::string type; } Channel;
	typedef struct { std::string name; double value; } Parameter;

	//FUNCTIONS:
	ProjectDefinition();
	bool Load(std::string filename);
	bool Write(std::string filename);
	void MakeDefaultNucleusSegmentation(int nucChan);
	void Clear(void);
	std::vector<Channel> ReadChannels(TiXmlElement * inputElement);
	std::vector<TaskType> ReadSteps(TiXmlElement * pipelineElement);
	std::vector<Parameter> ReadParameters(TiXmlElement * inputElement);
	TiXmlElement * GetParameterElement( Parameter param );

	int FindInputChannel(std::string name);				//Search inputs with this name
	std::string GetTaskString(TaskType task);

	//VARIABLES
	std::string name;

	std::vector<Channel> inputs;
	std::vector<TaskType> pipeline;

	std::vector<Parameter> nuclearParameters;
	std::vector<Parameter> cytoplasmParameters;
	std::vector<Parameter> classificationParameters;

	std::vector<std::string> associationRules;
	std::vector<std::string> intrinsicFeatures;
	std::vector<std::string> analyteMeasures;

	std::vector<std::string> classificationTrainingData;

protected:
private:
};

}  // end namespace ftk

#endif	// end __ftkProjectDefinition_h
