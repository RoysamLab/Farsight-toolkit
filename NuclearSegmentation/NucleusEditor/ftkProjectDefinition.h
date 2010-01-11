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
#ifndef __ftkProjectDefinition_h
#define __ftkProjectDefinition_h

#include <tinyxml/tinyxml.h>
#include <ftkCommon/ftkUtils.h>
#include <ftkCommon/ftkObjectAssociation.h>

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
	void ReadAssociationParameters(TiXmlElement * inputElement);
	std::vector<std::string> ParseText(TiXmlElement * element);
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

	std::vector<ftk::AssociationRule> associationRules;
	std::vector<std::string> intrinsicFeatures;
	std::vector<std::string> analyteMeasures;

	std::vector<std::string> classificationTrainingData;

protected:
private:
};

}  // end namespace ftk

#endif	// end __ftkProjectDefinition_h
