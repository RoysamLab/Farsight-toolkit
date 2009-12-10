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
#include "ftkProjectDefinition.h"

namespace ftk
{

//Constructor
ProjectDefinition::ProjectDefinition()
{
	this->Clear();
}

//************************************************************************
//************************************************************************
//************************************************************************
// LOAD TOOLS
//************************************************************************
//************************************************************************
bool ProjectDefinition::Load(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ProjectDefinition" ) != 0 )
		return false;

	name = rootElement->Attribute("name");

	TiXmlElement * parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "Inputs" ) == 0 )
		{
			inputs = this->ReadChannels(parentElement);
		}
		else if( strcmp( parent, "Pipeline" ) == 0 )
		{
			pipeline = this->ReadSteps(parentElement);
		}
		else if( strcmp( parent, "NuclearSegmentationParameters" ) == 0 )
		{
			nuclearParameters = this->ReadParameters(parentElement);
		}
		else if( strcmp( parent, "CytoplasmSegmentationParameters" ) == 0 )
		{
			 cytoplasmParameters = this->ReadParameters(parentElement);
		}
		else if( strcmp( parent, "ClassificationParameters" ) == 0 )
		{
			 classificationParameters = this->ReadParameters(parentElement);
		}
		else if( strcmp( parent, "AssociationRules" ) == 0 )
		{
		}
		else if( strcmp( parent, "IntrinsicFeatures" ) == 0 )
		{
		}
		else if( strcmp( parent, "AnalyteMeasures" ) == 0 )
		{
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	return true;
}

std::vector<ProjectDefinition::Channel> ProjectDefinition::ReadChannels(TiXmlElement * inputElement)
{
	std::vector<Channel> returnVector;

	TiXmlElement * channelElement = inputElement->FirstChildElement();
	while (channelElement)
	{
		const char * parent = channelElement->Value();
		if ( strcmp( parent, "channel" ) == 0 )
		{
			Channel channel;
			channel.number = atoi(channelElement->Attribute("number"));
			channel.name = channelElement->Attribute("name");
			channel.type = channelElement->Attribute("type");
			returnVector.push_back(channel);
		}
		channelElement = channelElement->NextSiblingElement();
	} // end while(channelElement)
	return returnVector;
}

std::vector<ProjectDefinition::TaskType> ProjectDefinition::ReadSteps(TiXmlElement * pipelineElement)
{
	std::vector<TaskType> returnVector;

	TiXmlElement * stepElement = pipelineElement->FirstChildElement();
	while (stepElement)
	{
		const char * parent = stepElement->Value();
		if ( strcmp( parent, "step" ) == 0 )
		{
			std::string step = stepElement->Attribute("name");
			if( step == "NUCLEAR_SEGMENTATION")
				returnVector.push_back(NUCLEAR_SEGMENTATION);
			else if(step == "CYTOPLASM_SEGMENTATION")
				returnVector.push_back(CYTOPLASM_SEGMENTATION);
			else if(step == "RAW_ASSOCIATIONS")
				returnVector.push_back(RAW_ASSOCIATIONS);
			else if(step == "CLASSIFY")
				returnVector.push_back(CLASSIFY);
			else if(step == "ANALYTE_MEASUREMENTS")
				returnVector.push_back(ANALYTE_MEASUREMENTS);
		}
		stepElement = stepElement->NextSiblingElement();
	} // end while(stepElement)
	return returnVector;
}

std::vector<ProjectDefinition::Parameter> ProjectDefinition::ReadParameters(TiXmlElement * inputElement)
{
	std::vector<Parameter> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parent = parameterElement->Value();
		if ( strcmp( parent, "channel" ) == 0 )
		{
			Parameter parameter;
			parameter.name = parameterElement->Attribute("name");
			parameter.value = atoi(parameterElement->Attribute("value"));
			returnVector.push_back(parameter);
		}
		parameterElement = parameterElement->NextSiblingElement();
	} // end while(parentElement)
	return returnVector;
}

//************************************************************************
//************************************************************************
//************************************************************************
// WRITE TOOLS
//************************************************************************
//************************************************************************
bool ProjectDefinition::Write(std::string filename)
{
	TiXmlDocument doc;   
	TiXmlElement * root = new TiXmlElement( "ProjectDefinition" );
	root->SetAttribute("name", name.c_str());
	doc.LinkEndChild( root );

	//INPUTS:
	if(inputs.size() > 0)
	{
		TiXmlElement * inputElement = new TiXmlElement("Inputs");
		for(int i=0; i<inputs.size(); ++i)
		{
			TiXmlElement *chElement = new TiXmlElement("channel");
			chElement->SetAttribute("number", ftk::NumToString(inputs.at(i).number));
			chElement->SetAttribute("name", inputs.at(i).name);
			chElement->SetAttribute("type", inputs.at(i).type);
			inputElement->LinkEndChild(chElement);
		}
		root->LinkEndChild(inputElement);
	}

	//PIPELINE:
	if(pipeline.size() > 0)
	{
		TiXmlElement * pipelineElement = new TiXmlElement("Pipeline");
		for(int i=0; i<pipeline.size(); ++i)
		{
			TiXmlElement * stepElement = new TiXmlElement("step");
			stepElement->SetAttribute("name", GetTaskString(pipeline.at(i)));
			pipelineElement->LinkEndChild(stepElement);
		}
		root->LinkEndChild(pipelineElement);
	}

	//NuclearSegmentationParameters:
	if(nuclearParameters.size() > 0)
	{
		TiXmlElement * paramsElement = new TiXmlElement("NuclearSegmentationParameters");
		for(int i=0; i<nuclearParameters.size(); ++i)
		{
			paramsElement->LinkEndChild( GetParameterElement(nuclearParameters.at(i)) );
		}
		root->LinkEndChild(paramsElement);
	}

	//CytoplasmSegmentationParameters:
	if(cytoplasmParameters.size() > 0)
	{
		TiXmlElement * paramsElement = new TiXmlElement("CytoplasmSegmentationParameters");
		for(int i=0; i<cytoplasmParameters.size(); ++i)
		{
			paramsElement->LinkEndChild( GetParameterElement(cytoplasmParameters.at(i)) );
		}
		root->LinkEndChild(paramsElement);
	}

	//ClassificationParameters:
	if(classificationParameters.size() > 0)
	{
		TiXmlElement * paramsElement = new TiXmlElement("ClassificationParameters");
		for(int i=0; i<classificationParameters.size(); ++i)
		{
			paramsElement->LinkEndChild( GetParameterElement(classificationParameters.at(i)) );
		}
		root->LinkEndChild(paramsElement);
	}

	//AssociationRules:
	//IntrinsicFeatures:
	///AnalyteMeasures:

	if(doc.SaveFile( filename.c_str() ))
		return true;
	else
		return false;
}

TiXmlElement * ProjectDefinition::GetParameterElement( Parameter param )
{
	TiXmlElement * returnElement = new TiXmlElement("parameter");
	returnElement->SetAttribute("name", param.name);
	returnElement->SetAttribute("value", ftk::NumToString(param.value));
	return returnElement;
}
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
void ProjectDefinition::MakeDefaultNucleusSegmentation(int nucChan)
{
	this->Clear();

	this->name = "default_nuclei";
	
	Channel c;
	c.name = "nuc";
	c.type = "NUCLEAR";
	c.number = nucChan;

	this->inputs.push_back(c);
	this->pipeline.push_back(NUCLEAR_SEGMENTATION);

	//Also need parameters:
}

//************************************************************************
//************************************************************************
//************************************************************************
//Search inputs for channel with this name:
int ProjectDefinition::FindInputChannel(std::string name)
{
	int retval = -1;

	//First look for the Nuclear Segmentation
	for(int i=0; i<inputs.size(); ++i)
	{
		if(inputs.at(i).type == name)
		{
			retval = inputs.at(i).number;
			break;
		}
	}
	return retval;
}

std::string ProjectDefinition::GetTaskString(TaskType task)
{
	std::string retText = "";
	switch(task)
	{
	case NUCLEAR_SEGMENTATION:
		retText = "NUCLEAR_SEGMENTATION";
		break;
	case CYTOPLASM_SEGMENTATION:
		retText = "CYTOPLASM_SEGMENTATION";
		break;
	case RAW_ASSOCIATIONS:
		retText = "RAW_ASSOCIATIONS";
		break;
	case CLASSIFY:
		retText = "CLASSIFY";
		break;
	case ANALYTE_MEASUREMENTS:
		retText = "ANALYTE_MEASUREMENTS";
		break;
	}
	return retText;
}

//************************************************************************
//************************************************************************
//************************************************************************
void ProjectDefinition::Clear(void)
{
	name.clear();
	inputs.clear();
	pipeline.clear();
	nuclearParameters.clear();
	cytoplasmParameters.clear();
	classificationParameters.clear();
	associationRules.clear();
	intrinsicFeatures.clear();
	analyteMeasures.clear();
	classificationTrainingData.clear();
}

}  // end namespace ftk
