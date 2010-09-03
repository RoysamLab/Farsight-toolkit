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
		else if( strcmp( parent, "PreprocessingParameters" ) == 0 )
		{
			preprocessingParameters = this->ReadPreprocessingParameters(parentElement);
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
			classificationParameters = this->ReadClassificationParameters(parentElement);
		}
		else if( strcmp( parent, "AssociationRules" ) == 0 )
		{
			associationRules = this->ReadAssociationRules(parentElement);
		}
		else if( strcmp( parent, "IntrinsicFeatures" ) == 0 )
		{
			intrinsicFeatures = this->ParseText(parentElement);
		}
		else if( strcmp( parent, "AnalyteMeasures" ) == 0 )
		{
		}
		else if( strcmp( parent, "PixelLevelAnalysis" ) == 0 )
		{
			pixelLevelRules = this->ReadPixelLevelRules(parentElement);
		}
		else if( strcmp( parent, "SqlQueryParameters" ) == 0 )
		{
			queryParameters = this->ReadQueryParameters(parentElement);
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
			if( step == "PREPROCESSING")
				returnVector.push_back(PREPROCESSING);
			else if( step == "NUCLEAR_SEGMENTATION")
				returnVector.push_back(NUCLEAR_SEGMENTATION);
			else if(step == "CYTOPLASM_SEGMENTATION")
				returnVector.push_back(CYTOPLASM_SEGMENTATION);
			else if(step == "RAW_ASSOCIATIONS")
				returnVector.push_back(RAW_ASSOCIATIONS);
			else if(step == "CLASSIFY")
				returnVector.push_back(CLASSIFY);
			else if(step == "ANALYTE_MEASUREMENTS")
				returnVector.push_back(ANALYTE_MEASUREMENTS);
			else if(step == "PIXEL_ANALYSIS")
				returnVector.push_back(PIXEL_ANALYSIS);
			else if(step == "QUERY")
				returnVector.push_back(QUERY);
		}
		stepElement = stepElement->NextSiblingElement();
	} // end while(stepElement)
	return returnVector;
}

std::vector<ProjectDefinition::preprocessParam > ProjectDefinition::ReadPreprocessingParameters(TiXmlElement * inputElement){

	std::vector<ProjectDefinition::preprocessParam> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parameter = parameterElement ->Value();
		ProjectDefinition::preprocessParam prepStep;
		if ( strcmp(parameter,"PreprocessingStep") == 0 )
		{
			prepStep.filterName = parameterElement->Attribute("Name");
			prepStep.channelName = parameterElement->Attribute("Channel");
			if( parameterElement->Attribute("Parameter1") != NULL ){
				prepStep.paramenter1 = parameterElement->Attribute("Parameter1");
				prepStep.value1 = atof(parameterElement->Attribute("Value1"));}
			if( parameterElement->Attribute("Parameter2") != NULL ){
				prepStep.paramenter2 = parameterElement->Attribute("Parameter2");
				prepStep.value2 = atof(parameterElement->Attribute("Value2"));}
			if( parameterElement->Attribute("Parameter3") != NULL ){
				prepStep.paramenter3 = parameterElement->Attribute("Parameter3");
				prepStep.value3 = atof(parameterElement->Attribute("Value3"));}
			if( parameterElement->Attribute("Parameter4") != NULL ){
				prepStep.paramenter4 = parameterElement->Attribute("Parameter4");
				prepStep.value4 = atof(parameterElement->Attribute("Value4"));}
			if( parameterElement->Attribute("Parameter5") != NULL ){
				prepStep.paramenter5 = parameterElement->Attribute("Parameter5");
				prepStep.value5 = atof(parameterElement->Attribute("Value5"));}
			if( parameterElement->Attribute("Parameter6") != NULL ){
				prepStep.paramenter6 = parameterElement->Attribute("Parameter6");
				prepStep.value6 = atof(parameterElement->Attribute("Value6"));}
		}
		returnVector.push_back(prepStep);
		parameterElement = parameterElement->NextSiblingElement();
	}
	return returnVector;
}

std::vector<ProjectDefinition::Parameter> ProjectDefinition::ReadParameters(TiXmlElement * inputElement)
{
	std::vector<Parameter> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parent = parameterElement->Value();
		if ( strcmp( parent, "parameter" ) == 0 )
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

std::vector<ftk::ProjectDefinition::QueryParameter> ProjectDefinition::ReadQueryParameters(TiXmlElement * inputElement)
{
	std::vector<ftk::ProjectDefinition::QueryParameter> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parent = parameterElement->Value();
		if ( strcmp( parent, "SqlQuery" ) == 0 )
		{
			ftk::ProjectDefinition::QueryParameter parameter;
			parameter.name = parameterElement->Attribute("name");
			parameter.value = (parameterElement->Attribute("value"));
			returnVector.push_back(parameter);
		}
		parameterElement = parameterElement->NextSiblingElement();
	} // end while(parentElement)
	return returnVector;
}

std::vector<ftk::AssociationRule> ProjectDefinition::ReadAssociationRules(TiXmlElement * inputElement)
{
	std::vector<ftk::AssociationRule> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parameter = parameterElement ->Value();
		ftk::AssociationRule assocRule("");
		if ( strcmp(parameter,"AssociationRule") == 0 )
		{
			assocRule.SetRuleName(parameterElement->Attribute("Name"));
			assocRule.SetSegmentationFileNmae(parameterElement->Attribute("SegmentationSource"));
			assocRule.SetTargetFileNmae(parameterElement->Attribute("Target_Image"));
			assocRule.SetOutDistance(atoi(parameterElement->Attribute("Outside_Distance")));
			assocRule.SetInDistance(atoi(parameterElement->Attribute("Inside_Distance")));

			if(strcmp(parameterElement->Attribute("Use_Whole_Object"),"True")==0)
				assocRule.SetUseWholeObject(true);
			else
				assocRule.SetUseWholeObject(false);

			if(strcmp(parameterElement->Attribute("Use_Background_Subtraction"),"True")==0)
				assocRule.SetUseBackgroundSubtraction(true);
			else
				assocRule.SetUseBackgroundSubtraction(false);

			if(strcmp(parameterElement->Attribute("Use_MultiLevel_Thresholding"),"True")==0){
				assocRule.SetUseMultiLevelThresholding(true);
				assocRule.SetNumberOfThresholds(atoi(parameterElement->Attribute("Number_Of_Thresholds")));
				assocRule.SetNumberIncludedInForeground(atoi(parameterElement->Attribute("Number_Included_In_Foreground")));
			}
			else{
				assocRule.SetUseMultiLevelThresholding(false);
				assocRule.SetNumberOfThresholds(1);
				assocRule.SetNumberIncludedInForeground(1);
			}

			if(strcmp(parameterElement->Attribute("Association_Type"),"MIN")==0)
				assocRule.SetAssocType(ASSOC_MIN);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"MAX")==0)
				assocRule.SetAssocType(ASSOC_MAX);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"TOTAL")==0)
				assocRule.SetAssocType(ASSOC_TOTAL);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"SURROUNDEDNESS")==0)
				assocRule.SetAssocType(ASSOC_SURROUNDEDNESS);
			else
				assocRule.SetAssocType(ASSOC_AVERAGE);
		}
		returnVector.push_back(assocRule);
		parameterElement = parameterElement->NextSiblingElement();
	}
	return returnVector;
}

std::vector<ProjectDefinition::ClassParam> ProjectDefinition::ReadClassificationParameters(TiXmlElement * inputElement){
	std::vector<ProjectDefinition::ClassParam> returnVector;
	TiXmlElement * parameterElement = inputElement->FirstChildElement();

	while (parameterElement){
		const char * parameter = parameterElement ->Value();
		if ( strcmp(parameter,"ClassificationParameter") == 0 ){
			ClassParam current_classification_parameter;
			current_classification_parameter.TrainingColumn = parameterElement->Attribute("TrainingColumn");
			std::string ClassColNames;
			ClassColNames = parameterElement->Attribute("ClassificationColumns");
			//Split by ","
			QString QClassColNames;
			QClassColNames = QString::fromStdString(ClassColNames);
			QStringList values = QClassColNames.split(",");
			for(int i=0; i<values.size(); ++i)
				current_classification_parameter.ClassificationColumns.push_back( values.at(i).toStdString() );

			returnVector.push_back(current_classification_parameter);
		}
		else if ( strcmp(parameter,"TrainingFile") == 0 ){
			classificationTrainingData = parameterElement->Attribute("Name");
		}
		parameterElement = parameterElement->NextSiblingElement();
	}
	return returnVector;
}


std::vector<ftk::PixelAnalysisDefinitions> ProjectDefinition::ReadPixelLevelRules(TiXmlElement * inputElement)
{

	std::vector<ftk::PixelAnalysisDefinitions> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parameter = parameterElement ->Value();
		ftk::PixelAnalysisDefinitions pixRule;
		if ( strcmp(parameter,"PixelLevelRule") == 0 ){
			pixRule.regionChannelName = parameterElement->Attribute("RoiImage");
			pixRule.targetChannelName = parameterElement->Attribute("TargetImage");
			pixRule.mode = atoi(parameterElement->Attribute("Mode"));
			pixRule.outputFilename = parameterElement->Attribute("OutputFilename");
			pixRule.radius = atoi(parameterElement->Attribute("Radius"));
		}
		returnVector.push_back(pixRule);
		parameterElement = parameterElement->NextSiblingElement();
	}
	return returnVector;
}



std::vector<std::string> ProjectDefinition::ParseText(TiXmlElement * element)
{
	std::vector<std::string> returnVector;

	std::string text = element->GetText();

	size_t begin = 0;
	size_t end = text.find_first_of(",");
	while (end != std::string::npos)
	{
		returnVector.push_back(text.substr(begin,end-begin));
		begin = end+1;
		end = text.find_first_of(",", begin);
	}
	returnVector.push_back(text.substr(begin,end-begin));
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
		for(int i=0; i<(int)inputs.size(); ++i)
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
		for(int i=0; i<(int)pipeline.size(); ++i)
		{
			TiXmlElement * stepElement = new TiXmlElement("step");
			stepElement->SetAttribute("name", GetTaskString(pipeline.at(i)));
			pipelineElement->LinkEndChild(stepElement);
		}
		root->LinkEndChild(pipelineElement);
	}

	//PreprocessParameters
	if(preprocessingParameters.size() > 0){
		TiXmlElement * paramElement = new TiXmlElement("PreprocessingParameters");
		for(int i=0; i<(int)preprocessingParameters.size(); ++i){
			paramElement->LinkEndChild( GetPreprocessingElement(preprocessingParameters.at(i)) );
		}
		root->LinkEndChild(paramElement);
	}

	//NuclearSegmentationParameters:
	if(nuclearParameters.size() > 0)
	{
		TiXmlElement * paramsElement = new TiXmlElement("NuclearSegmentationParameters");
		for(int i=0; i<(int)nuclearParameters.size(); ++i)
		{
			paramsElement->LinkEndChild( GetParameterElement(nuclearParameters.at(i)) );
		}
		root->LinkEndChild(paramsElement);
	}

	//CytoplasmSegmentationParameters:
	if(cytoplasmParameters.size() > 0){
		TiXmlElement * paramsElement = new TiXmlElement("CytoplasmSegmentationParameters");
		for(int i=0; i<(int)cytoplasmParameters.size(); ++i){
			paramsElement->LinkEndChild( GetParameterElement(cytoplasmParameters.at(i)) );
		}
		root->LinkEndChild(paramsElement);
	}

	//IntrinsicFeatures:
	if(intrinsicFeatures.size() > 0){
		TiXmlElement * featuresElement = new TiXmlElement("IntrinsicFeatures");
		std::string text = intrinsicFeatures.at(0);
		for(int i=1; i<(int)intrinsicFeatures.size(); ++i){
			text += "," + intrinsicFeatures.at(i);
		}
		featuresElement->LinkEndChild( new TiXmlText( text.c_str() ) );
		root->LinkEndChild(featuresElement);
	}

	//AssociationRules:
	if(associationRules.size() > 0){
		TiXmlElement * assocElement = new TiXmlElement("AssociationRules");
		for(int i=0; i<(int)associationRules.size(); ++i){
			assocElement->LinkEndChild( GetAssocRuleElement(associationRules.at(i)) );
		}
		root->LinkEndChild(assocElement);
	}

	//ClassificationParameters:
	if(classificationParameters.size() > 0){
		TiXmlElement * paramsElement = new TiXmlElement("ClassificationParameters");
		if( !classificationTrainingData.empty() )
			paramsElement->LinkEndChild(GetTrainingFileElement(classificationTrainingData));
		for(int i=0; i<(int)classificationParameters.size(); ++i)
			paramsElement->LinkEndChild(GetClassificationElement(classificationParameters.at(i)));
		root->LinkEndChild(paramsElement);
	}

	//ADD PIXEL DEFINITIONS & SQLQUERY

	//AnalyteMeasures:

	if(doc.SaveFile( filename.c_str() ))
		return true;
	else
		return false;
}

TiXmlElement * ProjectDefinition::GetPreprocessingElement( preprocessParam param )
{
	TiXmlElement * returnElement = new TiXmlElement("PreprocessingStep");
	returnElement->SetAttribute("Name", param.filterName );
	returnElement->SetAttribute("Channel", param.channelName );
	if( !param.paramenter1.empty() ){
		returnElement->SetAttribute("Parameter1", param.paramenter1 );
		returnElement->SetAttribute("Value1", ftk::NumToString(param.value1));
	}
	if( !param.paramenter2.empty() ){
		returnElement->SetAttribute("Parameter2", param.paramenter2 );
		returnElement->SetAttribute("Value2", ftk::NumToString(param.value2));
	}
	if( !param.paramenter3.empty() ){
		returnElement->SetAttribute("Parameter3", param.paramenter3 );
		returnElement->SetAttribute("Value3", ftk::NumToString(param.value3));
	}
	if( !param.paramenter4.empty() ){
		returnElement->SetAttribute("Parameter4", param.paramenter4 );
		returnElement->SetAttribute("Value4", ftk::NumToString(param.value4));
	}
	if( !param.paramenter5.empty() ){
		returnElement->SetAttribute("Parameter5", param.paramenter5 );
		returnElement->SetAttribute("Value5", ftk::NumToString(param.value5));
	}
	if( !param.paramenter6.empty() ){
		returnElement->SetAttribute("Parameter6", param.paramenter6 );
		returnElement->SetAttribute("Value6", ftk::NumToString(param.value6));
	}
	return returnElement;
}

TiXmlElement * ProjectDefinition::GetParameterElement( Parameter param )
{
	TiXmlElement * returnElement = new TiXmlElement("parameter");
	returnElement->SetAttribute("name", param.name);
	returnElement->SetAttribute("value", ftk::NumToString(param.value));
	return returnElement;
}

TiXmlElement * ProjectDefinition::GetAssocRuleElement( ftk::AssociationRule rule )
{
	TiXmlElement * returnElement = new TiXmlElement("AssociationRule");
	returnElement->SetAttribute("Name", rule.GetRuleName().c_str());
	returnElement->SetAttribute("SegmentationSource", rule.GetSegmentationFileName().c_str());
	returnElement->SetAttribute("Target_Image", rule.GetTargetFileNmae().c_str());
	returnElement->SetAttribute("Outside_Distance", ftk::NumToString(rule.GetOutDistance()).c_str());
	returnElement->SetAttribute("Inside_Distance", ftk::NumToString(rule.GetInDistance()).c_str());
	returnElement->SetAttribute("Use_Whole_Object", GetBoolString(rule.IsUseWholeObject()).c_str());
	returnElement->SetAttribute("Use_Background_Subtraction", GetBoolString(rule.IsUseBackgroundSubtraction()).c_str());
	returnElement->SetAttribute("Use_MultiLevel_Thresholding", GetBoolString(rule.IsUseMultiLevelThresholding()).c_str());
	returnElement->SetAttribute("Number_Of_Thresholds", ftk::NumToString(rule.GetNumberOfThresholds()).c_str());
	returnElement->SetAttribute("Number_Included_In_Foreground", ftk::NumToString(rule.GetNumberIncludedInForeground()).c_str());
	returnElement->SetAttribute("Association_Type", GetAssocTypeString(rule.GetAssocType()).c_str());
	return returnElement;
}

TiXmlElement * ProjectDefinition::GetTrainingFileElement( std::string file_name ){
	TiXmlElement * returnElement = new TiXmlElement("TrainingFile");
	returnElement->SetAttribute("Name", file_name.c_str());
	return returnElement;
}
TiXmlElement * ProjectDefinition::GetClassificationElement( ProjectDefinition::ClassParam ClassParameter ){
	TiXmlElement * returnElement = new TiXmlElement("ClassificationParameter");
	returnElement->SetAttribute("TrainingColumn", ClassParameter.TrainingColumn.c_str());
	std::string class_columns;
	class_columns = ClassParameter.ClassificationColumns.at(0);
	for(int i=1; i<(int)ClassParameter.ClassificationColumns.size(); ++i)
		class_columns = class_columns+","+ClassParameter.ClassificationColumns.at(i);
	returnElement->SetAttribute("ClassificationColumns", class_columns.c_str());
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
	for(int i=0; i<(int)inputs.size(); ++i)
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
	case PIXEL_ANALYSIS:
		retText = "PIXEL_ANALYSIS";
		break;
	}
	return retText;
}

std::string ProjectDefinition::GetAssocTypeString(ftk::AssociationType type)
{
	std::string retText = "";
	switch(type)
	{
	case ftk::ASSOC_AVERAGE:
		retText = "AVERAGE";
		break;
	case ftk::ASSOC_MAX:
		retText = "MAX";
		break;
	case ftk::ASSOC_MIN:
		retText = "MIN";
		break;
	case ftk::ASSOC_SURROUNDEDNESS:
		retText = "SURROUNDEDNESS";
		break;
	case ftk::ASSOC_TOTAL:
		retText = "TOTAL";
		break;
	}
	return retText;
}

std::string ProjectDefinition::GetBoolString(bool b)
{
	if(b)
		return "True";
	else
		return "False";
}

//************************************************************************
//************************************************************************
//************************************************************************
void ProjectDefinition::Clear(void)
{
	name.clear();
	inputs.clear();
	pipeline.clear();
	preprocessingParameters.clear();
	nuclearParameters.clear();
	cytoplasmParameters.clear();
	classificationParameters.clear();
	associationRules.clear();
	intrinsicFeatures.clear();
	analyteMeasures.clear();
	classificationTrainingData.clear();
}

}  // end namespace ftk
