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

#include <QtCore/qstring.h>
#include <QtCore/QStringList>

#include <tinyxml/tinyxml.h>
#include <ftkCommon/ftkUtils.h>
#include <ftkFeatures/ftkObjectAssociation.h>
#include <PixelAnalysis/ftkPixelLevelAnalysis.h>

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
	enum TaskType { PREPROCESSING, NUCLEAR_SEGMENTATION, FEATURE_COMPUTATION, CYTOPLASM_SEGMENTATION, RAW_ASSOCIATIONS, MULTI_MODEL_SEGMENTATION, CLASSIFY, CLASSIFY_MCLR, CLASS_EXTRACTION, ANALYTE_MEASUREMENTS, PIXEL_ANALYSIS, QUERY };
	typedef struct { int number; std::string name; std::string type; } Channel;
	typedef struct { std::string name; double value; } Parameter;
	typedef struct { std::string ClassColName; double ConfThreshold; std::string TrainingFileName; } ClassificationRule;
	typedef struct { std::string Class_Name; int Class; } ClassExtractionRule;
	typedef struct { std::string TrainingColumn; std::vector<std::string> ClassificationColumns; } ClassParam;
	typedef struct { std::string filterName; std::string channelName; std::string paramenter1; double value1; std::string paramenter2; double value2; std::string paramenter3; double value3; std::string paramenter4; double value4;
					 std::string paramenter5; double value5; std::string paramenter6; double value6; } preprocessParam;
	typedef struct { std::string name; std::string value; } QueryParameter;
	typedef struct { std::string type; std::string path; } mmSegFile;
	//FUNCTIONS:
	ProjectDefinition();
	bool Load(std::string filename);
	bool Write(std::string filename);
	void MakeDefaultNucleusSegmentation(int nucChan);
	void Clear(void);
	std::vector<Channel> ReadChannels(TiXmlElement * inputElement);
	std::vector<TaskType> ReadSteps(TiXmlElement * pipelineElement);
	std::vector<preprocessParam> ReadPreprocessingParameters(TiXmlElement * inputElement);
	std::vector<Parameter> ReadParameters(TiXmlElement * inputElement);
	std::vector<ftk::AssociationRule> ReadAssociationRules(TiXmlElement * inputElement);
	std::vector<mmSegFile> ReadMMSegFiles(TiXmlElement * inputElement);
	std::vector<std::string> ParseText(TiXmlElement * element);
	std::vector<ClassParam> ReadClassificationParameters(TiXmlElement * inputElement);
	std::vector<ClassificationRule> ReadClassificationMCLRRules(TiXmlElement * inputElement);
	std::vector<ClassExtractionRule> ReadClassExtractionRules(TiXmlElement * inputElement);
	std::vector<ftk::PixelAnalysisDefinitions> ReadPixelLevelRules(TiXmlElement * element);
	TiXmlElement * GetPreprocessingElement( preprocessParam param );
	std::vector<QueryParameter> ReadQueryParameters(TiXmlElement * inputElement);
	TiXmlElement * GetParameterElement( Parameter param );
	TiXmlElement * GetQueryParameterElement( QueryParameter param );
	TiXmlElement * GetAssocRuleElement( ftk::AssociationRule rule );
	TiXmlElement * GetClassificationElement( ClassParam ClassParameter );
	TiXmlElement * GetClassificationMCLRElement( ClassificationRule ClassRule );
	TiXmlElement * GetClassExtractRuleElement( ClassExtractionRule ClassExtractRule );
	TiXmlElement * GetTrainingFileElement( std::string file_name );
	TiXmlElement * GetPixelLevelParameterElement( ftk::PixelAnalysisDefinitions PixParameter );

	int FindInputChannel(std::string name);				//Search inputs with this name
	std::string GetTaskString(TaskType task);
	std::string GetAssocTypeString(ftk::AssociationType type);
	std::string GetBoolString(bool b);

	//VARIABLES
	std::string name;

	std::vector<Channel> inputs;
	std::vector<TaskType> pipeline;

	std::vector<preprocessParam> preprocessingParameters;
	std::vector<Parameter> nuclearParameters;
	std::vector<Parameter> cytoplasmParameters;
	std::vector<ClassParam> classificationParameters;

	std::vector<ftk::AssociationRule> associationRules;
	std::vector<mmSegFile> mmSegFiles;
	std::vector<ftk::PixelAnalysisDefinitions> pixelLevelRules;
	std::vector<std::string> intrinsicFeatures;
	std::vector<std::string> analyteMeasures;
	std::vector<ftk::ProjectDefinition::QueryParameter> queryParameters;

	std::vector<ClassificationRule> Classification_Rules;	//These rules are only for MCLR
	std::vector<ClassExtractionRule> Class_Extraction_Rules;	//This is to filter the label image 
	std::string classificationTrainingData;	//Training Data File

protected:
private:
};

}  // end namespace ftk

#endif	// end __ftkProjectDefinition_h
