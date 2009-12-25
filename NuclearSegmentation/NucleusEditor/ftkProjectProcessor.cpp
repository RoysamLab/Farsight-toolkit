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
#include "ftkProjectProcessor.h"

namespace ftk
{

//Constructor
ProjectProcessor::ProjectProcessor()
{
	inputImage = NULL;
	outputImage = NULL;
	definition = NULL;
	tasks.clear();
	table = NULL;
	numTasks = 0;
	lastTask = -1;
	resultIsEditable = false;
	inputTypeNeeded = 0;
}

void ProjectProcessor::Initialize(void)
{
	if(!definition) return;

	tasks.clear();

	for(int i=0; i<(int)definition->pipeline.size(); ++i)
	{
		Task t;
		t.type = definition->pipeline.at(i);
		switch(t.type)
		{
		case ProjectDefinition::NUCLEAR_SEGMENTATION:
			t.inputChannel1 = definition->FindInputChannel("NUCLEAR");
			break;
		case ProjectDefinition::CYTOPLASM_SEGMENTATION:
			t.inputChannel2 = definition->FindInputChannel("CYTOPLASM");
			t.inputChannel3 = definition->FindInputChannel("MEMBRANE");
			break;
		case ProjectDefinition::RAW_ASSOCIATIONS:
			break;
		case ProjectDefinition::CLASSIFY:
			break;
		case ProjectDefinition::ANALYTE_MEASUREMENTS:
			break;
		}
		t.done = false;
		tasks.push_back(t);
	}

	numTasks = (int)tasks.size();
	lastTask = -1;
}

void ProjectProcessor::ProcessNext(void)
{
	if(DoneProcessing()) return;

	int thisTask = lastTask + 1;
	bool taskDone = false;
	switch(tasks.at(thisTask).type)
	{
	case ProjectDefinition::NUCLEAR_SEGMENTATION:
		taskDone = SegmentNuclei(tasks.at(thisTask).inputChannel1);
		break;
	case ProjectDefinition::CYTOPLASM_SEGMENTATION:
		taskDone = SegmentCytoplasm(tasks.at(thisTask).inputChannel2, tasks.at(thisTask).inputChannel3);
		break;
	case ProjectDefinition::RAW_ASSOCIATIONS:
		taskDone = ComputeAssociations();
		break;
	case ProjectDefinition::ANALYTE_MEASUREMENTS:
		taskDone = false;
		break;
	case ProjectDefinition::CLASSIFY:
		taskDone = false;
		break;
	}
	
	if(taskDone)
	{
		tasks.at(thisTask).done = true;
		lastTask++;
	}
}

bool ProjectProcessor::SegmentNuclei(int nucChannel)
{
	if(!inputImage)
		return false;

	ftk::NuclearSegmentation * nucSeg = new ftk::NuclearSegmentation();
	nucSeg->SetInput(inputImage, "data_image", nucChannel);

	//Setup the Parameters:
	bool finalize = false;	//The Default
	for(int i=0; i<(int)definition->nuclearParameters.size(); ++i)
	{
		nucSeg->SetParameter(definition->nuclearParameters.at(i).name, int(definition->nuclearParameters.at(i).value));
		if(definition->nuclearParameters.at(i).name == "finalize_segmentation")
		{
			if(definition->nuclearParameters.at(i).value == 1)
				finalize = true;
		}
	}

	//Process:
	nucSeg->Binarize(false);
	nucSeg->DetectSeeds(false);
	if(finalize)
	{
		nucSeg->RunClustering(false);
		nucSeg->Finalize();
	}
	else
	{
		nucSeg->RunClustering(true);
	}
	nucSeg->ReleaseSegMemory();
	outputImage = nucSeg->GetLabelImage();

	//Update For params actually used:
	definition->nuclearParameters.clear();
	std::vector<std::string> paramNames = nucSeg->GetParameterNames();
	for(int i=0; i<(int)paramNames.size(); ++i)
	{	
		ProjectDefinition::Parameter p;
		p.name = paramNames.at(i);
		p.value = nucSeg->GetParameter(p.name);
		definition->nuclearParameters.push_back(p);
	}

	delete nucSeg;

	std::cout << "Done NucSeg\n";

	//Calc Features:
	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(inputImage,outputImage,nucChannel,0);
	if(definition->intrinsicFeatures.size() > 0)
		iCalc->SetFeaturesOn( GetOnIntrinsicFeatures() );
	//iCalc->SetFeaturePrefix("nuc_");
	table = iCalc->Compute();									//Create a new table
	delete iCalc;

	std::cout << "Done NucFeats\n";

	resultIsEditable = true;
	return true;
}

bool ProjectProcessor::SegmentCytoplasm(int cytChannel, int memChannel)
{
	if(!inputImage || !outputImage || !table)
		return false;

	ftk::CytoplasmSegmentation * cytoSeg = new ftk::CytoplasmSegmentation();
	cytoSeg->SetDataInput(inputImage, "data_image", cytChannel,memChannel);
	cytoSeg->SetNucleiInput(outputImage, "label_image");		//Will append the result to outputImage

	for(int i=0; i<(int)definition->cytoplasmParameters.size(); ++i)
		cytoSeg->SetParameter(definition->cytoplasmParameters.at(i).name, int(definition->cytoplasmParameters.at(i).value));

	cytoSeg->Run();

	definition->cytoplasmParameters.clear();
	std::vector<std::string> paramNames = cytoSeg->GetParameterNames();
	for(int i=0; i<(int)paramNames.size(); ++i)
	{	
		ProjectDefinition::Parameter p;
		p.name = paramNames.at(i);
		p.value = cytoSeg->GetParameter(p.name);
		definition->cytoplasmParameters.push_back(p);
	}


	delete cytoSeg;

	std::cout << "Done CytoSeg\n";

	//Calc Features:
	if( cytChannel > -1 ){
		ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
		iCalc->SetInputImages(inputImage,outputImage,cytChannel,1);
		if(definition->intrinsicFeatures.size() > 0)
			iCalc->SetFeaturesOn( GetOnIntrinsicFeatures() );
		iCalc->SetFeaturePrefix("cyto_");
		iCalc->Append(table); //Append features to the table
		delete iCalc;
	}

	std::cout << "Done CytoFeats\n";

	resultIsEditable = false;
	return true;
}

bool ProjectProcessor::ComputeAssociations(void)
{
	if(definition->associationRules.size() == 0)
	{
		inputTypeNeeded = 3;
		return false;
	}
	
	for(int i=0; i<(int)definition->associationRules.size(); ++i)
	{
		ftk::AssociativeFeatureCalculator * assocCal = new ftk::AssociativeFeatureCalculator();
		assocCal->SetInputFile( definition->associationRules.at(i) );
		if(table)
		{
			assocCal->Append(table);
		}
		delete assocCal;
	}

	resultIsEditable = false;
	std::cout << "Done Associations\n";
	return true;
}

std::set<int> ProjectProcessor::GetOnIntrinsicFeatures(void)
{
	std::set<int> retSet;

	for(int f=0; f<IntrinsicFeatures::N; ++f)
	{
		std::string name = IntrinsicFeatures::Info[f].name;
		for(int p=0; p<(int)definition->intrinsicFeatures.size(); ++p)
		{
			if( definition->intrinsicFeatures.at(p) == name )
			{
				retSet.insert(f);
			}
		}
	}
	return retSet;
}

//************************************************************************
//************************************************************************
//************************************************************************
}  // end namespace ftk
