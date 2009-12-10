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
}

void ProjectProcessor::Initialize(void)
{
	if(!definition) return;

	tasks.clear();

	for(int i=0; i<definition->pipeline.size(); ++i)
	{
		Task t;
		t.type = definition->pipeline.at(i);
		switch(t.type)
		{
		case ProjectDefinition::NUCLEAR_SEGMENTATION:
			t.inputChannel = definition->FindInputChannel("NUCLEAR");
			break;
		case ProjectDefinition::CYTOPLASM_SEGMENTATION:
			t.inputChannel = definition->FindInputChannel("CYTOPLASM");
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
		taskDone = SegmentNuclei(tasks.at(thisTask).inputChannel);
		break;
	case ProjectDefinition::CYTOPLASM_SEGMENTATION:
		taskDone = SegmentCytoplasm(tasks.at(thisTask).inputChannel);
		break;
	case ProjectDefinition::RAW_ASSOCIATIONS:
		break;
	case ProjectDefinition::CLASSIFY:
		break;
	case ProjectDefinition::ANALYTE_MEASUREMENTS:
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
	nucSeg->SetParameters("");
	nucSeg->Binarize(false);
	nucSeg->DetectSeeds(false);
	nucSeg->RunClustering(true);
	nucSeg->Finalize();
	nucSeg->ReleaseSegMemory();
	outputImage = nucSeg->GetLabelImage();
	delete nucSeg;

	std::cout << "Done NucSeg\n";

	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(inputImage,outputImage,nucChannel,0);
	//iCalc->SetFeaturePrefix("nuc_");
	table = iCalc->Compute();									//Create a new table
	delete iCalc;

	std::cout << "Done NucFeats\n";

	return true;
}

bool ProjectProcessor::SegmentCytoplasm(int cytChannel)
{
	if(!inputImage || !outputImage || !table)
		return false;

	ftk::CytoplasmSegmentation * cytoSeg = new ftk::CytoplasmSegmentation();
	cytoSeg->SetDataInput(inputImage, "data_image", cytChannel);
	cytoSeg->SetNucleiInput(outputImage, "label_image");		//Will append the result to outputImage
	cytoSeg->Run();
	delete cytoSeg;

	std::cout << "Done CytoSeg\n";

	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(inputImage,outputImage,cytChannel,1);
	iCalc->SetFeaturePrefix("cyto_");
	iCalc->Append(table);										//Append features to the table
	delete iCalc;

	std::cout << "Done CytoFeats\n";

	return true;
}

void ProjectProcessor::ComputeAssociations(void)
{
	/*
	ftk::AssociativeFeatureCalculator * assocCal = new ftk::AssociativeFeatureCalculator();
	assocCal->SetInputFile(fileName.toStdString());

	if(!table)
	{
		table = assocCal->Compute();
		CreateNewTableWindow();
	}
	else
		assocCal->Append(table);
	*/
}

//************************************************************************
//************************************************************************
//************************************************************************
}  // end namespace ftk
