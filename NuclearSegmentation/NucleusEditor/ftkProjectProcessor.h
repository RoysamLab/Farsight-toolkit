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
#ifndef __ftkProjectProcessor_h
#define __ftkProjectProcessor_h

#include "ftkProjectDefinition.h"
#include <ftkNuclearSegmentation.h>
#include <CytoplasmSegmentation/CytoplasmSegmentation.h>
#include <Nuclear_Association/ftkNuclearAssociationRules.h>
#include <ftkLabelImageToFeatures.h>
#include <ftkImage.h>
#include <ftkUtils.h>

namespace ftk
{

class ProjectProcessor
{
public:
	typedef struct { ftk::ProjectDefinition::TaskType type; int inputChannel; bool done; } Task;
	ProjectProcessor();

	//Steps to use
	void SetInputImage(ftk::Image::Pointer img){ inputImage = img; }; 
	void SetDefinition(ftk::ProjectDefinition * def){ definition = def; };
	void Initialize(void);
	void ProcessNext(void);						//Keep calling this until DoneProcessing is true;
	bool DoneProcessing(void){ return (lastTask == numTasks-1); };

	//Outputs
	ftk::Image::Pointer GetOutputImage(void){ return outputImage; };
	vtkSmartPointer<vtkTable> GetTable(void){ return table; };
	
protected:
	//Tasks I can do:
	bool SegmentNuclei(int nucChannel);					//Segment nucChannel as Nuclei & Compute Intrinsic Features
	bool SegmentCytoplasm(int cytChannel);				//Segment Cytoplasm Channel & Compute Intrinsic Features
	void ComputeAssociations(void);						//Compute Associative Measures
	void Classify(void);								//Classify Cells
	void ComputeAnalyteMeasures(void);					//Compute Analyte Measures by Class

	ftk::Image::Pointer inputImage;
	ftk::Image::Pointer outputImage;
	ftk::ProjectDefinition * definition;
	std::vector<Task> tasks;
	vtkSmartPointer<vtkTable> table;

private:
	int numTasks;
	int lastTask;
};

}  // end namespace ftk

#endif	// end __ftkProjectProcessor_h
