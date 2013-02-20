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
#ifndef __ftkProjectProcessor_h
#define __ftkProjectProcessor_h

#include "ftkProjectDefinition.h"
#include <ftkNuclearSegmentation.h>
#include <CytoplasmSegmentation/CytoplasmSegmentation.h>
#include <Nuclear_Association/ftkNuclearAssociationRules.h>
#include <ftkLabelImageToFeatures.h>
#include <ftkImage.h>
#include <ftkUtils.h>
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "ftkGUI/TrainingDialog.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGraphs/FTKgraph.h"
#include "ftkPreprocess2.h"
#include "PixelAnalysis/ftkPixelLevelAnalysis.h"
#include "SQLite/NESqlite/NESqliteFactory.h"
#include "PatternAnalysis/activeLearning/mclr.h"

#include "itkIntTypes.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMultiThreader.h"

// MODEL_SEG is defined as a compiler option in the Nucleus Editor's CMakeLists
#ifdef MODEL_SEG
#include "Statistical_Model_Segmentation/model_nucleus_seg.h"
#endif


namespace ftk
{

class ProjectProcessor
{
public:

	typedef unsigned char IPixelT;
	typedef unsigned short LPixelT;
	typedef itk::Image<unsigned char, 3> InputImageType;
	typedef itk::Image<unsigned short, 3> LabelImageType;

	typedef struct { ftk::ProjectDefinition::TaskType type; int inputChannel1; int inputChannel2; int inputChannel3; bool done; } Task;
	ProjectProcessor();

	//Steps to use
	void SetInputImage(ftk::Image::Pointer img){ inputImage = img; };
	void SetOutputImage(ftk::Image::Pointer img){ outputImage = img; };
	void SetTable(vtkSmartPointer<vtkTable> tab){ table = tab; };
	void SetDefinition(ftk::ProjectDefinition * def){ definition = def; };
	void SetPath( std::string path ){ save_path = path; };
	void Initialize(void);
	void SetExecPath( std::string inp_str ){ executable_path = inp_str; };
	void ProcessNext(void);						//Keep calling this until DoneProcessing is true;
	bool DoneProcessing(void){ return (lastTask == numTasks-1); };
	bool ReadyToEdit(void){ return resultIsEditable; };
	int NeedInput(void){ return inputTypeNeeded; };	//Return non-zero value indicating type of input needed
	void SetNumThreads( int threads );

	//Outputs
	ftk::Image::Pointer GetOutputImage(void){ return outputImage; };
	std::map< std::string, LabelImageType::Pointer > GetClassImageMap(void){ return classImageMap; };
	std::map< std::string, vtkSmartPointer<vtkTable> > GetClassCentroidMap(void){ return classCentroidMap; };
	vtkSmartPointer<vtkTable> GetTable(void){ return table; };
	vtkSmartPointer<vtkTable> GetNucAdjTable(void){ return NucAdjTable; }; 
	vtkSmartPointer<vtkTable> GetCellAdjTable(void){ return CellAdjTable; }; 

protected:
	//Tasks I can do:
	bool PreprocessImage( void );
	bool SegmentNuclei(int nucChannel);						//Segment nucChannel as Nuclei & Compute Intrinsic Features
	void mmSegmentation(int intChannel, int labChannel);
	bool ComputeFeatures(int nucChannel);
	template <typename InputPixelType, typename LabelPixelType>
			void ComputeMontageIntrinsicFeatures( int nucChannel );
	bool SegmentCytoplasm(int cytChannel, int memChannel);	//Segment Cytoplasm Channel & Compute Intrinsic Features
	bool ComputeAssociations(void);							//Compute Associative Measures
	bool PixLevAnalysis(void);								//If you must with all this nice object level machinery available
	bool Classify(void);									//Classify Cells
	bool Classify_mclr(void);
	bool Extract_Class(void);
	void ComputeAnalyteMeasures(void);						//Compute Analyte Measures by Class
	bool RunQuery(void);									//Run SQL Query

	std::set<int> GetOnIntrinsicFeatures(void);				//Return the list of intrinsic features to calculate

	ftk::Image::Pointer inputImage;
	ftk::Image::Pointer outputImage;
	std::map< std::string, LabelImageType::Pointer > classImageMap;
	std::map< std::string, vtkSmartPointer<vtkTable> > classCentroidMap;
	std::map<int, InputImageType::Pointer> original_image_map;
	ftk::ProjectDefinition * definition;
	std::vector<Task> tasks;
	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkTable> NucAdjTable;
	vtkSmartPointer<vtkTable> CellAdjTable;

private:
	char* stringToChar(std::string str);
	int numTasks;
	int lastTask;
	int n_thr;
	bool numThreadsSet;
	bool resultIsEditable;  //Only true when done nucleus segmentation!!
	int inputTypeNeeded;
	std::string save_path;
	std::string executable_path;
};

}  // end namespace ftk

#include "ftkProjectProcessor.txx"

#endif	// end __ftkProjectProcessor_h
