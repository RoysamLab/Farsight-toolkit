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

#include <time.h>

namespace ftk
{
//Constructor
ProjectProcessor::ProjectProcessor()
{
	inputImage = NULL;
	outputImage = NULL;
	definition = NULL;
	tasks.clear();
	classImageMap.clear();
	classCentroidMap.clear();
	original_image_map.clear();
	table = NULL;
	NucAdjTable = NULL;
	CellAdjTable = NULL;
	numTasks = 0;
	lastTask = -1;
	resultIsEditable = false;
	inputTypeNeeded = 0;
	save_path = ".";
	n_thr = 4;//Temporary number for number of omp threads
	numThreadsSet = false;	
}

void ProjectProcessor::SetNumThreads( int threads )
{
	numThreadsSet = true;
	n_thr = threads;
#ifdef _OPENMP
	std::cout<<"Number of threads set to:"<<threads<<"\n";
#else
	std::cout<<"Number of threads set to:"<<threads<<" but openMP is not enabled\n";
#endif
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
		case ProjectDefinition::PREPROCESSING:
			break;
		case ProjectDefinition::NUCLEAR_SEGMENTATION:
			t.inputChannel1 = definition->FindInputChannel("NUCLEAR");
			break;
		case ProjectDefinition::FEATURE_COMPUTATION:
			t.inputChannel1 = definition->FindInputChannel("NUCLEAR");
			break;
		case ProjectDefinition::CYTOPLASM_SEGMENTATION:
			t.inputChannel2 = definition->FindInputChannel("CYTOPLASM");
			t.inputChannel3 = definition->FindInputChannel("MEMBRANE");
			break;
		case ProjectDefinition::RAW_ASSOCIATIONS:
			break;
		//case ProjectDefinition::MULTI_MODEL_SEGMENTATION:
		//	break;
		case ProjectDefinition::CLASSIFY:
			break;
		case ProjectDefinition::CLASSIFY_MCLR:
			break;
		case ProjectDefinition::CLASS_EXTRACTION:
			break;
		case ProjectDefinition::ANALYTE_MEASUREMENTS:
			break;
		case ProjectDefinition::PIXEL_ANALYSIS:
			break;
		case ProjectDefinition::QUERY:
			break;
		}
		t.done = false;
		tasks.push_back(t);
	}

	numTasks = (int)tasks.size();
	std::cout<<"Number of tasks to be done:"<<numTasks<<std::endl;
	lastTask = -1;
}

void ProjectProcessor::ProcessNext(void)
{
	if(DoneProcessing()) return;

	int thisTask = lastTask + 1;
	bool taskDone = false;
	switch(tasks.at(thisTask).type)
	{
	case ProjectDefinition::PREPROCESSING:
		taskDone = PreprocessImage();
		break;
	case ProjectDefinition::NUCLEAR_SEGMENTATION:
		taskDone = SegmentNuclei(tasks.at(thisTask).inputChannel1);
		break;
	case ProjectDefinition::FEATURE_COMPUTATION:
		taskDone = ComputeFeatures(tasks.at(thisTask).inputChannel1);
		break;
	case ProjectDefinition::CYTOPLASM_SEGMENTATION:
		taskDone = SegmentCytoplasm(tasks.at(thisTask).inputChannel2, tasks.at(thisTask).inputChannel3);
		break;
	case ProjectDefinition::RAW_ASSOCIATIONS:
		taskDone = ComputeAssociations();
		break;
	//case ProjectDefinition::MULTI_MODEL_SEGMENTATION:
	//	taskDone = mmSegmentation(tasks.at(thisTask).inputChannel1);
	//	break;
	case ProjectDefinition::ANALYTE_MEASUREMENTS:
		taskDone = false;
		break;
	case ProjectDefinition::CLASSIFY:
		taskDone = Classify();
		break;
	case ProjectDefinition::CLASSIFY_MCLR:
		taskDone = Classify_mclr();
		break;
	case ProjectDefinition::CLASS_EXTRACTION:
		taskDone = Extract_Class();
		break;
	case ProjectDefinition::PIXEL_ANALYSIS:
		taskDone = PixLevAnalysis();
		break;
	case ProjectDefinition::QUERY:
		taskDone = RunQuery();
		break;
	}

	if(taskDone)
	{
		tasks.at(thisTask).done = true;
		lastTask++;
	}
}

//***********************************************************************************************************
//  PRE-PROCESSING
//***********************************************************************************************************
bool ProjectProcessor::PreprocessImage()
{
	if(!inputImage)
		return false;

	if(definition->preprocessingParameters.size() == 0)
	{
		inputTypeNeeded = 3;
		return false;
	}

	for(std::vector<ftk::ProjectDefinition::preprocessParam>::iterator ppit=definition->preprocessingParameters.begin(); ppit!=definition->preprocessingParameters.end(); ++ppit ){
		int chNum = definition->FindInputChannel( ppit->channelName );
		if( chNum == -1 ){
			std::cerr<<"ERROR: Cannot find the channel "<<ppit->channelName<<" for preprocessing\n";
			continue;
		}

		original_image_map[chNum] = inputImage->GetItkPtr<unsigned char>(0,chNum);

		ftk::Preprocess *prep = new ftk::Preprocess( inputImage->GetItkPtr<unsigned char>(0,chNum) );
		TiXmlDocument tempXML;   
		TiXmlElement *genericName = new TiXmlElement( "Preprocess" );
		tempXML.LinkEndChild( genericName );
		//Write Filter:
		TiXmlElement * filterElement = new TiXmlElement(ppit->filterName.c_str());
		if( !ppit->paramenter1.empty() ){
			filterElement->SetAttribute( ppit->paramenter1.c_str(), ftk::NumToString( ppit->value1 ) );
		}
		if( !ppit->paramenter2.empty() ){
			filterElement->SetAttribute( ppit->paramenter2.c_str(), ftk::NumToString( ppit->value2 ) );
		}
		if( !ppit->paramenter3.empty() ){
			filterElement->SetAttribute( ppit->paramenter3.c_str(), ftk::NumToString( ppit->value3 ) );
		}
		if( !ppit->paramenter4.empty() ){
			filterElement->SetAttribute( ppit->paramenter4.c_str(), ftk::NumToString( ppit->value4 ) );
		}
		if( !ppit->paramenter5.empty() ){
			filterElement->SetAttribute( ppit->paramenter5.c_str(), ftk::NumToString( ppit->value5 ) );
		}
		if( !ppit->paramenter6.empty() ){
			filterElement->SetAttribute( ppit->paramenter6.c_str(), ftk::NumToString( ppit->value6 ) );
		}
		genericName->LinkEndChild(filterElement);
		std::string tempFilename;
		tempFilename = save_path + ppit->channelName + ppit->filterName;
		tempFilename.append( ".xml" );
		tempXML.SaveFile( tempFilename.c_str() );
		prep->RunPipe( tempFilename );
		const ftk::Image::Info * Iinfo = inputImage->GetImageInfo();
		int byte_to_cpy = Iinfo->numZSlices * Iinfo->numRows * Iinfo->numColumns * Iinfo->bytesPerPix;
		memcpy( inputImage->GetDataPtr(0,chNum), prep->GetImage()->GetBufferPointer(), byte_to_cpy );
		delete prep;
	}
	return true;
}


//***********************************************************************************************************
//  NUCLEAR SEGMENTATION
//***********************************************************************************************************
bool ProjectProcessor::SegmentNuclei(int nucChannel)
{
	clock_t start_time = clock();
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
	const ftk::Image::Info *info = inputImage->GetImageInfo();
	if(info->numTSlices==1)
	{
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
	}
	else
	{
		nucSeg->SegmentAllTimes(finalize);
		outputImage = nucSeg->GetLabelImage();
	}

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

	if((int)definition->mmSegFiles.size() != 0)
	{		
		mmSegmentation(nucChannel, 0);
	}

	const ftk::Image::Info * image_info = inputImage->GetImageInfo();
	int byte_to_cpy = image_info->numZSlices * image_info->numRows * image_info->numColumns * image_info->bytesPerPix;
	std::map<int, InputImageType::Pointer>::iterator it;
	for ( it = original_image_map.begin() ; it != original_image_map.end(); ++it )
	{
		int chNum = (*it).first;
		memcpy( inputImage->GetDataPtr(0,chNum), ((*it).second)->GetBufferPointer(), byte_to_cpy );
	}
	original_image_map.clear();


	delete nucSeg;
	cout << "Total time to segmentation is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << endl;
	std::cout << "Done Nucleus Segmentation\nComputing intrisic features for the nuclei\n";

	if(info->numTSlices==1)
	{
		//Calc Features:
		ComputeFeatures( nucChannel );
	}

	resultIsEditable = true;
	
	return true;
}

//***********************************************************************************************************
//  NUCLEAR SEGMENTATION FOR MONTAGES
//***********************************************************************************************************



//***********************************************************************************************************
//  MULTI-MODEL SEGMENTATION
//***********************************************************************************************************
void ProjectProcessor::mmSegmentation(int intChannel, int labChannel)
{
	#ifdef MODEL_SEG
	{
		model_nucleus_seg *MNS = new model_nucleus_seg();
		MNS->SetImages(inputImage->GetItkPtr<IPixelT>(0,intChannel), outputImage->GetItkPtr<LPixelT>(0,labChannel));
		for(int i=0; i<(int)definition->mmSegFiles.size(); ++i)
		{
			if(definition->mmSegFiles.at(i).type == "MMSeg_Params")
				MNS->LoadSegParams(definition->mmSegFiles.at(i).path);
			if(definition->mmSegFiles.at(i).type == "Training_File")
				MNS->SetTrainingFile(stringToChar(definition->mmSegFiles.at(i).path));
		}

		if(MNS->getasscofeatnumb()>0)
			MNS->Associations_NucEd(inputImage, definition->associationRules);
		std::vector<unsigned short> US_Cell_ids;
		US_Cell_ids.clear();

		if(MNS->SPLIT)
			US_Cell_ids = MNS->Detect_undersegmented_cells();

		if(US_Cell_ids.size()>0)
		{
			MNS->bImage = MNS->SplitImage(US_Cell_ids,MNS->inputImage,MNS->bImage);
			MNS->splitflag = 1;
		}		

		typedef ftk::LabelImageToFeatures< unsigned char,  unsigned short, 3 > FeatureCalcType;
		typedef ftk::IntrinsicFeatures FeaturesType;

		//Calculates all the features and stores the individual 
		//labeled and raw images
		MNS->GetFeatsnImages();

		////Compute the associations after the split for split Ids only
		if(MNS->splitflag==1 && MNS->getasscofeatnumb()>0 )
			MNS->splitAssociations();			

		std::cout<<"Computing Scores"<<std::endl;

		//Iniscores contains the inital scores of all the fragments
		std::vector<double> IniScores = MNS->GetOriginalScores();

		std::cout<<"Finished Computing Scores"<<std::endl;

		//	 idlist will contain the list of all the ids for which the 
		//	 MergeTrees have to be built and analyzed.

		std::vector<unsigned short> idList = MNS->getListofIds(); 	
		std::vector< unsigned short > labelIndex = MNS->getLabelIndex();
		std::vector<FeaturesType> allFeat = MNS->getFeats();


		//Loop through the list of ids and build graphs.
		ftkgnt *sg1 = new ftkgnt();
		sg1->runLabFilter(MNS->inputImage,MNS->bImage);
		sg1->setFeats(allFeat,labelIndex);

		// 1e5 is the default volume I choose. With this
		// default value, the size of the tree is controlled 
		// only by the depth value and volume has no effect
		// If this limit is specified in the parameters file 
		// that value takes precedence

		if(MNS->MAX_VOL!=1e5)
			sg1->setmaxVol(MNS->MAX_VOL);

		// 6 is the default depth I choose. If specidied differently 
		// in the parameters file, then use that value
		if(MNS->MAX_DEPTH!=6)
			sg1->setmaxDepth(MNS->MAX_DEPTH);


		for (unsigned short r = 0;r < idList.size();r++)
		{		
			FeatureCalcType::LabelPixelType id = idList[r]; 
			cout<<"\rEvaluating Hypothesis for id "<<id;
			sg1->RAG = sg1->BuildRAG(id); 
			//Build the Merge Tree from the RAG and root information
			ftkgnt::MTreeType mTree =  sg1->BuildMergeTreeDcon(sg1->RAG,id,sg1->hypotheses);	
			MNS->GetScoresfromKPLS(mTree);
			//std::cout<<"\rr"<<r;
		}

		std::cout<<""<<std::endl;
		MNS->SelectHypothesis();

		//PERFORM MERGES 
		itk::Image< unsigned short, 3>::Pointer mmSegImage = MNS->PerformMerges_NucEd();

		const ftk::Image::Info * info = outputImage->GetImageInfo();
		int bpChunk = info->numZSlices * info->numRows * info->numColumns * info->bytesPerPix;
		memcpy( outputImage->GetDataPtr(0,labChannel), mmSegImage->GetBufferPointer(), bpChunk );

	}

	#endif

}

char* ProjectProcessor::stringToChar(std::string str)
{
	char * writable = new char[str.size() + 1];
	std::copy(str.begin(), str.end(), writable);
	writable[str.size()] = '\0'; 
	return writable;
}


//***********************************************************************************************************
//  INTRINSIC FEATURE COMPUTATION
//***********************************************************************************************************
bool ProjectProcessor::ComputeFeatures(int nucChannel)
{

#ifdef _OPENMP
	//Compute features in parallel instantiated for uchar and ushort intensity images
	//and ushort and uint label images
	if( inputImage->GetImageInfo()->dataType == itk::ImageIOBase::UCHAR )
	{
		if( outputImage->GetImageInfo()->dataType == itk::ImageIOBase::USHORT )
			ComputeMontageIntrinsicFeatures
				< unsigned char, unsigned short >( nucChannel );
		else if( outputImage->GetImageInfo()->dataType == itk::ImageIOBase::UINT )
			ComputeMontageIntrinsicFeatures
				< unsigned char, unsigned int >( nucChannel );
		else
		{
			std::cout<<"Feature computation can only handle ushort "
				<<"and uint labels";
			return false;
		}
	}
	else if( inputImage->GetImageInfo()->dataType == itk::ImageIOBase::USHORT )
	{
		if( outputImage->GetImageInfo()->dataType == itk::ImageIOBase::USHORT )
			ComputeMontageIntrinsicFeatures
				< unsigned short, unsigned short >( nucChannel );
		else if( outputImage->GetImageInfo()->dataType == itk::ImageIOBase::UINT )
			ComputeMontageIntrinsicFeatures
				< unsigned short, unsigned int >( nucChannel );
		else
		{
			std::cout<<"Feature computation can only handle ushort "
				<<"and uint labels";
			return false;
		}
	}
	else
	{
			std::cout<<"Feature computation can only handle uchar "
				<<"and ushort images";
			return false;
	}

#else
	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(inputImage,outputImage,nucChannel,0);
	if(definition->intrinsicFeatures.size() > 0)
		iCalc->SetFeaturesOn( GetOnIntrinsicFeatures() );
	//iCalc->SetFeaturePrefix("nuc_");
	table = iCalc->Compute();							//Create a new table
	delete iCalc;
	std::cout << "Done: Instrinsic Nuclear Features\n";
#endif
	resultIsEditable = true;
	return true;
}

//***********************************************************************************************************
//  CYTOPLASM SEGMENTATION
//***********************************************************************************************************
bool ProjectProcessor::SegmentCytoplasm(int cytChannel, int memChannel)
{
	if(!inputImage || !outputImage || !table)
		return false;

	std::cout << "Begin Cytoplasm Segmentation\n";

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

	std::cout << "Done: Cytoplasm Segmentation\nComputing intrisic features for the whole cell\n";

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

	std::cout << "Done: Intrisic features for the whole cell\n";

	resultIsEditable = false;
	return true;
}


//***********************************************************************************************************
//  ASSOCIATIVE FEATURE COMPUTATION
//***********************************************************************************************************
bool ProjectProcessor::ComputeAssociations(void)
{
	if(!inputImage || !outputImage || !table)
		return false;

	if(definition->associationRules.size() == 0)
	{
		inputTypeNeeded = 3;
		return false;
	}

	for(std::vector<ftk::AssociationRule>::iterator ascit=definition->associationRules.begin(); ascit!=definition->associationRules.end(); ++ascit )
	{
		int seg_channel_number=-1;
		int inp_channel_number=-1;
		if( strcmp( ascit->GetSegmentationFileName().c_str(), "NUCLEAR" ) == 0 )
		{
			if( !outputImage->GetImageInfo()->numChannels )
			{
				std::cout<<"Unable to access labeled nuclear image while computing associative features\n";
				return false;
			}
			seg_channel_number = 0;			
		} 
		else if( strcmp( ascit->GetSegmentationFileName().c_str(), "CYTOPLASM" ) == 0 )
		{
			if( outputImage->GetImageInfo()->numChannels < 2 )
			{
				std::cout<<"Unable to access labeled cytoplasmic image while computing associative features\n";
				return false;
			}
			seg_channel_number = 1;
		}
		else
		{
			std::cout<<"Please check region type for associative feature computation\n";
			return false;
		}

		for( int j=0; j<(int)inputImage->GetImageInfo()->channelNames.size(); ++j )
		{
			if( strcmp( inputImage->GetImageInfo()->channelNames.at(j).c_str(), ascit->GetTargetFileNmae().c_str() ) == 0 )
			{
				inp_channel_number=j;
				break;
			}
		}

		if( inp_channel_number == -1 )
		{
			std::cout<<"Unable to access grayscale image while computing associative feature: "<<ascit->GetRuleName()<<std::endl;
			return false;
		}

		(*ascit).set_path( save_path );

		ftk::AssociativeFeatureCalculator * assocCal = new ftk::AssociativeFeatureCalculator();
		assocCal->SetInputs(inputImage, inp_channel_number, outputImage, seg_channel_number, &(*ascit) );
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



//***********************************************************************************************************
//  PIXEL LEVEL ANALYSIS
//***********************************************************************************************************
bool ProjectProcessor::PixLevAnalysis(void){
	if( !inputImage )
		return false;

	if( !definition->pixelLevelRules.size() )
		return false;

	bool success_run=false;
	for(std::vector<ftk::PixelAnalysisDefinitions>::iterator pait=definition->pixelLevelRules.begin(); pait!=definition->pixelLevelRules.end(); ++pait ){
		ftk::PixelLevelAnalysis *PAn = new ftk::PixelLevelAnalysis();
		if( ((*pait).mode == 1) || ((*pait).mode == 5) ){
			PAn->SetInputs( (*pait).regionChannelName, (*pait).targetChannelName, (*pait).outputFilename, 0, (*pait).mode, (*pait).erodeRadius );
			success_run = PAn->RunAnalysis1();
		}
		else if( ((*pait).mode == 2) || ((*pait).mode == 6) ){
			PAn->SetInputs( (*pait).regionChannelName, (*pait).targetChannelName, (*pait).outputFilename, (*pait).radius, (*pait).mode, (*pait).erodeRadius );
			success_run = PAn->RunAnalysis2();
		}
		else if( ((*pait).mode == 3) || ((*pait).mode == 7) || ((*pait).mode == 4) || ((*pait).mode == 8) ){
			PAn->SetInputs( (*pait).regionChannelName, (*pait).targetChannelName, (*pait).outputFilename, (*pait).radius, (*pait).mode, (*pait).erodeRadius );
			success_run = PAn->RunAnalysis3();
		}
		else if( !success_run ){
			std::cerr<<"ERROR: Run Failed, Check Definitions\n";
		}
		else{
			std::cerr<<"ERROR: Check Pixel Anaysis Mode\n";
		}
		delete PAn;
	}
	return success_run;
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


//***********************************************************************************************************
//  CLASSIFICATION
//***********************************************************************************************************
bool ProjectProcessor::Classify(void){
	std::cout<<"Entered Classification step\n";
	if(!table)
		return false;
	std::cout<<"Table present strating classification\n";
	TrainingDialogNoGUI *d = new TrainingDialogNoGUI(table);
	d->loadModelFromFile1(definition->classificationTrainingData);
	delete d;
	std::cout<<"Model Loaded\n";
	for(int j=0; j<(int)definition->classificationParameters.size(); ++j){
		bool training_col_found = false;
		std::vector<int> KPLsColumnsToUse;
		std::string output_col_name;
		for( int i=0; i<table->GetNumberOfColumns(); ++i ){
			std::string current_column;
			current_column = table->GetColumnName(i);
			if( strcmp (current_column.c_str(),definition->classificationParameters.at(j).TrainingColumn.c_str()) == 0 ){
				std::string::iterator it;
				it=current_column.begin();
				current_column.erase ( current_column.begin(), current_column.begin()+6 );
				output_col_name = "prediction_" + current_column;
				training_col_found = true;
			}
			else{
				for( int k=0; k<(int)definition->classificationParameters.at(j).ClassificationColumns.size(); ++k ){
					if( strcmp (current_column.c_str(),definition->classificationParameters.at(j).ClassificationColumns.at(k).c_str()) == 0 ){
						KPLsColumnsToUse.push_back( i );
					}
				}
			}
		}
		std::cout<<"Running classifier:"<<j+1<<"\n";
		if( training_col_found && !KPLsColumnsToUse.empty() ){
			PatternAnalysisWizardNoGUI *p = new PatternAnalysisWizardNoGUI(this->table,definition->classificationParameters.at(j).TrainingColumn.c_str(), output_col_name.c_str() );
			p->KPLSrun1(KPLsColumnsToUse);
			delete p;
		}
		else{
			std::cerr<<"Check classification parameter:"<<j<<std::endl;
		}
	}
	return true;
}


/*void ProjectProcessor::sqliteCloseConnection(sqlite3 *dbConn)
{   
	//  sqlite3_close() routine is the destructor for the sqlite3 object. 
	//  Calls to sqlite3_close() return SQLITE_OK if the sqlite3 object is successfullly destroyed 
	//  and all associated resources are deallocated.
	
	sqlite3_close(dbConn);
	fprintf(stderr,"Execution complete!\n");

}
*/

//***********************************************************************************************************
//  CLASSIFICATION USING MCLR
//***********************************************************************************************************
bool ProjectProcessor::Classify_mclr(void)
{
	std::cout<<"Classification Started Using MCLR\n";
	if(!table)
		return false;

	for(int i=0; i<(int)definition->Classification_Rules.size(); ++i)
	{
		vtkSmartPointer<vtkTable> active_model_table = ftk::LoadTable(definition->Classification_Rules[i].TrainingFileName);
		
		// to generate the Active Learning Matrix
		vnl_matrix<double> act_learn_matrix;
		act_learn_matrix.set_size((int)active_model_table->GetNumberOfColumns() , (int)active_model_table->GetNumberOfRows() - 2);
		for(int row = 2; row<(int)active_model_table->GetNumberOfRows(); ++row)
		{
			for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
			{
				act_learn_matrix.put(col, row-2, active_model_table->GetValue(row,col).ToDouble());
			}
		}

		//to generate the std_deviation and the mean vectors
		vnl_vector<double> std_dev_vec, mean_vec; 
		std_dev_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
		mean_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
		for(int col=1; col<(int)active_model_table->GetNumberOfColumns(); ++col)
		{
			std_dev_vec.put(col-1, active_model_table->GetValue(0,col).ToDouble());
			mean_vec.put(col-1, active_model_table->GetValue(1,col).ToDouble());
		}
		active_model_table->RemoveRow(0);
		active_model_table->RemoveRow(0);
		active_model_table->RemoveColumn(0);
		
		std::string classification_name = definition->Classification_Rules[i].ClassColName;
		double confidence_thresh = definition->Classification_Rules[i].ConfThreshold;

		MCLR* mclr = new MCLR();
		// Number of features and classes needed in "add_bias" fuction of MCLR
		mclr->Set_Number_Of_Classes((int)active_model_table->GetNumberOfRows());
		mclr->Set_Number_Of_Features((int)active_model_table->GetNumberOfColumns());

		vtkSmartPointer<vtkTable> test_table  = vtkSmartPointer<vtkTable>::New();
		test_table->Initialize();
		test_table->SetNumberOfRows(table->GetNumberOfRows());
		for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName(active_model_table->GetColumnName(col));
			test_table->AddColumn(column);	
		}
		for(int row = 0; row < (int)table->GetNumberOfRows(); ++row)
		{		
			vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
			for(int c=0; c<(int)test_table->GetNumberOfColumns();++c)
				model_data1->InsertNextValue(table->GetValueByName(row,test_table->GetColumnName(c)));
			test_table->InsertNextRow(model_data1);
		}	
		
		////// Final Data  to classify from the model
		vnl_matrix<double> data_classify;
		data_classify =  mclr->Normalize_Feature_Matrix_w(mclr->tableToMatrix_w(test_table), std_dev_vec, mean_vec);
		data_classify = data_classify.transpose();

		vnl_matrix<double> currprob;
		currprob = mclr->Test_Current_Model_w(data_classify, act_learn_matrix);

		std::string prediction_col_name = "prediction_active_" + classification_name;
		std::string confidence_col_name = "confidence_" + classification_name;

		//// Add the Prediction Column 
		std::vector< std::string > prediction_column_names = ftk::GetColumsWithString(prediction_col_name, table );
		if(prediction_column_names.size() == 0)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName(prediction_col_name.c_str());
			column->SetNumberOfValues( table->GetNumberOfRows() );
			table->AddColumn(column);
		}

		// Add the confidence column
		std::vector< std::string > confidence_column_names = ftk::GetColumsWithString(confidence_col_name, table );
		if(confidence_column_names.size() == 0)
		{
			vtkSmartPointer<vtkDoubleArray> column_confidence = vtkSmartPointer<vtkDoubleArray>::New();
			column_confidence->SetName(confidence_col_name.c_str());
			column_confidence->SetNumberOfValues( table->GetNumberOfRows() );
			table->AddColumn(column_confidence);
		}
		
		for(int row = 0; row<(int)table->GetNumberOfRows(); ++row)  
		{
			vnl_vector<double> curr_col = currprob.get_column(row);
			table->SetValueByName(row, confidence_col_name.c_str(), vtkVariant(curr_col(curr_col.arg_max())));
			if(curr_col(curr_col.arg_max()) > confidence_thresh) 
			{
				table->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(curr_col.arg_max()+1));						
			}
			else
			{
				table->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(0));
			}
		}
	}

	std::cout<< "Classification done" << std::endl;
	
	return true;
}

//***********************************************************************************************************
//  EXTRACT CLASS
//***********************************************************************************************************
bool ProjectProcessor::Extract_Class(void)
{
	for(int i=0; i<(int)definition->Class_Extraction_Rules.size(); ++i)
	{
		bool found = false;
		std::string className = definition->Class_Extraction_Rules[i].Class_Name;
		std::string classColumnName = "prediction_active_" + className;
		int class_ext = definition->Class_Extraction_Rules[i].Class;
		std::cout<<"Extracting from classifier "<<classColumnName<<" class:"<<class_ext<<std::endl<<std::flush;

		for(int col=((int)table->GetNumberOfColumns())-1; col>=0; --col)
		{	
			std::string current_column = table->GetColumnName(col);
			if(current_column.find(classColumnName) != std::string::npos )
			{
				found = true;
				break;			
			}	
		}
		if(!found) continue;


		//Extract the Image of desired class
		LabelImageType::Pointer classImage = LabelImageType::New();
		LabelImageType::Pointer labelImage = outputImage->GetItkPtr<unsigned short>(0,0);
		LabelImageType::PixelType * labelArray = labelImage->GetBufferPointer();

		itk::Size<3> im_size = labelImage->GetBufferedRegion().GetSize();
		LabelImageType::IndexType start;
		start[0] =   0;  // first index on X
		start[1] =   0;  // first index on Y    
		start[2] =   0;  // first index on Z  
		LabelImageType::PointType origin;
		origin[0] = 0; 
		origin[1] = 0;    
		origin[2] = 0;    
		classImage->SetOrigin( origin );
		LabelImageType::RegionType region;
		region.SetSize( im_size );
		region.SetIndex( start );
		classImage->SetRegions( region );
		classImage->Allocate();
		classImage->FillBuffer(0);
		classImage->Update();
		LabelImageType::PixelType * classArray = classImage->GetBufferPointer();

		itk::SizeValueType slice_size = im_size[1] * im_size[0];
		itk::SizeValueType row_size = im_size[0];
		std::map<unsigned short, int> classMap;
		for(unsigned int row=0; row<table->GetNumberOfRows(); ++row)
		{
			classMap[table->GetValue(row,0).ToUnsignedShort()] = table->GetValueByName(row, classColumnName.c_str()).ToInt();
		}

		#pragma omp parallel for
		for(itk::IndexValueType i=0; i<im_size[2]; ++i)
		{
			for(itk::IndexValueType j=0; j<im_size[1]; ++j)
			{
				for(itk::IndexValueType k=0; k<im_size[0]; ++k)
				{
					itk::IndexValueType offset = (i*slice_size)+(j*row_size)+k;
					if(classMap[labelArray[offset]] == class_ext)
						classArray[offset] = labelArray[offset];
					else
						classArray[offset] = 0;
				}
			}
		}
		classImageMap[className] = classImage;

		//Extract the Centroids of desired class
		vtkSmartPointer<vtkTable> centroid_table  = vtkSmartPointer<vtkTable>::New();
		centroid_table->Initialize();
		vtkSmartPointer<vtkDoubleArray> column_x = vtkSmartPointer<vtkDoubleArray>::New();
		column_x->SetName("centroid_x");
		centroid_table->AddColumn(column_x);	
		vtkSmartPointer<vtkDoubleArray> column_y = vtkSmartPointer<vtkDoubleArray>::New();
		column_y->SetName("centroid_y");
		centroid_table->AddColumn(column_y);	
		vtkSmartPointer<vtkDoubleArray> column_z = vtkSmartPointer<vtkDoubleArray>::New();
		column_z->SetName("centroid_z");
		centroid_table->AddColumn(column_z);	
		for(unsigned int row = 0; row<table->GetNumberOfRows(); ++row)
		{		
			vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
			if(table->GetValueByName(row, classColumnName.c_str()).ToInt() == class_ext)
			{
				model_data1->InsertNextValue(table->GetValueByName(row, "centroid_x"));
				model_data1->InsertNextValue(table->GetValueByName(row, "centroid_y"));
				model_data1->InsertNextValue(table->GetValueByName(row, "centroid_z"));
				centroid_table->InsertNextRow(model_data1);
			}			
		}	
		classCentroidMap[className] = centroid_table;
	}


	return true;
}



//***********************************************************************************************************
//  RUN QUERIES
//***********************************************************************************************************
bool ProjectProcessor::RunQuery(void)
{
	int sql_db_img_id=-1;
	for(int j=0; j<(int)definition->queryParameters.size(); ++j){
		sqlite3 *dbConn;
	
		//  This sql string can be passed to this function from any UI module in Farsight.
		//  1) Open connection
		//  2) Execute query
		//  3) Dynamic binding od results to UI
		//  4) Close connection.

		char * sql;
    
		//   int sqlite3_open(
		//   const char *filename,  Database filename (UTF-8) 
		//   sqlite3 **ppDb         OUT: SQLite db handle )
		std::string db_name = executable_path + "\\NE.s3db";
		std::cout<<"Opening database: "<<db_name<<std::endl;
		dbConn = ftk::sqliteGetConnection( db_name.c_str() );
		//dbConn = ftk::sqliteOpenConnection();
		if( dbConn ){
			if(j==0){
				std::vector<std::string> column_names;
				//std::string temp1,temp2;
				//temp1 = "IMG_ID"; temp2 = "CELL_ID";
				//column_names.push_back( temp1 );
				//column_names.push_back( temp2 );
					for (int col = 1; col< table->GetNumberOfColumns(); ++col){
					std::string temp3=table->GetColumnName(col);
					column_names.push_back(temp3);
				}
				std::string table_name;
				table_name = "IMAGE_TEST";
				ftk::checkForUpdate( dbConn, column_names );

				std::string image_name;
				image_name = save_path + definition->name;
				char *im_nm_cstr = new char [image_name.size()+1];
				strcpy (im_nm_cstr, image_name.c_str());
				char *path_nm_cstr = new char [save_path.size()+1];
				strcpy (path_nm_cstr, save_path.c_str());
				std::vector< double > table_array;
				for (int row = 0; row< table->GetNumberOfRows(); ++row){
					for (int col = 0; col< table->GetNumberOfColumns(); ++col){
						table_array.push_back(table->GetValue(row,col).ToDouble());
					}
				}
				sql_db_img_id = ftk::GenericInsert( dbConn, im_nm_cstr, table_name.c_str(), path_nm_cstr, table_array,table->GetNumberOfColumns(), table->GetNumberOfRows(), column_names );	

				//sql = "select * from IMAGE;";
				//sql = "select * from IMAGE_TEST where IMG_ID = 1;";
				//sql = "select * from IMAGE where IMG_ID = 1;";
				//sql = "select * from IMAGE_TEST where IMG_ID = 1 and CELL_ID = 2;";
				//sql = "select * from IMAGE where IMG_ID = 2 and cell_id = 3 and eccentricity = 2.24 ;";
			}
			if( strcmp( definition->queryParameters.at(j).name.c_str(), "count" ) == 0 ){
				std::string temp = "SELECT COUNT IMG_ID FROM ( " + definition->queryParameters.at(j).value + " AND IMG_ID = " + ftk::NumToString(sql_db_img_id) + ");";
				sql = new char [temp.size()+1];
				strcpy (sql, temp.c_str());
			}
			else{
				std::string temp = definition->queryParameters.at(j).value + " AND IMG_ID = " + ftk::NumToString(sql_db_img_id) + ";";
					sql = new char [temp.size()+1];
				strcpy (sql, temp.c_str());
			}

			std::cout<<sql<<"\t---this is my query" <<std::endl;

			//Generic API.Can be called from any module
			ftk::sqliteExecuteQuery2(dbConn,sql,save_path,j+1);

			//Close the connetion.Free allocated resources
			ftk::sqliteCloseConnection(dbConn);
			return true;
		}
	}
	return false;
}


//************************************************************************
//************************************************************************
//************************************************************************
}  // end namespace ftk
