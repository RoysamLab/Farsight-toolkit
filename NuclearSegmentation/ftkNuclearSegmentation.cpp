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
#define Pi 3.14159265

#include "ftkNuclearSegmentation.h"
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "yousef_core/graphColLearn_3D/sequential_coloring.cpp"

namespace ftk 
{

//Constructor
NuclearSegmentation::NuclearSegmentation()
{
	NucleusSeg = NULL;
	this->ResetAll();

	paramNames.push_back("high_sensitivity");
	paramNames.push_back("adaptive_binarization");
	paramNames.push_back("LoG_size");
	paramNames.push_back("min_scale");
	paramNames.push_back("max_scale");
	paramNames.push_back("xy_clustering_res");
	paramNames.push_back("z_clustering_res");
	paramNames.push_back("finalize_segmentation");
	paramNames.push_back("sampling_ratio_XY_to_Z");
	paramNames.push_back("Use_Distance_Map");
	paramNames.push_back("refinement_range");
	paramNames.push_back("min_object_size");
	currentTime = 0;
}

NuclearSegmentation::~NuclearSegmentation()
{
	this->ReleaseSegMemory();
}

//*****************************************************************************
// Reset all variables and clear all memory (private)
//*****************************************************************************
void NuclearSegmentation::ResetAll(void)
{
	dataFilename.clear();
	dataImage = NULL;
	channelNumber = 0;
	labelImage = NULL;
	bBoxMap.clear();
	centerMap.clear();
	paramFilename.clear();
	myParameters.clear();
	EditsNotSaved = false;
	this->ReleaseSegMemory();
}

void NuclearSegmentation::ReleaseSegMemory()
{
	if(NucleusSeg)
	{
		delete NucleusSeg;
		lastRunStep = 0;
		NucleusSeg = NULL;
	}
}

bool NuclearSegmentation::LoadInput(std::string fname, int chNumber)
{
	ftk::Image::Pointer tmpImg = ftk::Image::New();
	if(!tmpImg->LoadFile(fname))	//Load for display
	{
		errorMessage = "Data Image failed to load";
		return false;
	}
	return this->SetInput(tmpImg,fname,chNumber);
}

bool NuclearSegmentation::SetInput(ftk::Image::Pointer inImg, std::string fname, int chNumber)
{
	if(chNumber > inImg->GetImageInfo()->numChannels)
	{
		errorMessage = "channel does not exist";
		return false;
	}
	if(inImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
	{
		errorMessage = "module only works for 8-bit unsigned char input data";
		return false;
	}
	this->dataFilename = fname;
	this->dataImage = inImg;
	this->channelNumber = chNumber;
	return true;
}

void NuclearSegmentation::SetParameters(std::string paramfile)
{
	this->paramFilename = paramfile;
}

void NuclearSegmentation::SetParameter(std::string name, int value)
{
	bool found = false;
	for(int i=0; i<(int)myParameters.size(); ++i)
	{
		if(myParameters.at(i).name == name)
		{
			myParameters.at(i).value = value;
			found = true;
			break;
		}
	}

	if(!found)		//Add the parameter
	{
		Parameter newP;
		newP.name = name;
		newP.value = value;
		myParameters.push_back(newP);
	}
}

int NuclearSegmentation::GetParameter(std::string name)
{
	for(int i=0; i<(int)myParameters.size(); ++i)
	{
		if(myParameters.at(i).name == name)
		{
			return myParameters.at(i).value;
		}
	}
	return -1;
}

void NuclearSegmentation::ConvertParameters(int params[12])
{
	//These are defaults:
	params[0]=0;
	params[1]=0;
	params[2]=30;
	params[3]=5;
	params[4]=8;
	params[5]=5;
	params[6]=2;
	params[7]=1;
	params[8]=2;
	params[9]=1;
	params[10]=6;
	params[11]=100;

	for(int i=0; i<(int)myParameters.size(); ++i)
	{
		if(myParameters.at(i).name == paramNames.at(0))
			params[0] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(1))
			params[1] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(2))
			params[2] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(3))
			params[3] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(4))
			params[4] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(5))
			params[5] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(6))
			params[6] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(7))
			params[7] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(8))
			params[8] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(9))
			params[9] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(10))
			params[10] = myParameters.at(i).value;
		if(myParameters.at(i).name == paramNames.at(11))
			params[11] = myParameters.at(i).value;
	}
}

//*****
// If you are not going to need to return the resulting image as ftk::Image, then pass false to this function
//*****
bool NuclearSegmentation::Binarize(bool getResultImg)
{
	if(!dataImage)
	{
		errorMessage = "No data loaded";
		return false;
	}

	const Image::Info *info = dataImage->GetImageInfo();
	int numStacks = info->numZSlices;
	int numRows = info->numRows;				//y-direction
	int numColumns = info->numColumns; 			//x-direction

	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,channelNumber,0);	//Expects grayscale image // Amin: should include a for loop

	if(NucleusSeg) delete NucleusSeg;
	NucleusSeg = new yousef_nucleus_seg();

	int adap_bin = 0;
	if(myParameters.size() == 0)
	{
		NucleusSeg->readParametersFromFile("");		//Will use automatic parameter detection	
	}
	else
	{
		int params[12];
		ConvertParameters(params);
		NucleusSeg->setParams(params); //Will use the parameters that have been set and defauls for the rest
		adap_bin = params[1];
	}

	NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
	if(adap_bin == 1)
	{
		itk::Image<unsigned char, 3>::Pointer binImage = Adaptive_Binarization(dataImage->GetItkPtr<unsigned char>(0,channelNumber));
		itk::Image<unsigned char, 3>::PixelType *binArray = binImage->GetBufferPointer();
		unsigned long long imgArrayLength = numStacks*numRows*numColumns;
		unsigned short *binPtr = new unsigned short[imgArrayLength];
		for(unsigned long long i=0; i<imgArrayLength; ++i)
		{
			binPtr[i] = binArray[i];
		}
		NucleusSeg->setBinImage(binPtr);
	}
	else
		NucleusSeg->runBinarization();
	lastRunStep = 1;

	//Get the output
	if(getResultImg)
		return GetResultImage();
	return true;
}

bool NuclearSegmentation::DetectSeeds(bool getResultImg)
{
	if(!NucleusSeg || lastRunStep < 1)
	{
		errorMessage = "No Binarization";
		return false;
	}
	NucleusSeg->runSeedDetection();
	lastRunStep = 2;

	//Get output image
	if(getResultImg)
		return GetResultImage();
	return true;
}

bool NuclearSegmentation::RunClustering(bool getResultImg)
{
	if(!NucleusSeg || lastRunStep < 2)
	{
		errorMessage = "No Seeds";
		return false;
	}
	NucleusSeg->runClustering();
	lastRunStep = 3;

	this->GetParameters();				//Get the parameters that were used from the module!!

	if(getResultImg)
		return this->GetResultImage();
	return true;
}

bool NuclearSegmentation::Finalize()
{
	if(!NucleusSeg || lastRunStep < 3)
	{
		errorMessage = "No Initial Clustering";
		return false;
	}
	NucleusSeg->runAlphaExpansion();
	lastRunStep = 4;
	return this->GetResultImage();
}

// This function segments nuclei of all time point images in the ftk input image
bool NuclearSegmentation::SegmentAllTimes(bool finalize)
{
	if(!dataImage)
	{
		errorMessage = "No data loaded";
		return false;
	}
	if(labelImage) labelImage = 0;
	labelImage = ftk::Image::New();


	const Image::Info *info = dataImage->GetImageInfo();
	int numStacks = info->numZSlices;			//z-direction
	int numRows = info->numRows;				//y-direction
	int numColumns = info->numColumns; 			//x-direction
	int numTSlices = info->numTSlices; 			//t-direction

	// Set Yousef Nuclear Segmentation:
	if(NucleusSeg) delete NucleusSeg;
	NucleusSeg = new yousef_nucleus_seg();

	if(myParameters.size() == 0)
	{
		NucleusSeg->readParametersFromFile("");		//Will use automatic parameter detection	
	}
	else
	{
		int params[12];
		ConvertParameters(params);
		NucleusSeg->setParams(params); //Will use the parameters that have been set and defauls for the rest
	}
	// Run the segmentation
	for (int t = 0; t<numTSlices; ++t)
	{
		unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(t,channelNumber,0);	//Expects grayscale image 
		NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
		NucleusSeg->runBinarization();
		NucleusSeg->runSeedDetection();
		NucleusSeg->runClustering();
		lastRunStep = 3;
		if(finalize)
		{
			NucleusSeg->runAlphaExpansion();
			lastRunStep = 4;
		}

		std::string filename = dataImage->GetTimeChannelFilenames().at(t).at(channelNumber);
		vector<int> size = NucleusSeg->getImageSize();
		unsigned short *outdptr = NULL;
		if(lastRunStep==3)
			outdptr = NucleusSeg->getClustImage();
		else if(lastRunStep==4)
			outdptr = NucleusSeg->getSegImage();
		if(outdptr)
			labelImage->AppendImageFromData3D(outdptr,itk::ImageIOBase::USHORT,sizeof(unsigned short), size[2], size[1], size[0],filename,true);
		else
		{
			errorMessage = "Error retrieving data pointer";
			return true;
		}
	}
	std::vector< std::vector <std::string> > tmp_filenames;
	for(int i = 0; i< numTSlices; ++i)
	{
		std::vector <std::string> tmp_file;
		std::string name = ftk::GetFilePath(dataImage->GetTimeChannelFilenames().at(i).at(channelNumber))+"\\labeled_"+ftk::GetFilenameFromFullPath(dataImage->GetTimeChannelFilenames().at(i).at(channelNumber));
		tmp_file.push_back(name);
		std::cout<<name<<std::endl;
		tmp_filenames.push_back(tmp_file);
	}
	labelImage->SetTimeChannelFilenames(tmp_filenames);


	this->GetParameters();				//Get the parameters that were used from the module!!
	return true;
		
}




void NuclearSegmentation::GetParameters()
{
	Parameter p1;
	p1.name = "high_sensitivity";
	p1.value = NucleusSeg->getShift();
	this->myParameters.push_back(p1);
	Parameter p2;
	p2.name = "adaptive_binarization";
	p2.value = NucleusSeg->isAdaptiveBinEnabled();
	this->myParameters.push_back(p2);
	Parameter p3;
	p3.name = "LoG_size";
	p3.value = NucleusSeg->getSigma();
	this->myParameters.push_back(p3);
	Parameter p4;
	p4.name = "min_scale";
	p4.value = NucleusSeg->getScaleMin();
	this->myParameters.push_back(p4);
	Parameter p5;
	p5.name = "max_scale";
	p5.value = NucleusSeg->getScaleMax();
	this->myParameters.push_back(p5);
	Parameter p6;
	p6.name = "xy_clustering_res";
	p6.value = NucleusSeg->getRegionXY();
	this->myParameters.push_back(p6);
	Parameter p7;
	p7.name = "z_clustering_res";
	p7.value = NucleusSeg->getRegionZ();
	this->myParameters.push_back(p7);
	Parameter p8;
	p8.name = "finalize_segmentation";
	p8.value = NucleusSeg->isSegmentationFinEnabled();
	this->myParameters.push_back(p8);
	Parameter p9;
	p9.name = "sampling_ratio_XY_to_Z";
	p9.value = NucleusSeg->getSamplingRatio();
	this->myParameters.push_back(p9);
	Parameter p10;
	p10.name = "Use_Distance_Map";
	p10.value = NucleusSeg->isUseDapEnabled();
	this->myParameters.push_back(p10);
	Parameter p11;
	p11.name = "refinement_range";
	p11.value = NucleusSeg->getRefineRange();
	this->myParameters.push_back(p11);
	Parameter p12;
	p12.name = "min_object_size";
	p12.value = NucleusSeg->getMinObjSize();
	this->myParameters.push_back(p12);
}

bool NuclearSegmentation::GetResultImage()
{
	if(!NucleusSeg)
	{
		errorMessage = "Nothing To Get";
		return false;
	}

	vector<int> size = NucleusSeg->getImageSize();
	unsigned short *dptr = NULL;

	switch(lastRunStep)
	{
	case 0:
		errorMessage = "Nothing To Get";
		return false;
		break;
	case 1:		
		dptr = NucleusSeg->getBinImage();
		break;
	case 2:	//Seeds:
		dptr = NucleusSeg->getSeedImage();
		break;
	case 3:
		dptr = NucleusSeg->getClustImage();
		break;
	case 4:
		dptr = NucleusSeg->getSegImage();
		break;
	}

	if(dptr)
	{
		std::vector<unsigned char> color;
		color.assign(3,255);
		if(labelImage) labelImage = 0;
		labelImage = ftk::Image::New();

		if(lastRunStep == 2)
		{
			Cleandptr(dptr,size); // Temporarily deletes the seeds in the background from dptr
			labelImage->AppendChannelFromData3D(dptr, itk::ImageIOBase::USHORT, sizeof(unsigned short), size[2], size[1], size[0], "nuc", color, true);		
			Restoredptr(dptr); // Adds the seeds to dptr which were deleted in Cleandptr
		}
		else
		{
			labelImage->AppendChannelFromData3D(dptr, itk::ImageIOBase::USHORT, sizeof(unsigned short), size[2], size[1], size[0], "nuc", color, true);
		}
	}
	else
	{
		errorMessage = "Error retrieving data pointer";
		return false;
	}
	return true;
}

bool NuclearSegmentation::ComputeAllGeometries(void)
{
	if(!labelImage)
	{
		errorMessage = "No label Image";
		return false;
	}
	if(!dataImage)
	{
		errorMessage = "No data Image";
		return false;
	}

	//Compute region features:
	typedef ftk::LabelImageToFeatures< IntrinsicFeatureCalculator::IPixelT,  IntrinsicFeatureCalculator::LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IntrinsicFeatureCalculator::IPixelT>(0,channelNumber), labelImage->GetItkPtr<IntrinsicFeatureCalculator::LPixelT>(0,0) );
	labFilter->SetLevel(3);
	labFilter->Update();

	bBoxMap.clear();
	centerMap.clear();

	//Now populate centroid and bounding box maps::
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();

		for (int i=0; i<(int)labels.size(); ++i)
		{
			FeatureCalcType::LabelPixelType id = labels.at(i);
			if(id == 0) continue;

			ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);

			Object::Point c;
			c.x = (int)features->Centroid[0];
			c.y = (int)features->Centroid[1];
			c.z = (int)features->Centroid[2];
			c.t = 0;

			Object::Box b;
			b.min.x = (int)features->BoundingBox[0];
			b.max.x = (int)features->BoundingBox[1];
			b.min.y = (int)features->BoundingBox[2];
			b.max.y = (int)features->BoundingBox[3];
			if( labelImage->GetImageInfo()->numZSlices > 2 )
			{
				b.min.z = (int)features->BoundingBox[4];
				b.max.z = (int)features->BoundingBox[5];
			}
			else
			{
				b.min.z = 0;
				b.max.z = 0;
			}
			b.min.t = 0;
			b.max.t = 0;

			bBoxMap[(int)id] = b;
			centerMap[(int)id] = c;

		}
	return true;
}


// For Time Series
bool NuclearSegmentation::ComputeAllGeometries(int ntimes)
{
	if(!labelImage)
	{
		errorMessage = "No label Image";
		return false;
	}
	if(!dataImage)
	{
		errorMessage = "No data Image";
		return false;
	}

	//Compute region features:
	std::vector<ftk::IntrinsicFeatures> featureVector;
	bBoxMap.clear();
	centerMap.clear();
	featureVector.clear();

	typedef ftk::LabelImageToFeatures< IntrinsicFeatureCalculator::IPixelT,  IntrinsicFeatureCalculator::LPixelT, 3 > FeatureCalcType;
	const Image::Info * imInfo = dataImage->GetImageInfo();

	

	for (int t=0; t< imInfo->numTSlices; t++)
	{
		FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
		labFilter->SetImageInputs( dataImage->GetItkPtr<IntrinsicFeatureCalculator::IPixelT>(t,channelNumber), labelImage->GetItkPtr<IntrinsicFeatureCalculator::LPixelT>(t,0) );
		labFilter->ComputeTexturesOn();
		//labFilter->ComputeHistogramOn();
		labFilter->SetLevel(3);
		labFilter->Update();
		std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();

			for (int i=0; i<(int)labels.size(); ++i)
			{
				FeatureCalcType::LabelPixelType id = labels.at(i);
				if(id == 0) continue;

				ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
				featureVector.push_back(*features);
				featureVector.back().num = id;
				featureVector.back().tag = 0;
				featureVector.back().time = t;
//				if(featureVector.back().spacing.size()==0)
//				 {
//					featureVector.back().spacing.push_back(1.0);
//					featureVector.back().spacing.push_back(1.0);
//					featureVector.back().spacing.push_back(5.0);
//				 }

				Object::Point c;
				c.x = (int)features->Centroid[0];
				c.y = (int)features->Centroid[1];
				c.z = (int)features->Centroid[2];
				c.t = t;

				Object::Box b;
				b.min.x = (int)features->BoundingBox[0];
				b.max.x = (int)features->BoundingBox[1];
				b.min.y = (int)features->BoundingBox[2];
				b.max.y = (int)features->BoundingBox[3];
				b.min.z = (int)features->BoundingBox[4];
				b.max.z = (int)features->BoundingBox[5];
				b.min.t = t;
				b.max.t = t+1;

				bBoxMap[(int)id] = b;
				centerMap[(int)id] = c;
		}

		centerMap4DImage.push_back(centerMap);
		bBoxMap4DImage.push_back(bBoxMap);
		featureVector4DImage.push_back(featureVector);
		
		// We can easily update tables compared to featureVectors 
		// when tables are loaded from files. Always better to refer to
		// vtktable than feature vector from other classes
#ifdef USE_TRACKING
		if(!nucsegTrackFeatures.empty())
			this->setCurrentTrackFeatures(t);
#endif
		table4DImage.push_back(featureVectorTovtkTable(featureVector));

		bBoxMap.clear();
		centerMap.clear();
		featureVector.clear();
	}
		centerMap = centerMap4DImage.at(currentTime);
		bBoxMap = bBoxMap4DImage.at(currentTime);
	//Create Mega Table : Concatenated table of all the tables
	createMegaTable();

	return true;
}
#ifdef USE_TRACKING
//This function sets the the current trackfeatures data for table
void NuclearSegmentation::setCurrentTrackFeatures(int time)
{
	currentTrackFeatures.clear();
	currentTrackFeatures = nucsegTrackFeatures.at(time);
}
#endif

// Loads a bunch of VTK tables and creates a mega table
// CONDITION: All the tables should have the same number of columns
void NuclearSegmentation::createMegaTable()
{	
	megaTable = vtkSmartPointer<vtkTable>::New();
	megaTable->Initialize();
	int cols = table4DImage.at(0)->GetNumberOfColumns();
	
	// Need to add columns first...
	// Just the way vtkTable works
    for ( unsigned int i = 0; i < cols; i++ )
    {
		vtkSmartPointer<vtkVariantArray> col = vtkSmartPointer<vtkVariantArray>::New();
		col->SetName(table4DImage.at(0)->GetColumnName(i));
		megaTable->AddColumn ( col );
    }

	for(int i=0;i<table4DImage.size();++i)
	{
		if(table4DImage.at(i)->GetNumberOfColumns()== cols)
			megaTable = ftk::AppendTables(megaTable,table4DImage.at(i));
		else
		{
			megaTable->Delete();
			return;
		}
	}

}


// Adds time to the megatable...
// This is used for correspondence when multiple images are loaded 
// and there are same ids present in multiple images
// We differentiate using the time value. 
// used in NucleusEditor::startActiveLearningMulti()
void NuclearSegmentation::AddTimeToMegaTable()
{
	vtkSmartPointer<vtkVariantArray> col = vtkSmartPointer<vtkVariantArray>::New();
	col->SetName("time");
	col->SetNumberOfValues(megaTable->GetNumberOfRows() );
	megaTable->AddColumn( col );

	int counter = 0;

	for(int i=0;i<table4DImage.size();++i)
	{
	  for(int j=0; j < table4DImage.at(i)->GetNumberOfRows() ;++j)
		{
			megaTable->SetValueByName(counter,"time",i);
			counter++;
		}
	}

}


void NuclearSegmentation::SetCurrentbBox(std::map<int, ftk::Object::Box> currentbBoxMap)
{
	bBoxMap = currentbBoxMap;
}


vtkSmartPointer<vtkTable> NuclearSegmentation::featureVectorTovtkTable(std::vector<ftk::IntrinsicFeatures> featurevector)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

	//Init the table (headers):
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	table->AddColumn(column);
	std::string fPrefix = "";
	for (int i=0; i < ftk::IntrinsicFeatures::N; ++i)
	{
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (fPrefix+ftk::IntrinsicFeatures::Info[i].name).c_str() );
			table->AddColumn(column);
	}

#ifdef USE_TRACKING
	//Add Track Point Features:
	if(!nucsegTrackFeatures.empty())
	{
		for (int i=0; i < ftk::TrackPointFeatures::M; ++i)
		{
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (fPrefix+ftk::TrackPointFeatures::Info[i].name).c_str() );
			table->AddColumn(column);
		}
	}
	// Add Time Feautes:
	//if(!nucsegTimeFeatures.empty())
	//{
	//	for (int i=0; i < ftk::TrackFeatures::NF; ++i)
	//	{
	//		if(!nucsegTimeFeatures.at(i).intrinsic_features.empty())
	//		{
	//			column = vtkSmartPointer<vtkDoubleArray>::New();
	//			column->SetName( (fPrefix+ftk::TrackFeatures::TimeInfo[i].name).c_str() );
	//			table->AddColumn(column);
	//		}
	//	}
	//}

#endif
	//Now populate the table:
	for (int i=0; i<(int)featurevector.size(); ++i)
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue(featurevector.at(i).num);
		for (int j=0; j<ftk::IntrinsicFeatures::N; ++j)
		{
			row->InsertNextValue( vtkVariant(featurevector.at(i).ScalarFeatures[j]) );
		}
#ifdef USE_TRACKING
		// Add Track Point Features:
		if(!nucsegTrackFeatures.empty())
		{
			for (int k=0; k<(int)currentTrackFeatures.size(); ++k)
			{
				if(currentTrackFeatures.at(k).num == featurevector.at(i).num)
				{
					for (int j=0; j<ftk::TrackPointFeatures::M; ++j)
					{
						row->InsertNextValue( vtkVariant(currentTrackFeatures.at(k).scalars[j]) );
					}
				}
			}
		}
		// Add Time Features:
		//if(!nucsegTimeFeatures.empty())
		//{
		//	for (int k=0; k<(int)nucsegTimeFeatures.size(); ++k)			// iterate through tracks
		//	{
		//		if (!nucsegTimeFeatures.at(k).intrinsic_features.empty())
		//		{
		//			if(nucsegTimeFeatures.at(k).intrinsic_features.at(0).num == featurevector.at(i).num)
		//			{
		//				for (int j=0; j<ftk::TrackFeatures::NF; ++j)
		//				{
		//					row->InsertNextValue( vtkVariant(nucsegTimeFeatures.at(i).scalars[j]) );
		//				}
		//			}
		//		}
		//	}
		//}
#endif
		table->InsertNextRow(row);
	}
	return table;
}

#ifdef USE_TRACKING
void NuclearSegmentation::SetTrackFeatures(std::vector<std::vector<ftk::TrackPointFeatures> > trackfeatures)
{
	nucsegTrackFeatures.clear();
	nucsegTrackFeatures = trackfeatures;
}
void NuclearSegmentation::SetTimeFeatures(std::vector<ftk::TrackFeatures> timefeatures)
{
	nucsegTimeFeatures.clear();
	nucsegTimeFeatures = timefeatures;
}
#endif
void NuclearSegmentation::updatetable4DImage(std::vector< vtkSmartPointer<vtkTable> > tableOfFeatures)
{
	table4DImage.clear();
	table4DImage = tableOfFeatures;
	createMegaTable();
}


bool NuclearSegmentation::LoadLabelImage(std::string fname)
{
	ftk::Image::Pointer tmpImg = ftk::Image::New();
	if(!tmpImg->LoadFile(fname))
	{
		errorMessage = "Label Image failed to load";
		return false;
	}
	return SetLabelImage(tmpImg, fname);
}

bool NuclearSegmentation::SetLabelImage(ftk::Image::Pointer labImg, std::string fname)
{
	if(dataImage)
	{
		if(labImg->Size() != this->dataImage->Size())
		{
			this->labelImage = labImg;
			return false;
		}
	}
	this->labelFilename = fname;
	this->labelImage = labImg;
	EditsNotSaved = false;
	return true;
}

bool NuclearSegmentation::SaveLabelImage(std::string fname)
{
	if(!labelImage)
	{
		errorMessage = "Nothing To Save";
		return false;
	}

	if(fname.size() == 0)
	{
		size_t pos = dataFilename.find_last_of(".");
		std::string base = dataFilename.substr(0,pos);
		std::string tag = "_label";
		std::string ext = "tif";
		labelFilename = base + tag + "." + ext;
	}

	size_t pos = labelFilename.find_last_of(".");
	std::string base = labelFilename.substr(0,pos);
	std::string ext = labelFilename.substr(pos+1);

	if(!labelImage->SaveChannelAs(0, base, ext ))
		return false;

	EditsNotSaved = false;
	return true;
}

//The .dat file is an old result from the idl-based segmentation
//This function loads from this format for changing to the new format
bool NuclearSegmentation::LoadFromDAT(std::string fname)
{
	if(!dataImage)
	{
		errorMessage = "Must set/load data image first";
		return false;
	}
	if( !FileExists(fname) )
	{
		errorMessage = "Could not find .dat file";
		return false;
	}

	ReleaseSegMemory();
	NucleusSeg = new yousef_nucleus_seg();
	const Image::Info *info = dataImage->GetImageInfo();
	int numStacks = info->numZSlices;
	int numRows = info->numRows;				//y-direction
	int numColumns = info->numColumns; 			//x-direction

	//We assume that the image is unsigned char, but just in case it isn't we make it so:
	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,channelNumber,0);		//Expects grayscale image	
	NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
	NucleusSeg->readFromIDLFormat(fname);
	this->lastRunStep = 4;
	this->GetResultImage();
	ReleaseSegMemory();
	EditsNotSaved = true;
	return true;
}

/*
bool NuclearSegmentation::SaveLabelByClass()
{
	if(!labelImage)
	{
		errorMessage = "Label Image has not be loaded";
		return false;
	}

	//Cast the label Image & Get ITK Pointer
	typedef unsigned short PixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	labelImage->Cast<PixelType>();
	ImageType::Pointer img = labelImage->GetItkPtr<PixelType>(0,0);

	//Find the possible classes:
	std::vector<int> classes;
	for(int i=0; i<(int)myObjects.size(); ++i)
	{
		int c = (int)myObjects.at(i).GetClass();
		bool found = false;
		for(int l=0; l<(int)classes.size(); ++l)
		{
			if( classes.at(l) == c )	//class in the set
			{
				found = true;
				break;
			}
		}
		if(!found)						//class wasn't found
			classes.push_back(c);		//so add it
	}
	int numClasses = (int)classes.size();
	std::vector< std::set<int> > objClass(numClasses);
	for(int i=0; i<(int)myObjects.size(); ++i)
	{
		int c = (int)myObjects.at(i).GetClass();
		int id = (int)myObjects.at(i).GetId();
		int p = 0;
		for(int j=0; j<numClasses; ++j)
		{
			if(c == classes.at(j))
				break;
			++p;
		}
		if(p < numClasses)
			objClass.at(p).insert(id);
	}
	//Create an image for each class:
	std::vector<ImageType::Pointer> outImgs;
	for(int i=0; i<numClasses; ++i)
	{
		ImageType::Pointer tmp = ImageType::New();   
		tmp->SetOrigin( img->GetOrigin() );
		tmp->SetRegions( img->GetLargestPossibleRegion() );
		tmp->SetSpacing( img->GetSpacing() );
		tmp->Allocate();
		tmp->FillBuffer(0);
		tmp->Update();

		outImgs.push_back(tmp);
	}

	//Iterate through Image & populate all of the other images
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
	IteratorType it(img,img->GetRequestedRegion());
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		int id = it.Get();
		for(int j=0; j<numClasses; ++j)
		{
			if( objClass.at(j).find(id) != objClass.at(j).end() )
			{
				outImgs.at(j)->SetPixel(it.GetIndex(), 1); 
			}
		}	
	}

	//Now Write All of the Images to File
	typedef itk::ImageFileWriter<ImageType> WriterType;
	for(int i=0; i<numClasses; ++i)
	{
		WriterType::Pointer writer = WriterType::New();
		size_t pos = dataFilename.find_last_of(".");
		std::string base = dataFilename.substr(0,pos);
		std::string ext = dataFilename.substr(pos);
		writer->SetFileName( base + "_class" + NumToString(classes.at(i)) + ext );
		writer->SetInput( outImgs.at(i) );
    
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << excp << std::endl;
			writer = 0;
			errorMessage = "Problem saving file to disk";
			return false;
		}
		
		writer = 0;
	}

	editsNotSaved = false;
	return true;
}
*/

//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
// EDITING UTILITES AND EDITING FUNCTIONS
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
// Find the max ID in the centroid map
//**********************************************************************************************************
long int NuclearSegmentation::maxID(void)
{
	long int max = -1;
	
	if (!centerMap4DImage.empty())
		centerMap = centerMap4DImage.at(currentTime);
	
	if(centerMap.size() == 0)
		return max;

	std::map<int, ftk::Object::Point>::iterator it = centerMap.end();
	it--;	//decrement
	max = (*it).first;
	return max;
}

//**********************************************************************************************************
//This function finds the bounding box of the object with id "fromId",
// and uses that area to change the pixels to "toId"
//**********************************************************************************************************
void NuclearSegmentation::ReassignLabel(int fromId, int toId)
{
	std::vector<int> fIds(0);
	fIds.push_back(fromId);
	ReassignLabels(fIds,toId);
}
//**********************************************************************************************************
//This function finds the bounding box of the objects with ids "fromIds",
// and uses that area to change the pixels to "toId" in 1 pass through the whole region
//**********************************************************************************************************
void NuclearSegmentation::ReassignLabels(vector<int> fromIds, int toId)
{
	int C = labelImage->Size()[3];
	int R = labelImage->Size()[2];
	int Z = labelImage->Size()[1];

	ftk::Object::Box region = ExtremaBox(fromIds);
	if(region.min.x < 0) region.min.x = 0;
	if(region.min.y < 0) region.min.y = 0;
	if(region.min.z < 0) region.min.z = 0;
	if(region.max.x >= C) region.max.x = C-1;
	if(region.max.y >= R) region.max.y = R-1;
	if(region.max.z >= Z) region.max.z = Z-1;

	if(labelImage->GetImageInfo()->numZSlices<2)
	{
		region.min.z = 0;
		region.max.z = 0;
	}
	if(labelImage->GetImageInfo()->numTSlices<2)
	{
		region.min.t = 0;
		region.max.t = 0;
	}

	for(int z = region.min.z; z <= region.max.z; ++z)
	{
		for(int r=region.min.y; r <= region.max.y; ++r)
		{
			for(int c=region.min.x; c <= region.max.x; ++c)
			{
				int pix = (int)labelImage->GetPixel(region.min.t,0,z,r,c); // Fix for channels
				for(int i = 0; i < (int)fromIds.size(); ++i)
				{
					if( pix == fromIds.at(i) )
					{
						labelImage->SetPixel(region.min.t,0,z,r,c,toId);  // Fix for channels
					}
				}
			}
		}
	}
}
void NuclearSegmentation::ReassignLabels(std::vector<int> times, std::vector<int> ids, std::vector<int> new_ids)
{
	int C = labelImage->Size()[3];
	int R = labelImage->Size()[2];
	int Z = labelImage->Size()[1];
	for(int i=0 ; i<times.size();++i) // maybe later add check of vector sizes  
	{
		int time = times.at(i);
		int id = ids.at(i);
		int new_id = new_ids.at(i);

		ftk::Object::Box region = bBoxMap4DImage.at(time)[id];
		if(region.min.x < 0) region.min.x = 0;
		if(region.min.y < 0) region.min.y = 0;
		if(region.min.z < 0) region.min.z = 0;
		if(region.max.x >= C) region.max.x = C-1;
		if(region.max.y >= R) region.max.y = R-1;
		if(region.max.z >= Z) region.max.z = Z-1;
		for(int z = region.min.z; z <= region.max.z; ++z)
		{
			for(int r=region.min.y; r <= region.max.y; ++r)
			{
				for(int c=region.min.x; c <= region.max.x; ++c)
				{
					int pix = (int)labelImage->GetPixel(time,0,z,r,c); // Fix for channels
					if( pix == id )
						labelImage->SetPixel(time,0,z,r,c,new_id);  // Fix for channels
				}
			}
		}

		// Update bBoxMap and cMap ids:
		bBoxMap4DImage.at(time).erase( id );
		bBoxMap4DImage.at(time)[new_id] = region;
		ftk::Object::Point point = centerMap4DImage.at(time)[id];
		centerMap4DImage.at(time).erase( id );		
		centerMap4DImage.at(time)[new_id] = point;
		// Update Table:
		for(int row = 0; row<table4DImage.at(time)->GetNumberOfRows(); ++row)
		{
			if(table4DImage.at(time)->GetValue(row,0) == id)
			{
				table4DImage.at(time)->SetValue(row,0,new_id);
				break;
			}
		}

	}
	centerMap = centerMap4DImage.at(currentTime);
	bBoxMap = bBoxMap4DImage.at(currentTime);


}

std::vector<int> NuclearSegmentation::GetNeighbors(int id)
{
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return std::vector<int>(0);
	}
	
	std::vector<int> in_id;
	in_id.push_back(id);
	ftk::Object::Box region = GrowBox( ExtremaBox(in_id), 1 ); //Get the extrema and grow it by 1.

	IntrinsicFeatureCalculator::LPixelT regionIndex[3];
	IntrinsicFeatureCalculator::LPixelT regionSize[3];
	regionIndex[0] = region.min.x;
	regionIndex[1] = region.min.y;
	regionIndex[2] = region.min.z;
	regionSize[0] = region.max.x - region.min.x + 1;
	regionSize[1] = region.max.y - region.min.y + 1;
	regionSize[2] = region.max.z - region.min.z + 1;

	typedef ftk::LabelImageToFeatures< IntrinsicFeatureCalculator::IPixelT,  IntrinsicFeatureCalculator::LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	itk::Image<IntrinsicFeatureCalculator::IPixelT, 3>::Pointer chImg = dataImage->GetItkPtr<IntrinsicFeatureCalculator::IPixelT>(currentTime,channelNumber);
	itk::Image<IntrinsicFeatureCalculator::LPixelT, 3>::Pointer lbImg = labelImage->GetItkPtr<IntrinsicFeatureCalculator::LPixelT>(currentTime,0);
//	labFilter->SetImageInputs( dataImage->GetItkPtr<IntrinsicFeatureCalculator::IPixelT>(0,channelNumber), labelImage->GetItkPtr<IntrinsicFeatureCalculator::LPixelT>(0,0), regionIndex, regionSize );
	labFilter->SetImageInputs( chImg, lbImg, regionIndex, regionSize );
	labFilter->SetLevel(3);	//Needed for neighbor information
	labFilter->Update();

	std::vector<IntrinsicFeatureCalculator::LPixelT> nbs = labFilter->GetContactNeighbors(id);
	std::vector<int> retVector;
	for(int i=0; i<(int)nbs.size(); ++i)
		retVector.push_back(nbs.at(i));
	return retVector;
}
//**********************************************************************************************************
//**********************************************************************************************************
ftk::Object::Box NuclearSegmentation::GrowBox(ftk::Object::Box b, int s)
{
	ftk::Object::Box nb;
	nb.min.x = b.min.x - s;
	nb.min.y = b.min.y - s;
	nb.min.z = b.min.z - s;
	nb.max.x = b.max.x + s;
	nb.max.y = b.max.y + s;
	nb.max.z = b.max.z + s;

	if(!labelImage)		//If no label image return the box
		return nb;

	//If there is a label image check to be sure the new box fits in it.
	if(nb.min.x < 0)
		nb.min.x = 0;
	if(nb.min.y < 0)
		nb.min.y = 0;
	if(nb.min.z < 0)
		nb.min.z = 0;

	if(nb.max.x >= labelImage->GetImageInfo()->numColumns)
		nb.max.x = labelImage->GetImageInfo()->numColumns - 1;
	if(nb.max.y >= labelImage->GetImageInfo()->numRows)
		nb.max.y = labelImage->GetImageInfo()->numRows - 1;
	if(nb.max.z >= labelImage->GetImageInfo()->numZSlices)
		nb.max.z = labelImage->GetImageInfo()->numZSlices - 1;

	return nb;

}

ftk::Object::Box NuclearSegmentation::ExtremaBox(std::vector<int> ids)
{
	ftk::Object::Box extreme;
	if(bBoxMap.size() == 0)
	{
		extreme.min.x=0;
		extreme.max.x=0;
		extreme.min.y=0;
		extreme.max.y=0;
		extreme.min.z=0;
		extreme.max.z=0;
		extreme.min.t=0;
		extreme.max.t=0;
		return extreme;	//Will return bad
	}

	extreme = bBoxMap[ids.at(0)];
	for(int i=1; i<(int)ids.size(); ++i)
	{
		ftk::Object::Box test = bBoxMap[ids.at(i)];

		extreme.min.x = min(extreme.min.x, test.min.x);
		extreme.min.y = min(extreme.min.y, test.min.y);
		extreme.min.z = min(extreme.min.z, test.min.z);
		extreme.min.t = min(extreme.min.t, test.min.t);
		extreme.max.x = max(extreme.max.x, test.max.x);
		extreme.max.y = max(extreme.max.y, test.max.y);
		extreme.max.z = max(extreme.max.z, test.max.z);
		extreme.max.t = max(extreme.max.t, test.max.t);
	}
	return extreme;
}

//**********************************************************************************************************
//**********************************************************************************************************
// EDITING FUNCTIONS:
//**********************************************************************************************************
std::vector< int > NuclearSegmentation::Split(ftk::Object::Point P1, ftk::Object::Point P2, vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> NucAdjTable)
{
	std::vector<int> ret_ids;

	//if no label (segmentation) image return
	if(!labelImage)
	{
		errorMessage = "label image doesn't exist";	
		return ret_ids;
	}

	//Check if the two points inside the same cell
	int id1 = (int)labelImage->GetPixel(0,0,P1.z,P1.y,P1.x);
	int id2 = (int)labelImage->GetPixel(0,0,P2.z,P2.y,P2.x);
	if( id1!=id2 || id1==0 || id2==0)
	{		
		errorMessage = "points are not within the same cell";
		return ret_ids;
	}
		
	int objID = id1;		//The ID of the object I am splitting!!

	//Remove the corresponding rows from the Nuclear Adjacency table
	if(NucAdjTable)
	{
		for(int row=0; row<(int)NucAdjTable->GetNumberOfRows(); ++row)
		{
			if((NucAdjTable->GetValue(row,0).ToInt()==objID)||(NucAdjTable->GetValue(row,1).ToInt()==objID))
			{
				NucAdjTable->RemoveRow(row);
				--row;
			}
		}
	}
	
	//Update the segmentation image
	//Now get the bounding box around the object
	std::vector <int> ids;
	ids.push_back(id1);
	ftk::Object::Box region = ExtremaBox(ids);

	//size of the bounding box:
	std::vector <int> sz;
	sz.push_back(region.max.x - region.min.x + 1);
	sz.push_back(region.max.y - region.min.y + 1);
	sz.push_back(region.max.z - region.min.z + 1);
	
	//get the indexes of the two seeds with respect to the beginning of the bounding box:
	int ind1 = ((P1.z- region.min.z)*sz[0]*sz[1]) + ((P1.y - region.min.y)*sz[0]) + (P1.x - region.min.x);
	int ind2 = ((P2.z- region.min.z)*sz[0]*sz[1]) + ((P2.y - region.min.y)*sz[0]) + (P2.x - region.min.x);

	//create two new itk images with the same size as the bounding box	 	 
	typedef float InputPixelType;
	typedef itk::Image< InputPixelType, 3 > InputImageType;
	InputImageType::Pointer sub_im1 = InputImageType::New();
	InputImageType::Pointer sub_im2 = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0.0; 
    origin[1] = 0.0;    
	origin[2] = 0.0;    
    sub_im1->SetOrigin( origin );
	sub_im2->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = sz[0];  // size along X
    size[1]  = sz[1];  // size along Y
	size[2]  = sz[2];  // size along Z
  
    InputImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
    sub_im1->SetRegions( rgn );
    sub_im1->Allocate();
    sub_im1->FillBuffer(0.0);
	sub_im1->Update();	
	sub_im2->SetRegions( rgn );
    sub_im2->Allocate();
    sub_im2->FillBuffer(0.0);
	sub_im2->Update();

	//set all the points in those images to zeros except for the two points corresponding to the two new seeds
	//notice that one seed is set in each image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(sub_im1,sub_im1->GetRequestedRegion());
	IteratorType iterator2(sub_im2,sub_im2->GetRequestedRegion());	
		
	for(int i=0; i<sz[0]*sz[1]*sz[2]; i++)
	{				
		if(i==ind1)
			iterator1.Set(255.0);
		else
			iterator1.Set(0.0);		
		if(i==ind2)
			iterator2.Set(255.0);
		else
			iterator2.Set(0.0);
		
		++iterator1;	
		++iterator2;	
	}
	
	//compute the distance transforms of those binary itk images
	typedef float OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 > OutputImageType;
	typedef itk::DanielssonDistanceMapImageFilter< InputImageType, OutputImageType > DTFilter;
	DTFilter::Pointer dt_obj1= DTFilter::New();
	DTFilter::Pointer dt_obj2= DTFilter::New();
	dt_obj1->UseImageSpacingOn();
	dt_obj1->SetInput(sub_im1) ;
	dt_obj2->UseImageSpacingOn();
	dt_obj2->SetInput(sub_im2) ;

	try
	{
		dt_obj1->Update() ;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error calculating distance transform: " << err << endl;
		return ret_ids;
	}
	try
	{
		dt_obj2->Update() ;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error calculating distance transform: " << err << endl;
		return ret_ids;
	}

	//Now, relabel the cell points into either newID1 or newID2 based on the distances to the seeds
	int max_id = maxID();		
	int newID1 = ++max_id;		
	int newID2 = ++max_id;
	IteratorType iterator3(dt_obj1->GetOutput(),dt_obj1->GetOutput()->GetRequestedRegion());
	IteratorType iterator4(dt_obj2->GetOutput(),dt_obj2->GetOutput()->GetRequestedRegion());	
	for(int k=region.min.z; k<=region.max.z; k++)
	{
		for(int i=region.min.y; i<=region.max.y; i++)
		{			
			for(int j=region.min.x; j<=region.max.x; j++)
			{
				int d1 = (int) fabs(iterator3.Get());	 
				int d2 = (int) fabs(iterator4.Get());	
				++iterator3;
				++iterator4;
				int pix = (int)labelImage->GetPixel(0,0,k,i,j);
				if(pix != objID)
					continue;
				if(d1>d2)
					labelImage->SetPixel(0,0,k,i,j,newID1);
				else
					labelImage->SetPixel(0,0,k,i,j,newID2);
			}
		}
	}	

	std::set<unsigned short> add_ids;
	add_ids.insert((unsigned short)newID1);
	add_ids.insert((unsigned short)newID2);
	//Add features for the two new objects
	this->addObjectsToMaps(add_ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z, table, NucAdjTable);
	this->removeObjectFromMaps(objID, table);
	EditsNotSaved = true;

	ret_ids.push_back(objID);
	ret_ids.push_back(newID1);
	ret_ids.push_back(newID2);
	return ret_ids;
}


std::vector< std::vector<int> > NuclearSegmentation::BatchSplit(std::vector<int> ids, int numObjs, vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> NucAdjTable)
{
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return std::vector< std::vector<int> >(0);
	}

	std::map<int,int> idToIndexMap;
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		idToIndexMap[table->GetValue(row,0).ToInt()] = row;
	}

	itk::Size<3> sz;
	sz[0] = dataImage->GetImageInfo()->numColumns;
	sz[1] = dataImage->GetImageInfo()->numRows;
	sz[2] = dataImage->GetImageInfo()->numZSlices;
	

	typedef ftk::LabelImageToFeatures< unsigned char, unsigned short, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<unsigned char>(0,0), labelImage->GetItkPtr<unsigned short>(0,0) );
	labFilter->SetLevel( 1 );
	labFilter->Update();

	unsigned char * dataArray = dataImage->GetItkPtr<unsigned char>(0,0)->GetBufferPointer();
	unsigned short * labelArray = labelImage->GetItkPtr<unsigned short>(0,0)->GetBufferPointer();

	std::vector< std::vector<int> > groups;
	groups.resize(ids.size());
	for(int i=0; i<ids.size(); ++i)
	{
		groups[i].push_back(ids[i]);
		int volume_ID = table->GetValueByName(idToIndexMap[ids[i]], "volume").ToInt();
		double estimatedScale;
		if(dataImage->GetImageInfo()->numZSlices > 1)
		{
			estimatedScale = pow((double)volume_ID/numObjs*(3/4)*(1/Pi),1/3);
		}
		else
		{
			estimatedScale = sqrt((double)volume_ID/(numObjs*3.1415));
		}

		int	deltaScale = 5;
		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(ids[i]);
		Object::Box b;
		b.min.x = (int)features->BoundingBox[0];
		b.max.x = (int)features->BoundingBox[1];
		b.min.y = (int)features->BoundingBox[2];
		b.max.y = (int)features->BoundingBox[3];
		b.min.z = (int)features->BoundingBox[4];
		b.max.z = (int)features->BoundingBox[5];

		int numRowsSplit = b.max.y - b.min.y + 1;
		int numColumnsSplit = b.max.x - b.min.x + 1;
		int numStacksSplit = b.max.z - b.min.z + 1;
		//double scaleMin = (double)max<int>((estimatedScale - deltaScale)/1.4142, 3);
		//double scaleMax = (double)(estimatedScale + deltaScale)/1.4142;
		double scaleMin = 1;
		double scaleMax = 50;

		unsigned short *binImagePtr = new unsigned short[numStacksSplit*numRowsSplit*numColumnsSplit];
		float *imgPtr = new float[numStacksSplit*numRowsSplit*numColumnsSplit];
		long ctr=0;
		for(int z=0; z<numStacksSplit; ++z)
		{
			for(int y=0; y<numRowsSplit; ++y)
			{
				for(int x=0; x<numColumnsSplit; ++x)
				{
					unsigned short label = labelArray[(sz[0] * sz[1] * (b.min.z+z)) + (sz[0] * (b.min.y+y)) + (b.min.x+x)];
					if(label == ids[i])
					{
						binImagePtr[ctr] = 255;
						imgPtr[ctr] = 255 - (float)dataArray[(sz[0] * sz[1] * (b.min.z+z)) + (sz[0] * (b.min.y+y)) + (b.min.x+x)];
					}
					else 
					{
						binImagePtr[ctr] = 0;
						imgPtr[ctr] = 255;
					}	
					++ctr;
				}
			}
		}

		typedef itk::Image<unsigned char, 3> UCharImageType;
		typedef itk::Image<unsigned short, 3> UShortImageType;

		UCharImageType::Pointer rawImage = UCharImageType::New();
		itk::Size<3> im_size;
		im_size[0] = numColumnsSplit;
		im_size[1] = numRowsSplit;
		im_size[2] = numStacksSplit;
		UCharImageType::IndexType start;
		start[0] =   0;  // first index on X
		start[1] =   0;  // first index on Y    
		start[2] =   0;  // first index on Z  
		UCharImageType::PointType origin;
		origin[0] = 0; 
		origin[1] = 0;    
		origin[2] = 0;    
		rawImage->SetOrigin( origin );
		UCharImageType::RegionType region;
		region.SetSize( im_size );
		region.SetIndex( start );
		rawImage->SetRegions( region );
		rawImage->Allocate();
		rawImage->FillBuffer(0);
		rawImage->Update();
		UCharImageType::PixelType * rawArray = rawImage->GetBufferPointer();

		ctr = 0;
		for(int z=0; z<numStacksSplit; ++z)
		{
			for(int y=0; y<numRowsSplit; ++y)
			{
				for(int x=0; x<numColumnsSplit; ++x)
				{
					rawArray[ctr] = imgPtr[ctr];
					++ctr;
				}
			}
		}


		itk::ImageFileWriter< itk::Image < unsigned char, 3 > >::Pointer writer1 = itk::ImageFileWriter< itk::Image< unsigned char, 3 > >::New();
		writer1->SetInput(rawImage);
		writer1->SetFileName("C:\\raw_window_1.tif");
		writer1->Update();

		UShortImageType::Pointer binImage = UShortImageType::New();
		UShortImageType::IndexType start1;
		start1[0] =   0;  // first index on X
		start1[1] =   0;  // first index on Y    
		start1[2] =   0;  // first index on Z  
		UShortImageType::PointType origin1;
		origin1[0] = 0; 
		origin1[1] = 0;    
		origin1[2] = 0;    
		binImage->SetOrigin( origin1 );
		UShortImageType::RegionType region1;
		region1.SetSize( im_size );
		region1.SetIndex( start1 );
		binImage->SetRegions( region1 );
		binImage->Allocate();
		binImage->FillBuffer(0);
		binImage->Update();
		UShortImageType::PixelType * binArray = binImage->GetBufferPointer();

		ctr = 0;
		for(int z=0; z<numStacksSplit; ++z)
		{
			for(int y=0; y<numRowsSplit; ++y)
			{
				for(int x=0; x<numColumnsSplit; ++x)
				{
					binArray[ctr] = binImagePtr[ctr];
					++ctr;
				}
			}
		}

		itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::Pointer writer2 = itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::New();
		writer2->SetInput(binImage);
		writer2->SetFileName("C:\\bin_window.tif");
		writer2->Update();
		//ucharToFloat(dataImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);

		float* logImagePtr;	// dont need to create outside of the function, memory allocated inside the function
		unsigned short* seedImagePtr; 

		double regionXY = 2;//(double)myParameters[4].value;
		double regionZ = 2;//(double)myParameters[5].value;
		double sampling_ratio_XY_to_Z = 2;//(double)myParameters[7].value;

		int useDistMap = 1;
		int minLoGImg = 10000;
		bool autoParamEstimation = false;

		//int ok = Seeds_Detection_2D( imgPtr, &logImagePtr, &seedImagePtr, numRowsSplit, numColumnsSplit, numStacksSplit, &scaleMin, &scaleMax, &regionXY, &regionZ, sampling_ratio_XY_to_Z, binImagePtr, useDistMap, &minLoGImg, autoParamEstimation);		
		
		int ok = 0;
		if (numStacksSplit == 1)
		{		
			seedImagePtr = new unsigned short[numStacksSplit*numRowsSplit*numColumnsSplit];		
			logImagePtr = new float[numStacksSplit*numRowsSplit*numColumnsSplit];
			for(int a=scaleMin; a<=scaleMax; ++a)
			{
				double aa = (double)a;
				ok = detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRowsSplit, numColumnsSplit, &aa, &aa, &regionXY, binImagePtr, autoParamEstimation );	
				std::stringstream ss; ss << a;

				itk::Image<float, 3>::Pointer logImage = itk::Image<float, 3>::New();
				itk::Image<float, 3>::IndexType start3;
				start3[0] =   0;  // first index on X
				start3[1] =   0;  // first index on Y    
				start3[2] =   0;  // first index on Z  
				itk::Image<float, 3>::PointType origin3;
				origin3[0] = 0; 
				origin3[1] = 0;    
				origin3[2] = 0;    
				logImage->SetOrigin( origin3 );
				itk::Image<float, 3>::RegionType region3;
				region3.SetSize( im_size );
				region3.SetIndex( start3 );
				logImage->SetRegions( region3 );
				logImage->Allocate();
				logImage->FillBuffer(0);
				logImage->Update();
				itk::Image<float, 3>::PixelType * logArray = logImage->GetBufferPointer();
				ctr = 0;
				for(int z=0; z<numStacksSplit; ++z)
				{
					for(int y=0; y<numRowsSplit; ++y)
					{
						for(int x=0; x<numColumnsSplit; ++x)
						{
							logArray[ctr] = logImagePtr[ctr];
							++ctr;
						}
					}
				}
				itk::ImageFileWriter< itk::Image < float, 3 > >::Pointer writer3 = itk::ImageFileWriter< itk::Image < float, 3 > >::New();
				writer3->SetInput(logImage);
				writer3->SetFileName("C:\\log_window_" + ss.str() + ".mhd");
				writer3->Update();

				UShortImageType::Pointer seedsImage = UShortImageType::New();
				UShortImageType::IndexType start2;
				start2[0] =   0;  // first index on X
				start2[1] =   0;  // first index on Y    
				start2[2] =   0;  // first index on Z  
				UShortImageType::PointType origin2;
				origin2[0] = 0; 
				origin2[1] = 0;    
				origin2[2] = 0;    
				seedsImage->SetOrigin( origin2 );
				UShortImageType::RegionType region2;
				region2.SetSize( im_size );
				region2.SetIndex( start2 );
				seedsImage->SetRegions( region2 );
				seedsImage->Allocate();
				seedsImage->FillBuffer(0);
				seedsImage->Update();
				UShortImageType::PixelType * seedsArray = seedsImage->GetBufferPointer();

				ctr = 0;
				for(int z=0; z<numStacksSplit; ++z)
				{
					for(int y=0; y<numRowsSplit; ++y)
					{
						for(int x=0; x<numColumnsSplit; ++x)
						{
							seedsArray[ctr] = seedImagePtr[ctr];
							++ctr;
						}
					}
				}
				itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::Pointer writer = itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::New();
				writer->SetInput(seedsImage);
				writer->SetFileName("C:\\seeds_window_" + ss.str() + ".mhd");
				writer->Update();

			}
		}
		else
		{	
			minLoGImg = 10000;
			//ok = Seeds_Detection_3D( imgPtr, &logImagePtr, &seedImagePtr, numRows, numColumns, numStacks, &scaleMin, &scaleMax, &regionXY, &regionZ, getSamplingRatio(), binImagePtr, useDistMap, &minLoGImg, autoParamEstimation );						
		}

		itk::Image<float, 3>::Pointer logImage = itk::Image<float, 3>::New();
		itk::Image<float, 3>::IndexType start3;
		start3[0] =   0;  // first index on X
		start3[1] =   0;  // first index on Y    
		start3[2] =   0;  // first index on Z  
		itk::Image<float, 3>::PointType origin3;
		origin3[0] = 0; 
		origin3[1] = 0;    
		origin3[2] = 0;    
		logImage->SetOrigin( origin3 );
		itk::Image<float, 3>::RegionType region3;
		region3.SetSize( im_size );
		region3.SetIndex( start3 );
		logImage->SetRegions( region3 );
		logImage->Allocate();
		logImage->FillBuffer(0);
		logImage->Update();
		itk::Image<float, 3>::PixelType * logArray = logImage->GetBufferPointer();
		ctr = 0;
		for(int z=0; z<numStacksSplit; ++z)
		{
			for(int y=0; y<numRowsSplit; ++y)
			{
				for(int x=0; x<numColumnsSplit; ++x)
				{
					logArray[ctr] = logImagePtr[ctr];
					++ctr;
				}
			}
		}
		itk::ImageFileWriter < itk::Image < float, 3 > >::Pointer writer3 = itk::ImageFileWriter < itk::Image < float, 3 > >::New();
		writer3->SetInput(logImage);
		writer3->SetFileName("C:\\log_window.mhd");
		writer3->Update();



		UShortImageType::Pointer seedsImage = UShortImageType::New();
		UShortImageType::IndexType start2;
		start2[0] =   0;  // first index on X
		start2[1] =   0;  // first index on Y    
		start2[2] =   0;  // first index on Z  
		UShortImageType::PointType origin2;
		origin2[0] = 0; 
		origin2[1] = 0;    
		origin2[2] = 0;    
		seedsImage->SetOrigin( origin2 );
		UShortImageType::RegionType region2;
		region2.SetSize( im_size );
		region2.SetIndex( start2 );
		seedsImage->SetRegions( region2 );
		seedsImage->Allocate();
		seedsImage->FillBuffer(0);
		seedsImage->Update();
		UShortImageType::PixelType * seedsArray = seedsImage->GetBufferPointer();

		ctr = 0;
		for(int z=0; z<numStacksSplit; ++z)
		{
			for(int y=0; y<numRowsSplit; ++y)
			{
				for(int x=0; x<numColumnsSplit; ++x)
				{
					seedsArray[ctr] = seedImagePtr[ctr];
					++ctr;
				}
			}
		}
		itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::Pointer writer = itk::ImageFileWriter< itk::Image < unsigned short, 3 > >::New();
		writer->SetInput(seedsImage);
		writer->SetFileName("C:\\seeds_window.tif");
		writer->Update();


	}

	return groups;

}


std::vector< int > NuclearSegmentation::SplitAlongZ(int objID, int cutSlice, vtkSmartPointer<vtkTable> table)
{
	std::vector <int> ret_ids;

	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";	
		return ret_ids;
	}

	std::vector<itk::SizeValueType> size = labelImage->Size();
	if(size[1] == 1)	//Only 1 z slice
	{
		errorMessage = "2D image cannot be split along z";
		return ret_ids;
	}
	if(bBoxMap.size() == 0)
	{
		errorMessage = "bounding boxes not known";
		return ret_ids;
	}
	
	//Get the bounding box around the object
	ftk::Object::Box region = bBoxMap[objID];	

	//Now, relabel the cell points into either newID1 or newID2 based on the z-slice
	int max_id = maxID();		
	int newID1 = ++max_id;		
	int newID2 = ++max_id;		
	for(int k=region.min.z; k<=region.max.z; k++)
	{
		for(int i=region.min.y; i<=region.max.y; i++)
		{			
			for(int j=region.min.x; j<=region.max.x; j++)
			{				
				int pix = (int)labelImage->GetPixel(0,0,k,i,j);
				if(pix != objID) continue;
				if(k<cutSlice)
					labelImage->SetPixel(0,0,k,i,j,newID1);
				else
					labelImage->SetPixel(0,0,k,i,j,newID2);
			}
		}
	}	

	std::set<unsigned short> ids;
	ids.insert((unsigned short)newID1);
	ids.insert((unsigned short)newID2);
	//Add features for the two new objects
	this->addObjectsToMaps(ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z, table);
	this->removeObjectFromMaps(objID, table);
	EditsNotSaved = true;

	//return the ids of the two cells resulting from spliting
	ret_ids.push_back(objID);
	ret_ids.push_back(newID1);
	ret_ids.push_back(newID2);
	return ret_ids;
}

std::vector< std::vector<int> > NuclearSegmentation::GroupMerge(std::vector<int> ids, vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> NucAdjTable)
{
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return std::vector< std::vector<int> >(0);
	}

	std::vector< std::vector<int> > groups;

	while(ids.size() != 0)
	{
		//Get neighbors of first object in vector:
		unsigned currentID = ids.at(0);
		std::vector<int> nbs = GetNeighbors(currentID);

		bool added = false;
		//iterate through groups, to see if I belong to a different group (if one of my neighbors is already in a group):
		for(unsigned n=0; n<nbs.size(); ++n)
		{
			int currentN = nbs.at(n);
			for(unsigned i=0; i<groups.size(); ++i)
			{
				for(unsigned j=0; j<groups.at(i).size(); ++j)
				{
					if(groups.at(i).at(j) == currentN) //Found match so add current ID to group
					{
						groups.at(i).push_back(currentID);
						added = true;
						break;
					}
				}
			}
		}

		if(!added)
		{
			std::vector<int> grp;

			//iterate though ids to find out which ones are neighbors and add to group.
			for(unsigned i=0; i<ids.size(); ++i)
			{
				for(unsigned j=0; j<nbs.size(); ++j)
				{
					if(ids.at(i) == nbs.at(j))
					{
						grp.push_back( ids.at(i) );
					}
				}
			}

			//Remove everything in the group from the input ids vector:
			for(unsigned g=0; g<grp.size(); ++g)
			{
				for(unsigned i=0; i<ids.size(); ++i)
				{
					if(grp.at(g) == ids.at(i))
					{
						ids.erase( ids.begin()+i );
						break;
					}
				}
			}

			//If I found neighbors add myself (currentID) to the group:
			if(grp.size() > 0) 
			{
				grp.push_back( currentID );
				groups.push_back( grp );
			}

		}

		//erase currentID from the input ids:
		for(unsigned i=0; i<ids.size(); ++i)
		{
			if(ids.at(i) == currentID)
			{
				ids.erase( ids.begin() + i );
			}
		}
	}

	//Iterate though groups and merge them together
	for(unsigned i=0; i<groups.size(); ++i)
	{
		int newID = Merge(groups.at(i), table, NucAdjTable);
		groups.at(i).push_back(newID);	//Put the new id at the end of the group
	}
	return groups;

}

int NuclearSegmentation::Merge(vector<int> ids, vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> NucAdjTable)
{
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return -1;
	}

	int newID = maxID() + 1;
	if(NucAdjTable)
	{
		for(unsigned j=0; j<ids.size(); ++j)
		{
			int OldID = ids.at(j);
			for(unsigned row=0; row<NucAdjTable->GetNumberOfRows(); ++row)
			{
				for(unsigned col=0; col<NucAdjTable->GetNumberOfRows(); ++col)
				{
					if(NucAdjTable->GetValue(row,col).ToInt() == OldID)
						NucAdjTable->SetValue(row,col,newID);
				}
			}
		}
		for(unsigned row=0; row<NucAdjTable->GetNumberOfRows(); ++row)
		{
			if((NucAdjTable->GetValue(row,0).ToInt()) == (NucAdjTable->GetValue(row,1).ToInt()))
			{
				NucAdjTable->RemoveRow(row);
				--row;
			}
		}
	}

	this->ReassignLabels(ids, newID);					//Assign all old labels to this new label
	ftk::Object::Box region = ExtremaBox(ids);
	this->addObjectToMaps(newID, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z, table);
	for(unsigned i=0; i<ids.size(); ++i)
	{
		this->removeObjectFromMaps(ids.at(i),table);
	}
	EditsNotSaved = true;

	return newID;
}

int NuclearSegmentation::AddObject(int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table)
{
	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return -1;
	}	
	std::vector<itk::SizeValueType> size = labelImage->Size();

	int sz_x = x2-x1+1;
	int sz_y = y2-y1+1;
	if(z1==z2)
	{
		//assume that the sampling ratio is 2
		int dz;
		if(sz_x > sz_y)
			dz = sz_x/4;
		else
			dz = sz_y/4;

		z1 -= dz;
		if(z1<0)
			z1=0;
		z2 += dz;
		if(z2>size[1]-1)
			z2=size[1]-1;
	}
	
	int sz_z = z2-z1+1;
	if(sz_x<1 || sz_y<1 || sz_z<1)
		return 0;

	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,channelNumber,0);		//Expects grayscale image
	unsigned short *lptr = labelImage->GetSlicePtr<unsigned short>(0,0,0);				//Expects grayscale image

	ReleaseSegMemory();	//If I'm in add mode must be done with segmentation!!!
	NucleusSeg = new yousef_nucleus_seg();

	//create std vectors of the points
	//I am doing that because I want to make yousef_seg isolated from ftk
	std::vector<int> p1;
	std::vector<int> p2;	
	p1.push_back(x1);
	p1.push_back(y1);
	p1.push_back(z1);
	p2.push_back(x2);
	p2.push_back(y2);
	p2.push_back(z2);
	
	int newID = 0;
	if(size[1] == 1)
		newID = NucleusSeg->AddObject2D(dptr, lptr, p1,p2,size, (int)maxID());	//these are "static" methods
	else
		newID = NucleusSeg->AddObject(dptr, lptr, p1,p2,size, (int)maxID());

	if(newID == 0) return 0;

	lastRunStep = 4;	//To make sure we retrieve the correct image!!!
	this->GetResultImage();
	this->addObjectToMaps(newID, x1, y1, z1, x2, y2, z2, table);
	EditsNotSaved = true;

	return newID;
}

bool NuclearSegmentation::Delete(std::vector<int> ids, vtkSmartPointer<vtkTable> table)
{
	if(!labelImage) return false;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		ReassignLabel(ids.at(i),0);				//Turn each label in list to zero
		removeObjectFromMaps(ids.at(i),table);
	}
	EditsNotSaved = true;

	return true;
}

bool NuclearSegmentation::Exclude(int l, int r, int t, int b, int z1, int z2, vtkSmartPointer<vtkTable> table)
{
	if(!labelImage) return false;

	//Find the bounds to use for exclusion margin:
	const ftk::Image::Info *info = labelImage->GetImageInfo();
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;
	int zSlices = (*info).numZSlices;

	int min_x = 0 + l;
	int min_y = 0 + t;
	int min_z = 0 + z1;
	int max_x = totalWidth - r - 1;
	int max_y = totalHeight - b - 1;
	int max_z = zSlices - z2 - 1;

	//Go through each object and exclude accordingly.
	std::vector<int> ids;
	std::map<int, ftk::Object::Point>::iterator it;
	for(it=centerMap.begin(); it!=centerMap.end(); ++it)
	{
		int id = (*it).first;
		ftk::Object::Point c = (*it).second;

		if(c.x < min_x || c.x > max_x ||
		   c.y < min_y || c.y > max_y ||
		   c.z < min_z || c.z > max_z)
		{
			ids.push_back(id);
		}
		
	}
	
	for(int i=0; i<(int)ids.size(); ++i)
	{
		ReassignLabel(ids.at(i),0);				//Turn each label in list to zero
		removeObjectFromMaps(ids.at(i),table);			
	}
	EditsNotSaved = true;

	return true;
}

//Used to fill objects with holes
bool NuclearSegmentation::FillObjects(std::vector<int> ids)
{
	if(!labelImage) return false;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		FillAnObject(ids.at(i));				//fill them one by one		
	}
	EditsNotSaved = true;

	return true;
}

//Used to fill a single object with holes
bool NuclearSegmentation::FillAnObject(int objID)
{
	//Start by getting the bounding box around the object of interest
	ftk::Object::Box region = bBoxMap[objID];	
	//Extend the bounding box by 10 from each size
	std::vector<itk::SizeValueType> SZ = labelImage->Size();
	int min_x = region.min.x-10;
	if(min_x<0)
		min_x = 0;
	int min_y = region.min.y-10;
	if(min_y<0)
		min_y = 0;
	int min_z = region.min.z-10;
	if(min_z<0)
		min_z = 0;
	int max_x = region.max.x+10;
	if(max_x > SZ[3]-1)
		max_x = SZ[3]-1;
	int max_y = region.max.y+10;
	if(max_y > SZ[2]-1)
		max_y = SZ[2]-1;
	int max_z = region.max.z+10;
	if(max_z > SZ[1]-1)
		max_z = SZ[1]-1;

	//create an ITK image of the same size as the bounding box	
	typedef unsigned short PixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	ImageType::Pointer img = ImageType::New();
	ImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    img->SetOrigin( origin );

    ImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    ImageType::SizeType  size;
 //   size[0]  = region.max.x-region.min.x+1;  // size along X
 //   size[1]  = region.max.y-region.min.y+1;  // size along Y
	//size[2]  = region.max.z-region.min.z+1;  // size along Z
	size[0] = max_x-min_x+1;
	size[1] = max_y-min_y+1;
	size[2] = max_z-min_z+1;
  
    ImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
    double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = 1; //spacing along z //leave it for 1 as of now

    img->SetRegions( rgn );
	img->SetSpacing(spacing);
    img->Allocate();
    img->FillBuffer(0);
	img->Update();	
		
	//Iterate through Image & fill in with the object of interest
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType it(img,img->GetRequestedRegion());
	/*for(int k=region.min.z; k<=region.max.z; k++)
	{
		for(int i=region.min.y; i<=region.max.y; i++)
		{			
			for(int j=region.min.x; j<=region.max.x; j++)*/
	for(int k=min_z; k<=max_z; k++)
	{
		for(int i=min_y; i<=max_y; i++)
		{			
			for(int j=min_x; j<=max_x; j++)
			{
				int pix = (int)labelImage->GetPixel(0,0,k,i,j);
				if(pix != objID)
					it.Set(0);
				else
					it.Set(255);
				++it;
			}
		}
	}

	//Fill the itk image
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();
	ImageType::SizeType indexRadius;
	indexRadius[0] = 5; // radius along x
	indexRadius[1] = 5; // radius along y
	indexRadius[2] = 5; // radius along z
	filter->SetRadius( indexRadius );
	filter->SetBackgroundValue( 0 );
	filter->SetForegroundValue( 255 );
	filter->SetMajorityThreshold( 2 );
	filter->SetMaximumNumberOfIterations( 5 );
	filter->SetInput( img );
	filter->Update();

	//Copy the filled itk image back to the label image	
	IteratorType it2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());
	/*for(int k=region.min.z; k<=region.max.z; k++)
	{
		for(int i=region.min.y; i<=region.max.y; i++)
		{			
			for(int j=region.min.x; j<=region.max.x; j++)*/
	for(int k=min_z; k<=max_z; k++)
	{
		for(int i=min_y; i<=max_y; i++)
		{			
			for(int j=min_x; j<=max_x; j++)
			{
				int pix = (int)labelImage->GetPixel(0,0,k,i,j);
				int pix2 = (int) it2.Get();
				if(pix ==0 && pix2>0) //just copy pixels that were zeros
					labelImage->SetPixel(0,0,k,i,j,objID);										
				++it2;
			}
		}
	}
	
	return true;
}

//*****************************************************************************
// Applied to initial segmentation, doesn't change table (table doesn't exist
//*****************************************************************************
bool NuclearSegmentation::DeleteInit(ftk::Object::Point P1)
{
	if(!NucleusSeg) return false;

	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";			
		return false;
	}	
	
	bool ids_ok = NucleusSeg->DeleteInit(P1);
	return ids_ok;
}

//This function is applied on the initial segmentation image and updates the LoG response image
std::vector< int > NuclearSegmentation::SplitInit(ftk::Object::Point P1, ftk::Object::Point P2)
{
	if(!NucleusSeg) return std::vector<int>(0);
	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";	
		std::vector <int> ids_err;
		ids_err.push_back(0);
		ids_err.push_back(0);
		return ids_err;
	}	
	
	//Apply the splitting
	std::vector <int> ids_ok = NucleusSeg->SplitInit(P1, P2);
		
	return ids_ok;
}

//this is used when we apply merging on the initial segmentation
ftk::Object::Point NuclearSegmentation::MergeInit(ftk::Object::Point P1, ftk::Object::Point P2, int* new_id)
{
	ftk::Object::Point newSeed;
	newSeed.t = newSeed.x = newSeed.y = newSeed.z = 0;
	if(!NucleusSeg) return newSeed; 
	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";			
		return newSeed;
	}	
	
	//Apply the splitting
	newSeed = NucleusSeg->MergeInit(P1, P2, new_id);
		
	return newSeed;
}
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
void NuclearSegmentation::removeObjectFromMaps(int ID, vtkSmartPointer<vtkTable> table)
{
	if(table)
	{
		for(int row = 0; row<table->GetNumberOfRows(); ++row)
		{
			if(table->GetValue(row,0) == ID)
			{
				table->RemoveRow( row );
				break;
			}
		}
	}

	//if(zernikeTable)
	//{
	//	for(int row = 0; row<zernikeTable->GetNumberOfRows(); ++row)
	//	{
	//		if(zernikeTable->GetValue(row,0) == ID)
	//		{
	//			zernikeTable->RemoveRow( row );
	//			break;
	//		}
	//	}
	//}
	
	centerMap.erase( ID );
	bBoxMap.erase( ID );
	if((!centerMap4DImage.empty())&& (!bBoxMap4DImage.empty()))
	{
		centerMap4DImage.at(currentTime).erase( ID );
		bBoxMap4DImage.at(currentTime).erase( ID );
	}
	
}

//Calculate the features within a specific region of the image for a specific ID, and update the table
bool NuclearSegmentation::addObjectToMaps(int ID, int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table)
{
	std::set<unsigned short> ids;
	ids.insert((unsigned short)ID);
	return addObjectsToMaps(ids,x1,y1,z1,x2,y2,z2,table);
}

bool NuclearSegmentation::addObjectsToMaps(std::set<unsigned short> IDs, int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> NucAdjTable)
{
	IntrinsicFeatureCalculator * calc = new IntrinsicFeatureCalculator();
	ftk::Image::Pointer chImg = ftk::Image::New();
	ftk::Image::Pointer lbImg = ftk::Image::New();
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(2); //use DEEP_COPY mode

	itk::Image<unsigned char, 3>::Pointer dataImageItk = dataImage->GetItkPtr<unsigned char>(currentTime,0, mode);
	itk::Image<unsigned short, 3>::Pointer labelImageItk = labelImage->GetItkPtr<unsigned short>(currentTime,0, mode);

	//chImg->CreateData3DFromItkPtr<unsigned char>(dataImageItk);
	//lbImg->CreateData3DFromItkPtr<unsigned short>(labelImageItk);

	ftk::Image::DataType dataType = dataImage->GetImageInfo()->dataType;
	ftk::Image::DataType labelType = labelImage->GetImageInfo()->dataType;
	unsigned char databpPix = dataImage->GetImageInfo()->bytesPerPix;
	unsigned char labelbpPix = labelImage->GetImageInfo()->bytesPerPix;
	unsigned short cs = dataImage->GetImageInfo()->numColumns;
	unsigned short rs = dataImage->GetImageInfo()->numRows;
	unsigned short zs = dataImage->GetImageInfo()->numZSlices;
	std::string name;
	if (dataImage->GetImageInfo()->numTSlices > 2)
	{
		std::vector< std::vector <std::string> > FileNames = dataImage->GetTimeChannelFilenames();
		name = FileNames.at(currentTime).at(0);
	}
	else
	{
		std::vector< std::string > FileNames = dataImage->GetFilenames();
		name = FileNames.at(0);
	}


	std::vector<unsigned char> color(3,255);

	chImg->AppendChannelFromData3D(dataImageItk->GetBufferPointer(), dataType, databpPix, cs, rs, zs, name, color, true);
	lbImg->AppendChannelFromData3D(labelImageItk->GetBufferPointer(), labelType, labelbpPix, cs, rs, zs, name, color, true);

	calc->SetInputImages(chImg, lbImg);
//	calc->SetInputImages(dataImage, labelImage);
	calc->SetRegion(x1,y1,z1,x2,y2,z2);
	calc->SetIDs(IDs);
	if((!centerMap4DImage.empty())&& (!bBoxMap4DImage.empty()))
	{
		if (centerMap.empty())
			centerMap = centerMap4DImage.at(currentTime);
		if (bBoxMap.empty())
			bBoxMap = bBoxMap4DImage.at(currentTime);
	}
	calc->Update(table, &centerMap, &bBoxMap, NucAdjTable, currentTime);
	//calc->UpdateZernike(zernikeTable);
	delete calc;


	/*
	FeatureCalcType::Pointer labFilter = computeGeometries(x1,y1,z1,x2,y2,z2);
	if(!labFilter) return false;

	//Add to the maps:
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);

		Object::Point c;
		c.x = (int)features->Centroid[0];
		c.y = (int)features->Centroid[1];
		c.z = (int)features->Centroid[2];
		c.t = 0;

		Object::Box b;
		b.min.x = (int)features->BoundingBox[0];
		b.max.x = (int)features->BoundingBox[1];
		b.min.y = (int)features->BoundingBox[2];
		b.max.y = (int)features->BoundingBox[3];
		b.min.z = (int)features->BoundingBox[4];
		b.max.z = (int)features->BoundingBox[5];
		b.min.t = 0;
		b.max.t = 0;

		bBoxMap[(int)id] = b;
		centerMap[(int)id] = c;
	}
	*/
	return true;
}

//Calculate the features within a specific region of the image and return the filter:
//bool NuclearSegmentation::computeFeatures(int x1, int y1, int z1, int x2, int y2, int z2)
/*
FeatureCalcType::Pointer NuclearSegmentation::computeGeometries(int x1, int y1, int z1, int x2, int y2, int z2)
{
	
	if(!dataImage)
	{
		errorMessage = "No Data Image";
		return NULL;
	}
	if(!labelImage)
	{
		errorMessage = "No Label Image";
		return NULL;		
	}

	//Calculate features using feature filter
	typedef itk::Image< IPixelT, 3 > IImageT;
	typedef itk::Image< LPixelT, 3 > LImageT;

	IImageT::Pointer itkIntImg = dataImage->GetItkPtr<IPixelT>(0,0);
	LImageT::Pointer itkLabImg = labelImage->GetItkPtr<LPixelT>(0,0);

	LPixelT index[3];
	index[0] = x1;
	index[1] = y1;
	index[2] = z1;

	LPixelT size[3];
	size[0] = x2 - x1 + 1;
	size[1] = y2 - y1 + 1;
	size[2] = z2 - z1 + 1;

	//Compute features:
	//typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( itkIntImg, itkLabImg, index, size );
	labFilter->SetLevel(1);
	labFilter->Update();
	return labFilter;
}
*/
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
void NuclearSegmentation::Cleandptr(unsigned short* p, vector<int> dim){
	int ctr =0;
	
if(dim.size() ==3) {
	for (int index1=0;index1<dim[2];index1++)
		{
			for(int index2=0;index2<dim[1];index2++)
				{
					for(int index3=0;index3<dim[0];index3++)
						{
						    //if(p[ctr]<0)klkl 
							if(p[ctr]==65535) 
							{
								p[ctr]=0;
								this->negativeseeds.push_back(ctr);									
							}
					ctr++;
						}
				}
		}
	}
else
{
for(int index1=0;index1<dim[1];index1++)
	{
	for(int index2=0;index2<dim[0];index2++)
		{
			if(p[ctr]==65535/*<0*/) 
				{
				p[ctr]=0;
				this->negativeseeds.push_back(ctr);									
				}
				ctr++;
		}
	}
}

}

void NuclearSegmentation::Restoredptr(unsigned short* p)
{
 	for(list<int>::iterator index =this->negativeseeds.begin();index!=this->negativeseeds.end();++index)
		{
	    		p[*index]=65535;//-1;					
		}
	this->negativeseeds.clear();
}

std::vector<Seed> NuclearSegmentation::getSeeds()
{
	if(NucleusSeg)
		return NucleusSeg->getSeedsList();
	else
		return std::vector<Seed>(0);
}
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//Added by Yousef on 04-08-2009
//This function will run graph coloring and will assign different colors for touching objects
//For now, it will just write the list of labels into a text file, but this should be relaxed later
//This function runs graph coloring 
//*******************************************************************************************************
std::vector<std::string> NuclearSegmentation::RunGraphColoring(std::string labelname, std::string filename)
{
	//get the label image (if not already done)
	std::cout<<"Loading Label Image ... ";
	labelImage = ftk::Image::New();
	labelImage->LoadFile(labelname);
	std::cout<<"done!"<<endl;

    int max_lab;
    int** RAG;    
    int* ColorOut;        
	int L, L1, L2, L3, L4, L5, L6, L7;
	std::vector<std::string> colorImages;
	
	int c = labelImage->Size()[3];
	int r = labelImage->Size()[2];
	int z = labelImage->Size()[1];	
	unsigned short* labs_vals = static_cast<unsigned short*> (labelImage->GetDataPtr(0,0));
	
	//get the maximum label
	std::cout<<"image size is "<<r<<"x"<<c<<"x"<<z<<std::endl;	
	max_lab = 0;
	for(int i=0; i<r-1; i++)
    {        		
		for(int j=0; j<c-1; j++)
		{						
			for(int k=0; k<z-1; k++)
			{	
				if((int)labs_vals[(k*r*c)+(j*r)+i]>max_lab)
					max_lab = (int)labs_vals[(k*r*c)+(j*r)+i];
			}
		}
	}
	std::cout<<"The maximum cell label is "<<max_lab<<std::endl;
	  
 //   //Build the region adjacency graph    
    std::cout<<"Building Region Adjacency Graph...";
    RAG = (int **) malloc(max_lab*sizeof(int*));
    for(int i=0; i<max_lab; i++)
    {        
		RAG[i] = (int *) malloc(max_lab*sizeof(int));
        for(int j=0; j<max_lab; j++)
            RAG[i][j] = 0;
    }
    for(int i=0; i<r-1; i++)
    {        
        for(int j=0; j<c-1; j++)
        {	
			for(int k=0; k<z-1; k++)
			{
				L = labs_vals[(k*r*c)+(j*r)+i];
				if( L == 0)
					continue;
				else
				{			
					L1 = labs_vals[(k*r*c)+(j*r)+(i+1)]; 
					L2 = labs_vals[(k*r*c)+((j+1)*r)+i];
					L3 = labs_vals[(k*r*c)+((j+1)*r)+(i+1)];
					L4 = labs_vals[((k+1)*r*c)+(j*r)+i];
					L5 = labs_vals[((k+1)*r*c)+((j+1)*r)+i];
					L6 = labs_vals[((k+1)*r*c)+(j*r)+(i+1)];
					L7 = labs_vals[((k+1)*r*c)+((j+1)*r)+(i+1)];

					if(L!=L1 && L1!=0)
						RAG[L-1][L1-1] = RAG[L1-1][L-1] = 1;
					if(L!=L2 && L2!=0)
						RAG[L-1][L2-1] = RAG[L2-1][L-1] = 1;
					if(L!=L3 && L3!=0)
						RAG[L-1][L3-1] = RAG[L3-1][L-1] = 1;
					if(L!=L4 && L4!=0)
						RAG[L-1][L4-1] = RAG[L4-1][L-1] = 1;
					if(L!=L5 && L5!=0)
						RAG[L-1][L5-1] = RAG[L5-1][L-1] = 1;
					if(L!=L6 && L6!=0)
						RAG[L-1][L6-1] = RAG[L6-1][L-1] = 1;
					if(L!=L7 && L7!=0)
						RAG[L-1][L7-1] = RAG[L7-1][L-1] = 1;
				}
            }                		
        }		
    }    
	std::cout<<"done!"<<endl;

    //copy the RAG into an std vector of vectors
	std::vector<std::vector<int> > MAP;
	MAP.resize(max_lab);
    std::vector<std::vector<int> > MAP2;
	MAP2.resize(max_lab);
	
	ColorOut = (int *) malloc(max_lab*sizeof(int));
	for(int i=0; i<max_lab; i++)
	{	
		ColorOut[i] = 0;
		int isIsolated = 1;
		for(int j=0; j<max_lab; j++)
		{
			if(RAG[i][j]==1)
            {
				MAP[i].push_back(j+1);   
				isIsolated = 0;
            }
		} 
		if(isIsolated==1)
			ColorOut[i] = 1;
		free(RAG[i]);
	}    
	free(RAG);
    
    //start the graph coloring using Sumit's sequential coloring code
    GVC* Gcol = new GVC(); 			
 	Gcol->sequential_coloring(max_lab,  max_lab, ColorOut, MAP );
	int numColors = 0;
	for(int i=0; i<max_lab; i++)
	{
		int c = ColorOut[i]+1;
		if(c>numColors)
			numColors=c;
	}
    std::cout<<"Graph Coloring Done"<<endl;
	//write the resulting colors into a file
	FILE *fp = fopen(filename.c_str(),"w");
	
	if(fp == NULL)
	{
		fprintf(stderr,"can't open %s for writing\n",filename.c_str());
		exit(1);
	}
	for(int i=0; i<max_lab; i++)
	{
		fprintf(fp,"%d\n",ColorOut[i]+1);
		
	}
	fclose(fp);
	//Try this: save the colors into the classes list
	std::vector<int> classes;
	classes.clear();
	for(int i=0; i<numColors; i++)
		classes.push_back(i+1);

	//Cast the label Image & Get ITK Pointer
	typedef unsigned short PixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	labelImage->Cast<PixelType>();
	ImageType::Pointer img = labelImage->GetItkPtr<PixelType>(0,0);

	//Create an image for each class:	
	std::cout<<"Creating an image for each class...";
	std::vector<ImageType::Pointer> outImgs;
	for(int i=0; i<numColors; ++i)
	{
		ImageType::Pointer tmp = ImageType::New();   
		tmp->SetOrigin( img->GetOrigin() );
		tmp->SetRegions( img->GetLargestPossibleRegion() );
		tmp->SetSpacing( img->GetSpacing() );
		tmp->Allocate();
		tmp->FillBuffer(0);
		tmp->Update();

		outImgs.push_back(tmp);
	}

	//create lists of object ids in each class:
	std::vector< std::set<int> > objClass(numColors);
	for(int i=0; i<max_lab; ++i)
	{
		int c = ColorOut[i]+1;
		int id = i+1;
		int p = 0;
		for(int j=0; j<numColors; ++j)
		{
			if(c == classes.at(j))
				break;
			++p;
		}
		if(p < numColors)
			objClass.at(p).insert(id);
	}

	//Iterate through Image & populate all of the other images
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
	IteratorType it(img,img->GetRequestedRegion());
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		int id = it.Get();
		for(int j=0; j<numColors; ++j)
		{
			if( objClass.at(j).find(id) != objClass.at(j).end() )
			{
				outImgs.at(j)->SetPixel(it.GetIndex(), 1); 
			}
		}	
	}

	//Now Write All of the Images to File
	typedef itk::ImageFileWriter<ImageType> WriterType;
	for(int i=0; i<numColors; ++i)
	{
		WriterType::Pointer writer = WriterType::New();
		size_t pos = filename.find_last_of(".");
		std::string base = filename.substr(0,pos);
		std::string ext = filename.substr(pos);
		std::string colorImage;
		writer->SetFileName( base + "_class" + NumToString(classes.at(i)) + ext );
		writer->SetInput( outImgs.at(i) );
		colorImage = writer->GetFileName();
		colorImages.push_back(colorImage);

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << excp << std::endl;
			writer = 0;
			errorMessage = "Problem saving file to disk";
			return std::vector<std::string>(0);
		}
		
		writer = 0;
	}
	std::cout<<"done!"<<std::endl;
	return colorImages;
}
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
} //END NAMESPACE FTK
