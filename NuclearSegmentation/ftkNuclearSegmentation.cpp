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
#include "ftkNuclearSegmentation.h"
#include <itkImageRegionConstIteratorWithIndex.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "yousef_core/graphColLearn_3D/sequential_coloring.cpp"

namespace ftk 
{

//Constructor
NuclearSegmentation::NuclearSegmentation()
{
	NucleusSeg = NULL;
	this->ResetAll();
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
	labelFilename.clear();
	paramFilename.clear();
	featureFilename.clear();
	headerFilename.clear();

	featureTable = NULL;
	dataImage = NULL;
	channelNumber = 0;
	labelImage = NULL;
	editsNotSaved = false;

	bBoxMap.clear();
	centerMap.clear();
	//idToRowMap.clear();
	myParameters.clear();
	myEditRecords.clear();

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

//*************************************************************************************************************
// This function will initialize the pipeline for a new nuclear segmentation.
// It should be called first.  If this class already contains data and the datafilenames do not match, we delete
// all data and prepare for a new segmentation.
//*************************************************************************************************************
bool NuclearSegmentation::SetInputs(std::string datafile, std::string paramfile)
{
	if( this->dataFilename.compare(datafile) != 0 ||
		this->paramFilename.compare(paramfile) != 0 )	//The names do not match so clear things and start over:
	{
		this->ResetAll();
	}

	this->dataFilename = datafile;
	this->paramFilename = paramfile;
	return true;
}

//***********************************************************************************************************
// Will load the data image into memory (using ftk::Image)
//***********************************************************************************************************
bool NuclearSegmentation::LoadData()
{
	dataImage = ftk::Image::New();
	if(!dataImage->LoadFile(dataFilename))	//Load for display
	{
		errorMessage = "Data Image failed to load";
		dataImage = 0;
		return false;
	}
	return true;
}

//***********************************************************************************************************
// Will load the label image into memory (using ftk::Image)
//***********************************************************************************************************
bool NuclearSegmentation::LoadLabel(bool updateMaps)
{
	labelImage = ftk::Image::New();
	if(!labelImage->LoadFile(labelFilename))	//Load for display
	{
		errorMessage = "Label Image failed to load";
		labelImage = 0;
		return false;
	}
	editsNotSaved = false;

	if(!updateMaps)
		return true;

	//Compute region features:
	//typedef ftk::LabelImageToFeatures< unsigned char, unsigned short, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IPixelT>(0,channelNumber), labelImage->GetItkPtr<LPixelT>(0,0) );
	labFilter->SetLevel(1);
	labFilter->Update();

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
		b.min.z = (int)features->BoundingBox[4];
		b.max.z = (int)features->BoundingBox[5];
		b.min.t = 0;
		b.max.t = 0;

		bBoxMap[(int)id] = b;
		centerMap[(int)id] = c;
		//idToRowMap[(int)id] = r++;
	}
	return true;
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
	int numChannels = info->numChannels;
	if(channelNumber >= numChannels)
		channelNumber = 0;

	//We assume that the image is unsigned char, but just in case it isn't we make it so:
	dataImage->Cast<unsigned char>();
	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,channelNumber,0);	//Expects grayscale image

	if(NucleusSeg) delete NucleusSeg;
	NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(paramFilename.c_str());
	NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
	NucleusSeg->runBinarization();
	lastRunStep = 1;

	//Get the output
	if(getResultImg)
		return GetResultImage();

	return true;
}

bool NuclearSegmentation::DetectSeeds(bool getResultImg)
{
	if(!NucleusSeg)
	{
		errorMessage = "No Binarization";
		return false;
	}
	NucleusSeg->runSeedDetection();
	lastRunStep = 2;
	if(getResultImg)
		return GetResultImage();
	return true;
}

bool NuclearSegmentation::RunClustering()
{
	if(!NucleusSeg || lastRunStep < 2)
	{
		errorMessage = "No Seeds";
		return false;
	}
	NucleusSeg->runClustering();
	lastRunStep = 3;
	this->GetParameters();
	return this->GetResultImage();
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

void NuclearSegmentation::GetParameters()
{
	Parameter p;
	p.name = "sampling_ratio";
	p.value = NucleusSeg->getSamplingRatio();
	this->myParameters.push_back(p);
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
			labelImage->AppendChannelFromData3D(dptr, itk::ImageIOBase::INT, sizeof(unsigned short), size[2], size[1], size[0], "cyan", color, true);		
			Restoredptr(dptr); // Adds the seeds to dptr which were deleted in Cleandptr
		}
		else
		{
			labelImage->AppendChannelFromData3D(dptr, itk::ImageIOBase::INT, sizeof(unsigned short), size[2], size[1], size[0], "cyan", color, true);		
		}
		labelImage->Cast<unsigned short>();
	}
	else
	{
		errorMessage = "Error retrieving data pointer";
		return false;
	}
	return true;
}

bool NuclearSegmentation::SaveResultImage()
{
	if(!labelImage)
	{
		errorMessage = "Nothing To Save";
		return false;
	}

	if(labelFilename.size() == 0)
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

	return true;
}

bool NuclearSegmentation::ComputeFeatures()
{
	if(!dataImage)
	{
		errorMessage = "No Data Image";
		return false;
	}
	if(!labelImage)
	{
		errorMessage = "No Label Image";
		return false;		
	}
	//Compute features:
	//typedef ftk::LabelImageToFeatures< unsigned char, unsigned short, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IPixelT>(0,channelNumber), labelImage->GetItkPtr<LPixelT>(0,0) );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->ComputeTexturesOn();
	labFilter->Update();

	//Init the table (headers):
	featureTable = vtkSmartPointer<vtkTable>::New();		//Start with a new table
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	featureTable->AddColumn(column);
	for (int i=0; i < IntrinsicFeatures::N; ++i)
	{
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (IntrinsicFeatures::Info[i].name).c_str() );
		featureTable->AddColumn(column);
	}

	//Add class column
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "class" );
	featureTable->AddColumn(column); 

	//Add visited column
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "visited?" );
	featureTable->AddColumn(column);

	//Now populate the table:
	//int r = 0;
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant(id) );
		for (int i=0; i< IntrinsicFeatures::N; ++i)
		{
			row->InsertNextValue( vtkVariant(features->ScalarFeatures[i]) );
		}
		int numExtraRows = featureTable->GetNumberOfColumns() - row->GetNumberOfValues();
		for (int i=0; i<numExtraRows; ++i)
		{
			row->InsertNextValue( vtkVariant(-1) );
		}
		featureTable->InsertNextRow(row);
		

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
		//idToRowMap[(int)id] = r++;
	}
	editsNotSaved = true;
	return true;
}

bool NuclearSegmentation::SaveFeaturesTable()
{
	//This function writes the features to a text file that can be read be MetaNeural program
	ofstream outFile; 
	outFile.open(featureFilename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the features
	for(unsigned int row = 0; row < featureTable->GetNumberOfRows(); ++row)
	{
		for(int c=0; c < featureTable->GetNumberOfColumns(); ++c)
		{
			outFile << NumToString( featureTable->GetValue(row,c).ToFloat() ) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();

	//Now print header file:
	outFile.open(headerFilename.c_str(), std::ios::out | std::ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	for(int c=0; c<featureTable->GetNumberOfColumns(); ++c)
	{
		outFile << featureTable->GetColumnName(c) << "\n";
	}
	outFile.close();
	return true;
}

//Save the Edit record
bool NuclearSegmentation::SaveEditRecords()
{
	//This function writes the features to a text file that can be read be MetaNeural program
	ofstream outFile; 
	outFile.open(editFilename.c_str(), ios::out | ios::app );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the edits
	for(unsigned int rec = 0; rec < (int)myEditRecords.size(); ++rec)
	{
		outFile << myEditRecords.at(rec).date << "\t\t" << myEditRecords.at(rec).description << "\n";
	}
	outFile.close();
	return true;
}

//This function will take the data and result files and compute the features
//
bool NuclearSegmentation::LoadFromImages(std::string dfile, std::string rfile)
{
	this->ResetAll();

	//See if data file exists:
	if( !FileExists(dfile) )
	{
		errorMessage = "Could not find data file";
		return 0;
	}
	//See if result file exists
	if( !FileExists(rfile) )
	{
		errorMessage = "Could not find result file";
		return 0;
	}
	//Save the filenames
	dataFilename = dfile;
	labelFilename = rfile;

	if(!LoadData()) return false;
	if(!LoadLabel()) return false;

	return ComputeFeatures();
}

//The .dat file is an old result from the idl-based segmentation
//This function loads from this format for changing to the new format
bool NuclearSegmentation::LoadFromDAT(std::string dfile, std::string rfile)
{
	this->ResetAll();

	if( !FileExists(dfile) )
	{
		errorMessage = "Could not find data file";
		return 0;
	}
	if( !FileExists(rfile) )
	{
		errorMessage = "Could not find result file";
		return 0;
	}
	dataFilename = dfile;
	if(!LoadData()) return false;
	
	ReleaseSegMemory();
	NucleusSeg = new yousef_nucleus_seg();
	const Image::Info *info = dataImage->GetImageInfo();
	int numStacks = info->numZSlices;
	int numRows = info->numRows;				//y-direction
	int numColumns = info->numColumns; 			//x-direction

	//We assume that the image is unsigned char, but just in case it isn't we make it so:
	dataImage->Cast<unsigned char>();
	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,channelNumber,0);		//Expects grayscale image	
	NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
	NucleusSeg->readFromIDLFormat(rfile);
	ReleaseSegMemory();
	this->GetResultImage();

	return ComputeFeatures();
}

//The function will read the given file and load the association measures into the table:
bool NuclearSegmentation::LoadAssociationsFromFile(std::string fName)
{
	if(!FileExists(fName.c_str()))
	{
		errorMessage = "File does not exist";
		return false;
	}

	TiXmlDocument doc;
	if ( !doc.LoadFile( fName.c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return false;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ObjectAssociationRules" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return false;
	}

	std::string source = rootElement->Attribute("SegmentationSource");
	//int numMeasures = atoi( rootElement->Attribute("NumberOfAssociativeMeasures") );
	unsigned int numObjects = atoi( rootElement->Attribute("NumberOfObjects") );

	vtkSmartPointer<vtkTable> assocTable = vtkSmartPointer<vtkTable>::New();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	assocTable->AddColumn(column);

	bool firstRun = true;
	//Parents we know of: Object
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "Object" ) == 0 )
		{
			int id = atoi( parentElement->Attribute("ID") );
			vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
			row->InsertNextValue(vtkVariant(id));

			TiXmlElement *association = parentElement->FirstChildElement();
			while (association)
			{
				std::string name = association->Attribute("Name");
				float value = atof( association->Attribute("Value") );

				if(firstRun)
				{
					column = vtkSmartPointer<vtkDoubleArray>::New();
					column->SetName( name.c_str() );
					assocTable->AddColumn(column);
				}
				
				row->InsertNextValue( vtkVariant(value) );
				association = association->NextSiblingElement();
			}
			assocTable->InsertNextRow(row);
			firstRun = false;
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return false;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)

	//I've created an association table, merge it with featureTable: NEED TO COMPARE IDS and UPDATE IDTOROW MAPPING!!
	if(!featureTable)
	{
		featureTable = assocTable;
	}
	else
	{
		for(int c=1; c<assocTable->GetNumberOfColumns(); ++c)
		{
			featureTable->AddColumn( assocTable->GetColumn(c) );
		}
	}

	return true;
}

bool NuclearSegmentation::SaveLabelByClass()
{
	/*
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
	*/
	return true;
	
}

std::vector<std::string> NuclearSegmentation::GetFeatureNames()
{
	if(!featureTable) 
		return std::vector<std::string>(0);

	std::vector<std::string> names;
	for(int c=1; c<featureTable->GetNumberOfColumns(); ++c)
	{
		names.push_back( featureTable->GetColumnName(c) );
	}
	return names;
}

//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// LOAD UP FORMER RESULTS FROM FILES:
//*******************************************************************************************
bool NuclearSegmentation::LoadAll(std::string filename)
{
	this->ResetAll();

	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return 0;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "NuclearSegmentation" ) != 0 && strcmp( docname, "SegmentationResult" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return 0;
	}
	
	const char* version = rootElement->Attribute("version");
	if(version == 0)
		return RestoreFromXML(filename);

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "datafile" ) == 0 )
		{
			dataFilename = parentElement->GetText();
			channelNumber = atoi(parentElement->Attribute("ch"));
		}
		else if ( strcmp( parent, "resultfile" ) == 0 )
		{
			labelFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "featurefile" ) == 0 )
		{
			featureFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "headerfile" ) == 0 )
		{
			headerFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "paramfile" ) == 0 )
		{
			paramFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "editfile" ) == 0 )
		{
			editFilename = parentElement->GetText();
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return 0;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)

	//doc.close();

	if(!LoadData())
		return false;

	if(!LoadLabel(true))	//load label and update the centroid and bounding boxes
		return false;

	if(!LoadFeatures())
		return false;

	editsNotSaved = false;
	return true;
}

bool NuclearSegmentation::LoadFeatures(void)
{
	if(!FileExists(featureFilename.c_str()))
	{
		errorMessage = "No Features File";
		return false;
	}
	if(!FileExists(headerFilename.c_str()))
	{
		errorMessage = "No Features File";
		return false;
	}

	const int MAXLINESIZE = 1024;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//LOAD THE HEADER INFO:
	ifstream headerFile; 
	headerFile.open( headerFilename.c_str() );
	if ( !headerFile.is_open() )
	{
		errorMessage = "Failed to Load Document";
		return false;
	}

	std::vector< std::string > header; 
	headerFile.getline(line, MAXLINESIZE);
	while ( !headerFile.eof() ) //Get all values
	{
		std::string h;
		char * pch = strtok (line," \t");
		while (pch != NULL)
		{
			h = pch;
			pch = strtok (NULL, " \t");
		}
		header.push_back( h );
		headerFile.getline(line, MAXLINESIZE);
	}
	headerFile.close();

	if ( !header.at(0).compare("ID") && !header.at(0).compare("id") && !header.at(0).compare("Id") )
	{
		errorMessage = "First column must be ID";
		return false;
	}

	//SET UP THE TABLE:
	featureTable = vtkSmartPointer<vtkTable>::New();		//Start with a new table
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	for (int i=0; i < (int)header.size(); ++i)
	{
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( header.at(i).c_str() );
		featureTable->AddColumn(column);
	}

	//LOAD ALL OF THE FEATURES INFO:
	ifstream featureFile; 
	featureFile.open( featureFilename.c_str() );
	if ( !featureFile.is_open() )
	{
		errorMessage = "Failed to Load Document";
		return false;
	}

	featureFile.getline(line, MAXLINESIZE);
	while ( !featureFile.eof() ) //Get all values
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		char * pch = strtok (line," \t");
		while (pch != NULL)
		{
			row->InsertNextValue( vtkVariant( atof(pch) ) );
			pch = strtok (NULL, " \t");
		}
		featureTable->InsertNextRow(row);
		featureFile.getline(line, MAXLINESIZE);
	}
	featureFile.close();
	return true;
}

//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// These functions write the objects to an XML file using particular XML Format!!
//*******************************************************************************************
bool NuclearSegmentation::SaveChanges(std::string filename)
{
	size_t pos = dataFilename.find_last_of(".");
	std::string base = dataFilename.substr(0,pos);

	if(featureFilename.size() == 0)
	{
		featureFilename = base + "_features.txt";
	}
	if(headerFilename.size() == 0)
	{
		headerFilename = base + "_header.txt";
	}
	if(editFilename.size() == 0)
	{
		editFilename = base + "_edits.txt";
	}

	if(this->editsNotSaved)
	{
		SaveResultImage();
		SaveFeaturesTable();
		SaveEditRecords();
		myEditRecords.clear();
	}

	//Write the XML file:
	TiXmlDocument doc;   
 
	TiXmlElement * root = new TiXmlElement( "NuclearSegmentation" );  
	doc.LinkEndChild( root );  
	root->SetAttribute("program", "Yousef_Nucleus_Seg");
	root->SetAttribute("version", "2"); 
 
	if(dataFilename.size() > 0)
	{
		TiXmlElement * dfile = new TiXmlElement("datafile");
		dfile->SetAttribute("ch", NumToString(channelNumber).c_str());
		dfile->LinkEndChild( new TiXmlText( dataFilename.c_str() ) );
		root->LinkEndChild(dfile);
	}
	if(labelFilename.size() > 0)
	{
		TiXmlElement *rfile = new TiXmlElement("resultfile");
		rfile->LinkEndChild( new TiXmlText( labelFilename.c_str() ) );
		root->LinkEndChild(rfile);
	}
	if(headerFilename.size() > 0)
	{
		TiXmlElement *hfile = new TiXmlElement("headerfile");
		hfile->LinkEndChild( new TiXmlText( headerFilename.c_str() ) );
		root->LinkEndChild(hfile);
	}
	if(featureFilename.size() > 0)
	{
		TiXmlElement *ffile = new TiXmlElement("featurefile");
		ffile->LinkEndChild( new TiXmlText( featureFilename.c_str() ) );
		root->LinkEndChild(ffile);
	}
	if(paramFilename.size() > 0)
	{
		TiXmlElement *pfile = new TiXmlElement("featurefile");
		pfile->LinkEndChild( new TiXmlText( paramFilename.c_str() ) );
		root->LinkEndChild(pfile);
	}
	if(paramFilename.size() > 0)
	{
		TiXmlElement *efile = new TiXmlElement("editfile");
		efile->LinkEndChild( new TiXmlText( paramFilename.c_str() ) );
		root->LinkEndChild(efile);
	}

	doc.SaveFile( filename.c_str() );

	this->editsNotSaved = false;
	return true;
}


//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// ADDITIONAL WRITERS:
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// Write out the features to Meta Format!!!
//*******************************************************************************************
bool NuclearSegmentation::WriteToMETA(std::string filename)
{
	//This function writes the features to a text file that can be read be MetaNeural program
	ofstream outFile; 
	outFile.open(filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the features
	for(unsigned int row = 0; row < featureTable->GetNumberOfRows(); ++row)
	{
		for(int c=1; c<featureTable->GetNumberOfColumns(); ++c)
		{
			outFile << NumToString( featureTable->GetValue(row,c).ToFloat() ) << "\t";
		}
		outFile << NumToString( -1 ) << "\t";	//class
		outFile << NumToString( featureTable->GetValue(row,0).ToFloat() ) << endl;    //ID
	}
	outFile.close();

	//Now print header file:
	filename.append(".header");
	outFile.open(filename.c_str(), std::ios::out | std::ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	for(int c=1; c<featureTable->GetNumberOfColumns(); ++c)
	{
		outFile << featureTable->GetColumnName(c) << "\n";
	}
	outFile << "CLASS" << "\n";
	outFile << "ID" << std::endl;

	outFile.close();
	return true;
}

bool NuclearSegmentation::WriteToLibSVM(std::string filename)
{
	if(!featureTable) return false;
	//This function writes the features to a text file that can be read by libsvm program

	ofstream outFile; 
	outFile.open(filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the features
	for(unsigned int row = 0; row < featureTable->GetNumberOfRows(); ++row)
	{
		outFile << NumToString( featureTable->GetValue(row,0).ToFloat() ) << " "; //This should be the class			
		for(int c=1; c<featureTable->GetNumberOfColumns(); ++c)
		{
			outFile << c << ":" << NumToString( featureTable->GetValue(row,c).ToFloat() ) << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
	return true;
}

//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
//**********************************************************************************************************
// A FEW UTILITIES
//**********************************************************************************************************
//Check to see if the file will filename fname exists in 
// the project path.
//**********************************************************************************************************
bool NuclearSegmentation::FileExists(std::string filename)
{
	FILE * pFile = fopen (filename.c_str(),"r");
	if (pFile==NULL)
	{
		return false;
	}
	fclose (pFile);
	return true;
}

std::string NuclearSegmentation::NumToString(double d)
{
	stringstream out;
	out << setprecision(2) << fixed << d;	//Default is to use 2 decimal places
	return out.str();
}

std::string NuclearSegmentation::NumToString(int i)
{
	stringstream out;
	out << i ;	 
	return out.str();
}

std::string NuclearSegmentation::NumToString(double d, int p)
{
	stringstream out;
	out << setprecision(p) << fixed << d;	
	return out.str();
}

std::string NuclearSegmentation::TimeStamp()
{
	time_t rawtime;
	struct tm *timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string dt = asctime(timeinfo);
	size_t end = dt.find('\n');
	dt.erase(end);
	return dt;
}
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
// Find the max ID in the table
//**********************************************************************************************************
long int NuclearSegmentation::maxID(void)
{
	long int max = -1;
	if(featureTable)
	{
		for(int row = 0; row < featureTable->GetNumberOfRows(); ++row)
		{
			long int id = featureTable->GetValue(row,0).ToLong();
			if( id > max ) max = id;
		}
	}
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

	for(int z = region.min.z; z <= region.max.z; ++z)
	{
		for(int r=region.min.y; r <= region.max.y; ++r)
		{
			for(int c=region.min.x; c <= region.max.x; ++c)
			{
				int pix = (int)labelImage->GetPixel(0,0,z,r,c);
				for(int i = 0; i < (int)fromIds.size(); ++i)
				{
					if( pix == fromIds.at(i) )
						labelImage->SetPixel(0,0,z,r,c,toId);
				}
			}
		}
	}
}
//**********************************************************************************************************
//**********************************************************************************************************
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
bool NuclearSegmentation::SetClass(vector<int> ids, int clss)
{
	for(int i=0; i<ids.size(); ++i)
	{
		featureTable->SetValueByName( rowForID(ids.at(i)), "class", vtkVariant(clss) );
	}
	this->editsNotSaved = true;
	return true;
}

bool NuclearSegmentation::MarkAsVisited(vector<int> ids, int val)
{
	for(int i=0; i<ids.size(); ++i)
	{
		featureTable->SetValueByName( rowForID(ids.at(i)), "visited?", vtkVariant(val) );
	}
	this->editsNotSaved = true;
	return true;
}

std::vector< int > NuclearSegmentation::Split(ftk::Object::Point P1, ftk::Object::Point P2)
{
	std::vector <int> ret_ids;
	ret_ids.push_back(0);
	ret_ids.push_back(0);

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
		std::cerr << "Error calculating distance transform: " << err << endl ;
	}
	try
	{
		dt_obj2->Update() ;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error calculating distance transform: " << err << endl ;
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

	std::set<int> add_ids;
	add_ids.insert(newID1);
	add_ids.insert(newID2);
	//Add features for the two new objects
	this->addObjectsToTable(add_ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	this->removeFeatures(objID);				//Remove features of the old object

	//Add an edit record:
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	record.description = std::string("S") + "\t" + NumToString(objID) + "\t" + NumToString(newID1) + "," + NumToString(newID2);
	myEditRecords.push_back(record);

	ret_ids.at(0) = newID1;
	ret_ids.at(1) = newID2;
	return ret_ids;
}


std::vector< int > NuclearSegmentation::SplitAlongZ(int objID, int cutSlice)
{
	std::vector <int> ret_ids;
	ret_ids.push_back(0);
	ret_ids.push_back(0);

	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";	
		return ret_ids;
	}

	std::vector<unsigned short> size = labelImage->Size();
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

	std::set<int> ids;
	ids.insert(newID1);
	ids.insert(newID2);
	//Add features for the two new objects
	this->addObjectsToTable(ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	this->removeFeatures(objID);				//Remove features of the old object

	//ALSO NEED AN EDIT LOG:
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	record.description = std::string("S") + "\t" + NumToString(objID) + "\t" + NumToString(newID1) + "," + NumToString(newID2);
	myEditRecords.push_back(record);

	//return the ids of the two cells resulting from spliting
	ret_ids.at(0) = newID1;
	ret_ids.at(1) = newID2;
	return ret_ids;
}
int NuclearSegmentation::Merge(vector<int> ids)
{
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return 0;
	}

	int newID = maxID() + 1;
	ReassignLabels(ids, newID);					//Assign all old labels to this new label
	ftk::Object::Box region = ExtremaBox(ids);
	this->addObjectToTable(newID, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	for(int i=0; i<ids.size(); ++i)
	{
		removeFeatures(ids.at(i));
	}

	//NEED TO ADD AN EDIT RECORD:
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	record.description = std::string("M") + "\t" + NumToString(ids.at(0));
	for(int i=1; i<ids.size(); ++i)
	{
		record.description.append("," + NumToString(ids.at(i)));
	}
	record.description.append("\t" + NumToString(newID));
	myEditRecords.push_back(record);

	return newID;
}

int NuclearSegmentation::AddObject(int x1, int y1, int z1, int x2, int y2, int z2)
{
	//if no label (segmentation) or no data image is available then return
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";					
		return 0;
	}	
	std::vector<unsigned short> size = labelImage->Size();

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
	this->addObjectToTable(newID, x1, y1, z1, x2, y2, z2);
	
	//NEED TO ADD AN EDIT RECORD:
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	record.description = std::string("A") + "\t" + NumToString(newID);
	myEditRecords.push_back(record);

	return newID;
}

bool NuclearSegmentation::Delete(vector<int> ids)
{
	if(!labelImage) return false;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		ReassignLabel(ids.at(i),0);				//Turn each label in list to zero
		removeFeatures(ids.at(i));				//Remove Features From Table:
		//Create Edit Record:
		ftk::Object::EditRecord record;
		record.date = TimeStamp();
		record.description = std::string("D") + "\t" + NumToString(ids.at(i));
		myEditRecords.push_back(record);
	}

	return true;
}

void NuclearSegmentation::removeFeatures(int ID)
{
	if(featureTable)
	{
		featureTable->RemoveRow( rowForID(ID) );
		//idToRowMap.erase( ID );
		centerMap.erase( ID );
		bBoxMap.erase( ID );
	}
	editsNotSaved = true;
}

int NuclearSegmentation::rowForID(int id)
{
	if(featureTable)
	{
		for(int row = 0; row < featureTable->GetNumberOfRows(); ++row)
		{
			if( id == featureTable->GetValue(row,0).ToInt() )
				return row;
		}
	}
	return -1;
}

//Calculate the features within a specific region of the image for a specific ID, and update the table
bool NuclearSegmentation::addObjectToTable(int ID, int x1, int y1, int z1, int x2, int y2, int z2)
{
	std::set<int> ids;
	ids.insert(ID);
	return addObjectsToTable(ids,x1,y1,z1,x2,y2,z2);
}

bool NuclearSegmentation::addObjectsToTable(std::set<int> IDs, int x1, int y1, int z1, int x2, int y2, int z2)
{
	if(!featureTable)
	{
		errorMessage = "Please Compute All Features First";
		return false;
	}

	FeatureCalcType::Pointer labFilter = computeFeatures(x1,y1,z1,x2,y2,z2);
	if(!labFilter) return false;

	//Add to the table:
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant(id) );
		for (int i=0; i< IntrinsicFeatures::N; ++i)
		{
			row->InsertNextValue( vtkVariant(features->ScalarFeatures[i]) );
		}
		int numExtraRows = featureTable->GetNumberOfColumns() - row->GetNumberOfValues();
		for (int i=0; i<numExtraRows; ++i)
		{
			row->InsertNextValue( vtkVariant(-1) );
		}
		featureTable->InsertNextRow(row);

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
		//idToRowMap[(int)id] = featureTable->GetNumberOfRows()-1;
	}
	editsNotSaved = true;
	return true;

}

//Calculate the features within a specific region of the image and return the filter:
//bool NuclearSegmentation::computeFeatures(int x1, int y1, int z1, int x2, int y2, int z2)
FeatureCalcType::Pointer NuclearSegmentation::computeFeatures(int x1, int y1, int z1, int x2, int y2, int z2)
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

	IImageT::RegionType intRegion;
	IImageT::SizeType intSize;
	IImageT::IndexType intIndex;
	LImageT::RegionType labRegion;
	LImageT::SizeType labSize;
	LImageT::IndexType labIndex;

	intIndex[0] = x1;
	intIndex[1] = y1;
	intIndex[2] = z1;
	intSize[0] = x2 - x1 + 1;
	intSize[1] = y2 - y1 + 1;
	intSize[2] = z2 - z1 + 1;

	labIndex[0] = x1;
	labIndex[1] = y1;
	labIndex[2] = z1;
	labSize[0] = x2 - x1 + 1;
	labSize[1] = y2 - y1 + 1;
	labSize[2] = z2 - z1 + 1;

	intRegion.SetSize(intSize);
    intRegion.SetIndex(intIndex);
    itkIntImg->SetRequestedRegion(intRegion);

    labRegion.SetSize(labSize);
    labRegion.SetIndex(labIndex);
    itkLabImg->SetRequestedRegion(labRegion);

	//Compute features:
	//typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( itkIntImg, itkLabImg );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->ComputeTexturesOn();
	labFilter->Update();

	return labFilter;
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
	  
    //Build the region adjacency graph    
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
// BELOW YOU WILL FIND LEGACY LOADERS TO KEEP SOME MEASURE OF BACKWARDS COMPATABILITY
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
// THESE FUNCTIONS RESTORE FROM AN XML FILE:
//***********************************************************************************************************
bool NuclearSegmentation::RestoreFromXML(std::string filename)
{
	this->ResetAll();

	size_t pos = filename.find_last_of("/\\");
	std::string path = filename.substr(0,pos);

	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return 0;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "NuclearSegmentation" ) != 0 && strcmp( docname, "SegmentationResult" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return 0;
	}

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "datafile" ) == 0 )
		{
			std::string dname = parentElement->GetText();
			dataFilename = path + "/" + dname;
		}
		else if ( strcmp( parent, "resultfile" ) == 0 )
		{
			std::string rname = parentElement->GetText();
			labelFilename = path + "/" + rname;
		}
		else if ( strcmp( parent, "parameter" ) == 0 )
		{
			Parameter p;
			p.name = parentElement->Attribute("name");
			p.value = atoi( parentElement->GetText() );
			myParameters.push_back( p );
		}
		else if ( strcmp( parent, "object" ) == 0 )
		{
			//Object o = parseObject(parentElement);	//no need to store object
			parseObject(parentElement);
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return 0;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)

	//doc.close();

	if(!LoadData())
		return false;

	if(!LoadLabel())
		return false;

	editsNotSaved = true;
	return true;
}

Object NuclearSegmentation::parseObject(TiXmlElement *objectElement)
{
	if ( strcmp( objectElement->Value(), "object" ) != 0 )
	{
		errorMessage = "This is not an object: ";
		errorMessage.append(objectElement->Value());
		return Object("null");
	}

	Object object( objectElement->Attribute("type") );

	int id;
	TiXmlElement *member = objectElement->FirstChildElement();
	while(member)
	{
		const char* memberName = member->Value();
		if ( strcmp( memberName, "id" ) == 0 )
		{
			id = atoi( member->GetText() );
			object.SetId( id );
		}
		else if ( strcmp( memberName, "validity" ) == 0 )
		{
			object.SetValidity( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "duplicated" ) == 0 )
		{
			object.SetDuplicated( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "class" ) == 0 )
		{
			object.SetClass( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "center") == 0 )
		{
			centerMap[id] = parseCenter(member);
			object.SetCentroid( centerMap[id] );
		}
		else if ( strcmp( memberName, "bound") == 0 )
		{
			bBoxMap[id] = parseBound(member);
			object.SetBoundingBox( bBoxMap[id] );
		}
		else if ( strcmp( memberName, "features") == 0 )
		{
		//	object.SetFeatures( parseFeatures(member) );
			parseFeaturesToTable(id, member);
		}
		else if ( strcmp( memberName, "EditHistory") == 0 )
		{
			TiXmlElement *record = member->FirstChildElement();
			while(record)
			{
				if ( strcmp( record->Value(), "record" ) == 0 )
				{
					Object::EditRecord r;
					r.date = record->Attribute("date");
					r.description = record->GetText();
					object.AddEditRecord( r );
					myEditRecords.push_back(r);
				}
				record = record->NextSiblingElement();
			}
		}
		else
		{
			errorMessage = "Unrecognized object member: ";
			errorMessage.append(memberName);
			return Object("null");
		}
		member = member->NextSiblingElement();
	}

	return object;
}

Object::Point NuclearSegmentation::parseCenter(TiXmlElement *centerElement)
{
	Object::Point center;

	TiXmlElement *coord = centerElement->FirstChildElement();
	while(coord)
	{
		const char* coordName = coord->Value();
		if ( strcmp( coordName, "x" ) == 0 )
		{
			center.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "y" ) == 0 )
		{
			center.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "z" ) == 0 )
		{
			center.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "t" ) == 0 )
		{
			center.t = atoi( coord->GetText() );
		}
		coord = coord->NextSiblingElement();
	}
	return center;
}

Object::Box NuclearSegmentation::parseBound(TiXmlElement *boundElement)
{
	Object::Box bound;

	TiXmlElement *coord = boundElement->FirstChildElement();
	while(coord)
	{
		const char* coordName = coord->Value();
		if ( strcmp( coordName, "xmin" ) == 0 )
		{
			bound.min.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "ymin" ) == 0 )
		{
			bound.min.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "zmin" ) == 0 )
		{
			bound.min.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "tmin" ) == 0 )
		{
			bound.min.t = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "xmax" ) == 0 )
		{
			bound.max.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "ymax" ) == 0 )
		{
			bound.max.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "zmax" ) == 0 )
		{
			bound.max.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "tmax" ) == 0 )
		{
			bound.max.t = atoi( coord->GetText() );
		}
		coord = coord->NextSiblingElement();
	}
	return bound;
}

void NuclearSegmentation::parseFeaturesToTable(int id, TiXmlElement *featureElement)
{
	vector< float > tempFeatures(0);
	vector< string > tempNames(0); 

	//Extract all of the features for this object
	TiXmlElement *feature = featureElement->FirstChildElement();
	while(feature)
	{
		tempNames.push_back( feature->Value() );
		tempFeatures.push_back( float( atof( feature->GetText() ) ) );
		
		feature = feature->NextSiblingElement();
	}

	//Initialize the table with headers:
	if( !featureTable )
	{
		featureTable = vtkSmartPointer<vtkTable>::New();
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( "ID" );
		featureTable->AddColumn(column);
		for (int i=0; i < (int)tempNames.size(); ++i)
		{
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( tempNames.at(i).c_str() );
			featureTable->AddColumn(column);
		}
	}
	//Populate Table with features:
	vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
	row->InsertNextValue( vtkVariant(id) );
	for (int i=0; i< IntrinsicFeatures::N; ++i)
	{
		row->InsertNextValue( vtkVariant(tempFeatures.at(i)) );
	}
	featureTable->InsertNextRow(row);
	//idToRowMap[id] = featureTable->GetNumberOfRows()-1;
}

/* OLD METHOD:
vector< float > NuclearSegmentation::parseFeatures(TiXmlElement *featureElement)
{
	vector< float > tempFeatures(0);
	vector< string > tempNames(0); 

	//Extract all of the features for this object
	TiXmlElement *feature = featureElement->FirstChildElement();
	while(feature)
	{
		tempNames.push_back( feature->Value() );
		tempFeatures.push_back( float( atof( feature->GetText() ) ) );
		
		feature = feature->NextSiblingElement();
	}

	//Now test for existing featureNames
	if (featureNames.size() == 0) //Meaning hasn't been set yet
	{
		featureNames = tempNames; //Set my feature names to these names
		return tempFeatures;
	}

	//Have already set the names, so need to make sure we return the features in the correct order
	if( featureNames.size() != tempNames.size() )
	{
		errorMessage = "Features do not match";
		return vector< float >(0);
	}
		
	vector<float> retFeatures( tempFeatures.size() );
	unsigned int numMatches = 0;
	for (unsigned int i=0; i<tempFeatures.size(); ++i)
	{
		string featName = tempNames[i];
		for (unsigned int j=0; j<featureNames.size(); ++j)
		{
			if( featName == featureNames[j] )
			{
				retFeatures[j] = tempFeatures[i];
				++numMatches;
			}
		}
	}

	if ( numMatches != tempFeatures.size() )
	{
		errorMessage = "Features do not match";
		return vector< float >(0);
	}

	return retFeatures;
}
*/
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************


} //END NAMESPACE FTK