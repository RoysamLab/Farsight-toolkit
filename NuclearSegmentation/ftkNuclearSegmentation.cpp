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

#include "yousef_core/graphColLearn_3D/sequential_coloring.cpp"

namespace ftk 
{

//Constructor
NuclearSegmentation::NuclearSegmentation()
{
	dataFilename.clear();		
	labelFilename.clear();
	paramFilename.clear();

	myParameters.clear();
	myObjects.clear();
	featureNames.clear();
	maxID = 0;
	errorMessage.clear();
	classes.clear();

	dataImage = NULL;
	labelImage = NULL;
	NucleusSeg = NULL;
	lastRunStep = 0;
	editsNotSaved = false;
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
		this->dataFilename.clear();
		this->paramFilename.clear();
		dataImage = 0;
		if(NucleusSeg)
		{
			delete NucleusSeg;
			lastRunStep = 0;
			NucleusSeg = 0;
		}
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
	if(dataImage)
	{
		errorMessage = "Data already loaded";
		return false;
	}

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
bool NuclearSegmentation::LoadLabel()
{
	if(labelImage)
	{
		errorMessage = "Label already loaded";
		return false;
	}

	labelImage = ftk::Image::New();
	if(!labelImage->LoadFile(labelFilename))	//Load for display
	{
		errorMessage = "Label Image failed to load";
		labelImage = 0;
		return false;
	}
	return true;
}

bool NuclearSegmentation::Binarize()
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

	//We assume that the image is unsigned char, but just in case it isn't we make it so:
	dataImage->Cast<unsigned char>();
	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(0,0,0);	//Expects grayscale image

	if(NucleusSeg)
		delete NucleusSeg;
	NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(paramFilename.c_str());
	NucleusSeg->setDataImage( dptr, numColumns, numRows, numStacks, dataFilename.c_str() );
	NucleusSeg->runBinarization();
	lastRunStep = 1;

	//Get the output
	GetResultImage();
	return true;
}

bool NuclearSegmentation::DetectSeeds()
{
	if(!NucleusSeg)
	{
		errorMessage = "No Binarization";
		return false;
	}
	NucleusSeg->runSeedDetection();
	lastRunStep = 2;
	GetResultImage();
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
	GetResultImage();
	return true;
}

bool NuclearSegmentation::Finalize()
{
	if(!NucleusSeg || lastRunStep < 3)
	{
		errorMessage = "No Initial Clustering";
		return false;
	}
	NucleusSeg->runAlphaExpansion3D();
	lastRunStep = 4;
	GetResultImage();
	return true;
}

bool NuclearSegmentation::GetResultImage()
{
	if(!NucleusSeg)
	{
		errorMessage = "Nothing To Get";
		return false;
	}

	vector<int> size = NucleusSeg->getImageSize();
	int *dptr = NULL;

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
		if(labelImage)
			labelImage = 0;
		labelImage = ftk::Image::New();
		labelImage->AppendChannelFromData3D(dptr, itk::ImageIOBase::INT, sizeof(int), size[2], size[1], size[0], "gray", color, true);
		labelImage->Cast<unsigned short>();
	}
	return true;
}

bool NuclearSegmentation::SaveOutput()
{
	if(!NucleusSeg || !labelImage)
	{
		errorMessage = "Nothing To Save";
		return false;
	}

	std::string tag;

	switch(lastRunStep)
	{
	case 1:
		tag = "_bin";
		break;
	case 2:	//Seeds:
		tag = "_seed";
		break;
	case 3:
		tag = "_label";
		break;
	case 4:
		tag = "_label";
		break;
	default:
		errorMessage = "Nothing To Save";
		return false;
		break;
	}

	size_t pos = dataFilename.find_last_of(".");
	std::string base = dataFilename.substr(0,pos);
	std::string ext = dataFilename.substr(pos);
	labelImage->SaveChannelAs(0, base + tag, ext );

	return true;
}

//Calculate the object information from the data/result images:
bool NuclearSegmentation::LabelsToObjects(void)
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

	typedef ftk::LabelImageToFeatures< unsigned char, unsigned short, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<unsigned char>(0,0), labelImage->GetItkPtr<unsigned short>(0,0) );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->ComputeTexturesOn();
	labFilter->Update();

	//Set Feature Names
	featureNames.clear();
	for (int i=0; i < IntrinsicFeatures::N; ++i)
	{
		featureNames.push_back( IntrinsicFeatures::Info[i].name );
	}

	//Now populate the objects
	myObjects.clear();
	IdToIndexMap.clear();
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;

		if(id > maxID) maxID = id;

		myObjects.push_back( GetNewObject(id, labFilter->GetFeatures(id) ) );
		IdToIndexMap[id] = (int)myObjects.size() - 1;
	}

	/*
	if(associationFile.size())
		LoadAssociationsFromFile(associationFile);
	if(classFile.size())
		LoadClassInfoFromFile(classFile);
	*/

	return true;
}


//This function will take the data and result filenames, and create the NuclearSegmentation structure
// file should be filename only.  Path should come from path set for project
bool NuclearSegmentation::LoadFromImages(std::string dfile, std::string rfile)
{
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

	LoadData();
	LoadLabel();


	editsNotSaved = false;

	return LabelsToObjects();
}



//The function will read the given file and load the association measures into the object info;
bool NuclearSegmentation::LoadAssociationsFromFile(std::string fName)
{
	if(myObjects.size() == 0)
	{
		errorMessage = "No Objects";
		return false;
	}
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

	if( numObjects != myObjects.size() )
	{
		errorMessage = "Number of Objects does not match";
		return false;
	}

	bool firstRun = true;
	//Parents we know of: Object
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "Object" ) == 0 )
		{
			int id = atoi( parentElement->Attribute("ID") );
			std::vector<float> feats = myObjects.at( IdToIndexMap[id] ).GetFeatures();

			TiXmlElement *association = parentElement->FirstChildElement();
			while (association)
			{
				std::string name = association->Attribute("Name");
				float value = atof( association->Attribute("Value") );

				if(firstRun)
				{
					featureNames.push_back( name );
				}
				feats.push_back(value);

				association = association->NextSiblingElement();
			}
			firstRun = false;
			myObjects.at( IdToIndexMap[id] ).SetFeatures(feats);
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return false;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	return true;
}

bool NuclearSegmentation::LoadClassInfoFromFile( std::string fName )
{
	if(myObjects.size() == 0)
	{
		errorMessage = "No Objects";
		return false;
	}

	if(!FileExists(fName.c_str()))
	{
		errorMessage = "No Class File";
		return false;
	}

	ifstream inFile; 
	inFile.open( fName.c_str() );
	if ( !inFile.is_open() )
	{
		errorMessage = "Failed to Load Document";
		return false;
	}

	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//std::map< int, int > classNumber; 
	std::vector<int> classNumber;
	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all rows
	{
		char * pch = strtok (line," \t\n");
		//int id = (int)atof(pch);
		//pch = strtok (NULL, " \t\n");
		int clss = (int)atof(pch);

		//classNumber[id] = clss;
		classNumber.push_back(clss);
	
		inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();

	/*
	std::map< int, int >::iterator it;
	for ( it=classNumber.begin() ; it != classNumber.end(); it++ )

	{
		ftk::Object * obj = GetObjectPtr( (*it).first );
		obj->SetClass( (*it).second );
	}
	*/

	std::vector<Object> *objects = GetObjectsPtr();
	if(classNumber.size() != objects->size())
	{
		errorMessage = "File sizes do not match";
		return false;
	}

	std::set<int> classList;
	std::set<int>::iterator it;
	for(unsigned int i=0; i<objects->size(); ++i)
	{
		int c = classNumber.at(i);

		it=classList.find(c);
		if(it==classList.end())
			classList.insert(c);

		objects->at(i).SetClass( char(c) );
	}

	classes.clear();
	for(it=classList.begin(); it!=classList.end(); ++it)
	{
		classes.push_back(*it);
	}
	return true;
}

bool NuclearSegmentation::LoadFromMETA(std::string META_file, std::string header_file, std::string data_file, std::string label_file)
{
	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Save the filenames & Load the Images
	dataFilename=data_file;
	labelFilename=label_file;
	dataImage = ftk::Image::New();
	dataImage->LoadFile(data_file);	//Assume there is just one data file and one result file
	labelImage = ftk::Image::New();
	labelImage->LoadFile(label_file);

	//NOW LOAD THE HEADER INFO:
	ifstream headerFile; 
	headerFile.open( header_file.c_str() );
	if ( !headerFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << headerFile << std::endl;
		return 0;
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

	//SEARCH FOR CLASS/RESPONSE AND ID COLUMNS
	//IF I DO NOT FIND ID I ASSUME ORDERED BY MAGNITUDE, AND GET IDS FROM LABEL IMAGE
	int classColumn = -1;
	int idColumn = -1;
	for (int i=0; i<(int)header.size(); ++i)
	{
		if ( !header.at(i).compare("CLASS") || !header.at(i).compare("class")|| !header.at(i).compare("RESPONSE") )
			classColumn = i;
		else if ( !header.at(i).compare("ID") )
			idColumn = i;
	}

	//NOW LOAD ALL OF THE FEATURES INFO:
	ifstream metaFile; 
	metaFile.open( META_file.c_str() );
	if ( !metaFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << metaFile << std::endl;
		return 0;
	}
	std::vector< std::vector<double> > meta; 
	metaFile.getline(line, MAXLINESIZE);
	while ( !metaFile.eof() ) //Get all values
	{
		std::vector<double> row;
		char * pch = strtok (line," \t");
		while (pch != NULL)
		{
			row.push_back( atof(pch) );
			pch = strtok (NULL, " \t");
		}
		meta.push_back( row );
		metaFile.getline(line, MAXLINESIZE);
	}
	metaFile.close();

	//NOW RUN THE FEATURES FILTER TO GET BOUNDING BOXES:
	typedef unsigned char IPixelT;
	dataImage->Cast<IPixelT>();
	typedef unsigned short LPixelT;
	labelImage->Cast<LPixelT>();

	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IPixelT>(0,0), labelImage->GetItkPtr<LPixelT>(0,0) );
	labFilter->SetLevel(1);
	labFilter->ComputeHistogramOff();
	labFilter->ComputeTexturesOff();
	labFilter->Update();

	//Set Feature Names
	featureNames.clear();
	for (int i=0; i<(int)header.size(); ++i)
	{
		if (i == classColumn || i == idColumn) continue;
		featureNames.push_back( header.at(i) );
	}

	//Now populate the objects
	myObjects.clear();
	IdToIndexMap.clear();
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for(int row=0; row < (int)labels.size(); row++)
	{
		FeatureCalcType::LabelPixelType id = labels.at(row);
		if(id == 0) continue;

		if(id > maxID) maxID = id;

		int metaRow = -1;
		if(idColumn != -1)
		{
			//Search meta for the current id
			for(int i=0; i<(int)meta.size(); ++i)
			{
				if(meta.at(i).at(idColumn) == id)
					metaRow = i;
			}
		}
		else
		{
			metaRow = row-1;
		}

		Object object("nucleus");
		object.SetId(id);
		object.SetValidity(1);
		object.SetDuplicated(0);
		if( classColumn != -1 && metaRow != -1 )
			object.SetClass( char(meta.at(metaRow).at(classColumn)) );
		else
			object.SetClass(-1);

		ftk::IntrinsicFeatures *features = labFilter->GetFeatures(id);

		Object::Point c;
		c.x = (int)features->Centroid[0];
		c.y = (int)features->Centroid[1];
		c.z = (int)features->Centroid[2];
		c.t = 0;
		object.AddCenter(c);

		Object::Box b;
		b.min.x = (int)features->BoundingBox[0];
		b.max.x = (int)features->BoundingBox[1];
		b.min.y = (int)features->BoundingBox[2];
		b.max.y = (int)features->BoundingBox[3];
		b.min.z = (int)features->BoundingBox[4];
		b.max.z = (int)features->BoundingBox[5];
		b.min.t = 0;
		b.max.t = 0;
		object.AddBound(b);

		vector< float > f(0);
		for (int i=0; i<(int)header.size(); ++i)
		{
			if (i == classColumn || i == idColumn) continue;
			if(metaRow == -1)
				f.push_back(0.0);
			else
				f.push_back( (float)meta.at(metaRow).at(i) );
		}
		object.SetFeatures( f );
		myObjects.push_back( object );
		IdToIndexMap[id] = (int)myObjects.size() - 1;
	}

	return 1;
}

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

	//Create an image for each class:
	int numClasses = (int)classes.size();
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

	//create lists of object ids in each class:
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

	return true;
	
}

bool NuclearSegmentation::SaveLabel()
{
	if(!labelImage)
	{
		errorMessage = "No Label Image to Save";
		return false;
	}

	labelFilename.clear();
	size_t pos = dataFilename.find_last_of(".");
	std::string base = dataFilename.substr(0,pos);
	std::string ext = dataFilename.substr(pos+1);
	labelFilename = base + "_label" + "." + ext;

	labelImage->Cast<unsigned short>();		//Cannot Save as int type to tiff
	if(!labelImage->SaveChannelAs(0, base + "_label", ext))
		std::cerr << "FAILED TO SAVE LABEL IMAGE" << std::endl;

	//Added by Yousef on 1/18/2009: save results into a format readable by the IDL farsight
	if(NucleusSeg)
		NucleusSeg->saveIntoIDLFormat(base);

	return true;
}


ftk::Object::Box NuclearSegmentation::ExtremaBox(vector<int> ids)
{
	ftk::Object::Box extreme = myObjects.at( GetObjectIndex(ids.at(0),"nucleus") ).GetBounds().at(0);
	for(int i=1; i<(int)ids.size(); ++i)
	{
		ftk::Object::Box test = myObjects.at( GetObjectIndex(ids.at(i),"nucleus") ).GetBounds().at(0);;

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

int NuclearSegmentation::Merge(vector<int> ids)
{
	if(!labelImage)
		return 0;

	//Merge is difficult because we are going to invalidate all merged objects,
	// create a new label and assign the old object labels to it,
	// calculate the features for this new label
	// then the model should also be updated to display properly
	// Return the new ID

	int newID = ++maxID;
	for(int i=0; i<(int)ids.size(); ++i)
	{
		int index = GetObjectIndex(ids.at(i),"nucleus");
		if(index < 0 ) return 0;

		myObjects.at(index).SetValidity(false);
		ftk::Object::EditRecord record;
		record.date = TimeStamp();
		std::string msg = "MERGED TO BECOME ";
		msg.append(NumToString(newID));
		record.description = msg;
		myObjects.at(index).AddEditRecord(record);
	}
	ftk::Object::Box region = ExtremaBox(ids);
	ReassignLabels(ids, newID, region);		//Assign all old labels to this new label

	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";
		return 0;
	}

	//Calculate features using feature filter
	typedef unsigned char IPixelT;
	typedef unsigned short LPixelT;
	typedef itk::Image< IPixelT, 3 > IImageT;
	typedef itk::Image< LPixelT, 3 > LImageT;

	dataImage->Cast<IPixelT>();
	labelImage->Cast<LPixelT>();

	IImageT::Pointer itkIntImg = dataImage->GetItkPtr<IPixelT>(0,0);
	LImageT::Pointer itkLabImg = labelImage->GetItkPtr<LPixelT>(0,0);

	IImageT::RegionType intRegion;
	IImageT::SizeType intSize;
	IImageT::IndexType intIndex;
	LImageT::RegionType labRegion;
	LImageT::SizeType labSize;
	LImageT::IndexType labIndex;

	intIndex[0] = region.min.x;
	intIndex[1] = region.min.y;
	intIndex[2] = region.min.z;
	intSize[0] = region.max.x - region.min.x + 1;
	intSize[1] = region.max.y - region.min.y + 1;
	intSize[2] = region.max.z - region.min.z + 1;

	labIndex[0] = region.min.x;
	labIndex[1] = region.min.y;
	labIndex[2] = region.min.z;
	labSize[0] = region.max.x - region.min.x + 1;
	labSize[1] = region.max.y - region.min.y + 1;
	labSize[2] = region.max.z - region.min.z + 1;

	intRegion.SetSize(intSize);
    intRegion.SetIndex(intIndex);
    itkIntImg->SetRequestedRegion(intRegion);

    labRegion.SetSize(labSize);
    labRegion.SetIndex(labIndex);
    itkLabImg->SetRequestedRegion(labRegion);

	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( itkIntImg, itkLabImg );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->ComputeTexturesOn();
	labFilter->Update();

	myObjects.push_back( GetNewObject(newID, labFilter->GetFeatures(newID) ) );
	
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	std::string msg = "MERGED TO FROM: ";
	msg.append(NumToString(ids.at(0)));
	for(int i=1; i<(int)ids.size(); ++i)
	{
		msg.append(", ");
		msg.append(NumToString(ids.at(i)));
	}
	record.description = msg;
	myObjects.back().AddEditRecord(record);
	IdToIndexMap[newID] = (int)myObjects.size() - 1;

	editsNotSaved = true;
	return newID;
}

bool NuclearSegmentation::Delete(vector<int> ids)
{
	if(!labelImage)
		return false;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		//Find the object with the ID
		int index = GetObjectIndex(ids.at(i),"nucleus");
		if (index < 0) return false;

		// 1. Invalidate
		myObjects.at( index ).SetValidity(false);
		// 2. Add to Edit Record
		ftk::Object::EditRecord record;
		record.date = TimeStamp();
		record.description = "DELETED";
		myObjects.at( index ).AddEditRecord(record);
		ReassignLabel(ids.at(i),0);
	}

	editsNotSaved = true;
	return true;
}

ftk::Object NuclearSegmentation::GetNewObject(int id, IntrinsicFeatures *features )
{
	Object object("nucleus");
	object.SetId(id);
	object.SetValidity(1);
	object.SetDuplicated(0);
	object.SetClass(-1);

	if(features == NULL)
		return object;

	Object::Point c;
	c.x = (int)features->Centroid[0];
	c.y = (int)features->Centroid[1];
	c.z = (int)features->Centroid[2];
	c.t = 0;
	object.AddCenter(c);

	Object::Box b;
	b.min.x = (int)features->BoundingBox[0];
	b.max.x = (int)features->BoundingBox[1];
	b.min.y = (int)features->BoundingBox[2];
	b.max.y = (int)features->BoundingBox[3];
	b.min.z = (int)features->BoundingBox[4];
	b.max.z = (int)features->BoundingBox[5];
	b.min.t = 0;
	b.max.t = 0;
	object.AddBound(b);

	vector< float > f(0);
	for (int i=0; i< IntrinsicFeatures::N; ++i)
	{
		f.push_back( features->ScalarFeatures[i] );
	}

	object.SetFeatures( f );

	return object;
}


//This function finds the bounding box of the object with id "fromId",
// and uses that area to change the pixels to "toId"
void NuclearSegmentation::ReassignLabel(int fromId, int toId)
{
	int C = labelImage->Size()[3];
	int R = labelImage->Size()[2];
	int Z = labelImage->Size()[1];

	ftk::Object::Box region = myObjects.at( GetObjectIndex(fromId,"nucleus") ).GetBounds().at(0);

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
				if( pix == fromId )
					labelImage->SetPixel(0,0,z,r,c,toId);
			}
		}
	}
}

//This function finds the bounding box of the objects with ids "fromIds",
// and uses that area to change the pixels to "toId" in 1 pass through the whole region
void NuclearSegmentation::ReassignLabels(vector<int> fromIds, int toId, ftk::Object::Box region)
{
	int C = labelImage->Size()[3];
	int R = labelImage->Size()[2];
	int Z = labelImage->Size()[1];

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

string NuclearSegmentation::TimeStamp()
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

int NuclearSegmentation::GetObjectIndex(int objectID, string type)
{
	//Find the object with the ID
	for( int i=0; i<(int)myObjects.size(); ++i )
	{
		ftk::Object obj = myObjects[i];
		if (obj.GetId() == objectID && obj.GetType() == type)
		{
			return i;
		}
	}
	return -1;
}



//This function runs graph coloring
bool NuclearSegmentation::RunGraphColoring(std::string labelname, std::string filename)
{
	//get the label image (if not already done)
	std::cout<<"Loading Label Image ... ";
	labelImage = ftk::Image::New();
	labelImage->LoadFile(labelname);
	std::cout<<"done!"<<endl;

    //int*** labs_im;
    int max_lab;
    //int ncolors;
    int** RAG;    
    int* ColorOut;        
	int L, L1, L2, L3, L4, L5, L6, L7;
	
	int c = labelImage->Size()[3];
	int r = labelImage->Size()[2];
	int z = labelImage->Size()[1];	
	unsigned short* labs_vals = static_cast<unsigned short*> (labelImage->GetDataPtr(0,0));

	//get the maximum label
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
		fprintf(stderr,"can't open %s for writing\n",filename);
		exit(1);
	}
	for(int i=0; i<max_lab; i++)
	{
		fprintf(fp,"%d\n",ColorOut[i]+1);
		
	}
	fclose(fp);
	//Try this: save the colors into the classes list
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
	std::cout<<"done!"<<std::endl;
	return 1;
}

bool NuclearSegmentation::RestoreFromXML(std::string filename)
{
	dataFilename.clear();		
	labelFilename.clear(); 
	myParameters.clear();
	myObjects.clear();
	featureNames.clear();

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
			Object o = parseObject(parentElement);
			if ( o.GetType() != "null" )
			{
				myObjects.push_back( o );
				int id = o.GetId();
				IdToIndexMap[id] = (int)myObjects.size() - 1;
				if(id > maxID) maxID = id;
			}
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

	TiXmlElement *member = objectElement->FirstChildElement();
	while(member)
	{
		const char* memberName = member->Value();
		if ( strcmp( memberName, "id" ) == 0 )
		{
			object.SetId( atoi( member->GetText() ) );
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
			object.AddCenter( parseCenter(member) );
		}
		else if ( strcmp( memberName, "bound") == 0 )
		{
			object.AddBound( parseBound(member) );
		}
		else if ( strcmp( memberName, "features") == 0 )
		{
			object.SetFeatures( parseFeatures(member) );
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

//*******************************************************************************************
// This functions writes the objects to an XML file using a particular XML Format!!
//*******************************************************************************************
bool NuclearSegmentation::WriteToXML(std::string filename)
{
	TiXmlDocument doc;   
 
	TiXmlElement * root = new TiXmlElement( "NuclearSegmentation" );  
	doc.LinkEndChild( root );  
	root->SetAttribute("program", "Yousef_Nucleus_Seg");

	TiXmlComment * comment = new TiXmlComment();
	comment->SetValue(" Segmentation Results/Parameters/Features/Edits " );  
	root->LinkEndChild( comment );  
 
	TiXmlElement * dfile = new TiXmlElement("datafile");
	dfile->LinkEndChild( new TiXmlText( dataFilename.c_str() ) );
	root->LinkEndChild(dfile);

	TiXmlElement *rfile = new TiXmlElement("resultfile");
	rfile->LinkEndChild( new TiXmlText( labelFilename.c_str() ) );
	root->LinkEndChild(rfile);

	//Attach parameters
	for (unsigned int pnum = 0; pnum < myParameters.size(); ++pnum)
	{
		TiXmlElement *element = new TiXmlElement("parameter");
		element->SetAttribute( "name", myParameters[pnum].name );
		element->LinkEndChild( new TiXmlText( NumToString(myParameters[pnum].value) ) );
		root->LinkEndChild(element);

	}

	//Attach Objects
	for (unsigned int onum=0; onum<myObjects.size(); ++onum)
	{
		root->LinkEndChild( GetObjectElement(myObjects[onum]) );
	}

	doc.SaveFile( filename.c_str() );

	return true;
}

TiXmlElement* NuclearSegmentation::GetObjectElement(Object object)
{
	TiXmlElement *objectElement = new TiXmlElement("object");
	objectElement->SetAttribute("type", object.GetType());

	if ( object.GetId() != -1 )
	{
		TiXmlElement *idElement = new TiXmlElement("id");
		idElement->LinkEndChild( new TiXmlText( NumToString(object.GetId()) ) );
		objectElement->LinkEndChild(idElement);
	}
	if ( object.GetValidity() != -1 )
	{
		TiXmlElement *vElement = new TiXmlElement("validity");
		vElement->LinkEndChild( new TiXmlText( NumToString(object.GetValidity()) ) );
		objectElement->LinkEndChild(vElement);
	}
	if ( object.GetDuplicated() != -1 )
	{
		TiXmlElement *dElement = new TiXmlElement("duplicated");
		dElement->LinkEndChild( new TiXmlText( NumToString(object.GetDuplicated()) ) );
		objectElement->LinkEndChild(dElement);
	}
	if ( object.GetClass() != -1 )
	{
		TiXmlElement *cElement = new TiXmlElement("class");
		cElement->LinkEndChild( new TiXmlText( NumToString(object.GetClass()) ) );
		objectElement->LinkEndChild(cElement);
	}
	vector<Object::Point> centers = object.GetCenters();
	for (unsigned int c=0; c<centers.size(); ++c)
	{
		TiXmlElement *oElement = new TiXmlElement("center");
		TiXmlElement *xElement = new TiXmlElement("x");
		xElement->LinkEndChild( new TiXmlText( NumToString(centers[c].x) ) );
		oElement->LinkEndChild(xElement);
		TiXmlElement *yElement = new TiXmlElement("y");
		yElement->LinkEndChild( new TiXmlText( NumToString(centers[c].y) ) );
		oElement->LinkEndChild(yElement);
		TiXmlElement *zElement = new TiXmlElement("z");
		zElement->LinkEndChild( new TiXmlText( NumToString(centers[c].z) ) );
		oElement->LinkEndChild(zElement);
		TiXmlElement *tElement = new TiXmlElement("t");
		tElement->LinkEndChild( new TiXmlText( NumToString(centers[c].t) ) );
		oElement->LinkEndChild(tElement);
		objectElement->LinkEndChild(oElement);
	}
	vector<Object::Box> bounds = object.GetBounds();
	for (unsigned int b=0; b<bounds.size(); ++b)
	{
		TiXmlElement *bElement = new TiXmlElement("bound");
		TiXmlElement *xminElement = new TiXmlElement("xmin");
		xminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.x) ) );
		bElement->LinkEndChild(xminElement);
		TiXmlElement *yminElement = new TiXmlElement("ymin");
		yminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.y) ) );
		bElement->LinkEndChild(yminElement);
		TiXmlElement *zminElement = new TiXmlElement("zmin");
		zminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.z) ) );
		bElement->LinkEndChild(zminElement);
		TiXmlElement *tminElement = new TiXmlElement("tmin");
		tminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.t) ) );
		bElement->LinkEndChild(tminElement);
		TiXmlElement *xmaxElement = new TiXmlElement("xmax");
		xmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.x) ) );
		bElement->LinkEndChild(xmaxElement);
		TiXmlElement *ymaxElement = new TiXmlElement("ymax");
		ymaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.y) ) );
		bElement->LinkEndChild(ymaxElement);
		TiXmlElement *zmaxElement = new TiXmlElement("zmax");
		zmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.z) ) );
		bElement->LinkEndChild(zmaxElement);
		TiXmlElement *tmaxElement = new TiXmlElement("tmax");
		tmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.t) ) );
		bElement->LinkEndChild(tmaxElement);
		objectElement->LinkEndChild(bElement);
	}
	if(featureNames.size() > 0)
	{
		TiXmlElement *fsElement = new TiXmlElement("features");
		vector<float> features = object.GetFeatures();
		for(unsigned int f=0; f<features.size(); ++f)
		{
			TiXmlElement *fElement = new TiXmlElement( featureNames[f].c_str() );
			fElement->LinkEndChild( new TiXmlText( NumToString(features[f]) ) );
			fsElement->LinkEndChild(fElement);
		}
		objectElement->LinkEndChild(fsElement);
	}
	vector<Object::EditRecord> records = object.getHistory();
	if(records.size() > 0)
	{
		TiXmlElement *histElement = new TiXmlElement("EditHistory");
		for (unsigned int r=0; r<records.size(); ++r)
		{
			TiXmlElement *recElement = new TiXmlElement("record");
			recElement->SetAttribute("date",records[r].date);
			recElement->LinkEndChild( new TiXmlText( records[r].description ) );
			histElement->LinkEndChild(recElement);
		}
		objectElement->LinkEndChild(histElement);
	}
	return objectElement;
}

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
	for(unsigned int obj = 0; obj < myObjects.size(); ++obj)
	{
		if( myObjects.at(obj).GetValidity() == false )
			continue;

		vector<float> feats = myObjects.at(obj).GetFeatures();
		for(unsigned int f = 0; f < feats.size(); ++f)
		{
			outFile << NumToString(feats.at(f)) << "\t";
		}
		outFile << (int)myObjects.at(obj).GetClass() << "\t";
		outFile << myObjects.at(obj).GetId() << endl;
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

	for(int i=0; i<(int)featureNames.size(); ++i)
	{
		outFile << featureNames.at(i) << "\n";
	}

	outFile << "CLASS" << "\n";
	outFile << "ID" << std::endl;

	outFile.close();

	return true;
}

bool NuclearSegmentation::WriteToLibSVM(std::string filename)
{
	//This function writes the features to a text file that can be read by libsvm program

	ofstream outFile; 
	outFile.open(filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the features
	for(unsigned int obj = 0; obj < myObjects.size(); ++obj)
	{
		if( myObjects.at(obj).GetValidity() == false )
			continue;

		outFile << myObjects.at(obj).GetId() << " ";			//This should be the class

		vector<float> feats = myObjects.at(obj).GetFeatures();
		for(unsigned int f = 0; f < feats.size(); ++f)
		{
			outFile << f+1 << ":" << NumToString(feats.at(f)) << " ";	//FeatureNumber:FeatureValue
		}
		outFile << std::endl;
	}
	outFile.close();
	return true;
}


Object* NuclearSegmentation::GetObjectPtr(int id)
{
	int index = IdToIndexMap[id];
	return &myObjects.at(index);
}

//Check to see if the file will filename fname exists in 
// the project path.
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

string NuclearSegmentation::NumToString(double d)
{
	stringstream out;
	out << setprecision(2) << fixed << d;	//Default is to use 2 decimal places
	return out.str();
}

string NuclearSegmentation::NumToString(int i)
{
	stringstream out;
	out << i ;	 
	return out.str();
}

string NuclearSegmentation::NumToString(double d, int p)
{
	stringstream out;
	out << setprecision(p) << fixed << d;	
	return out.str();
}

} //end namespace ftk

