/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkNuclearSegmentation.h"
#include <ctime>

namespace ftk 
{

//Constructor
NuclearSegmentation::NuclearSegmentation(string projpath, string projname):SegmentationResult(projpath, projname)
{
	programName = "Yousef_Nucleus_Seg";

	initConstants();

	dataImage = NULL;
	labelImage = NULL;
	NucleusSeg = NULL;
	editsNotSaved = false;
}

//This function will take the data and result filenames, and create the SegmentationResult structure
// file should be filename only.  Path should come from path set for project
bool NuclearSegmentation::LoadFromResult(const char* dfile, const char* rfile)
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
	dataFilenames.push_back(dfile);
	resultFilenames.push_back(rfile);

	editsNotSaved = false;

	return LabelsToObjects();
}

//Calculate the object information from the data/result images:
bool NuclearSegmentation::LabelsToObjects(void)
{
	if(!dataImage)
	{
		dataImage = new ftk::Image();
		dataImage->LoadFile(PrependProjectPath(dataFilenames[0]));	//Assume there is just one data file and one result file
	}
	if(!labelImage)
	{
		labelImage = new ftk::Image();
		labelImage->LoadFile(PrependProjectPath(resultFilenames[0]));
	}

	//Calculation
	if(!labelImage || !dataImage)
	{
		errorMessage = "label image or data image doesn't exist";
		return false;
	}

	typedef unsigned char IPixelT;
	dataImage->Cast<IPixelT>();
	typedef unsigned short LPixelT;
	labelImage->Cast<LPixelT>();

	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IPixelT>(0,0), labelImage->GetItkPtr<LPixelT>(0,0) );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->ComputeAdvancedOn();
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

	LoadAssociationsFromFile(projectName + "_Associations.txt");
	LoadClassInfoFromFile(projectName + "_class.txt");

	return true;
}

//The function will read the given file and load the association measures into the object info;
void NuclearSegmentation::LoadAssociationsFromFile(std::string fName)
{
	if(myObjects.size() == 0)
		return;

	if(!FileExists(fName.c_str()))
		return;

	ifstream inFile; 
	inFile.open( PrependProjectPath(fName).c_str() );
	if ( !inFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << fName << std::endl;
		return;
	}

	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	std::vector< std::vector< float > > vals(0); 
	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all values
	{
		std::vector< float > lvector(0);
		char * pch = strtok (line," \t");
		while (pch != NULL)
		{
			lvector.push_back( atof(pch) );
			pch = strtok (NULL, " \t");
		}

		vals.push_back( lvector );
		inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();

	std::vector< std::string > possibleNames;
	possibleNames.push_back( "Nissl_sig" );
	possibleNames.push_back( "Iba1_sig" );
	possibleNames.push_back( "GFAP_sig" );
	possibleNames.push_back( "EBA_sig" );

	for (int i=0; i<vals.at(0).size(); ++i)
	{
		featureNames.push_back( possibleNames.at(i) );
	}

	//Now add these features to the objects that we have:
	int n = vals.size() > myObjects.size() ? myObjects.size() : vals.size();

	for (int i=0; i<n; ++i)
	{
		std::vector<float> feats = myObjects.at(i).GetFeatures();
		for (int j=0; j<vals.at(i).size(); ++j)
		{
			feats.push_back(vals.at(i).at(j));
		}
		myObjects.at(i).SetFeatures(feats);
	}
}

void NuclearSegmentation::LoadClassInfoFromFile( std::string fName )
{
	if(myObjects.size() == 0)
		return;

	if(!FileExists(fName.c_str()))
		return;

	ifstream inFile; 
	inFile.open( PrependProjectPath(fName).c_str() );
	if ( !inFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << fName << std::endl;
		return;
	}

	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	std::map< int, int > classNumber; 
	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all rows
	{
		char * pch = strtok (line," \t\n");
		int id = (int)atof(pch);
		pch = strtok (NULL, " \t\n");
		int clss = (int)atof(pch);

		classNumber[id] = clss;
	
		inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();

	std::map< int, int >::iterator it;
	for ( it=classNumber.begin() ; it != classNumber.end(); it++ )

	{
		ftk::Object * obj = GetObjectPtr( (*it).first );
		obj->SetClass( (*it).second );
	}
}

bool NuclearSegmentation::LoadData()
{
	if(dataImage)
	{
		errorMessage = "Data already loaded";
		return 0;
	}
	if( dataFilenames.size() <= 0 )
	{
		errorMessage = "No Data Image specified in XML";
		return 0;
	}

	dataImage = new ftk::Image();
	if(!dataImage->LoadFile(PrependProjectPath(dataFilenames[0]), true))	//Load for display
	{
		delete dataImage;
		errorMessage = "Data Image failed to load";
		return 0;
	}
	return true;
}

bool NuclearSegmentation::LoadLabel()
{
	if(labelImage)
	{
		errorMessage = "Label already loaded";
		return 0;
	}
	if( resultFilenames.size() <= 0 )
	{
		errorMessage = "No Label Image specified in XML";
		return 0;
	}


	labelImage = new ftk::Image();
	if(!labelImage->LoadFile(PrependProjectPath(resultFilenames[0])))
	{
		delete labelImage;
		errorMessage = "Label Image failed to load";
		return 0;
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

	//First check and make sure resultFilename is there
	if(resultFilenames.size() == 0)
	{
		string labelFilename = projectName;
		labelFilename.append("_label.tiff");
		resultFilenames.push_back(labelFilename);
	}

	//Write the labelImage to file
	size_t pos = resultFilenames[0].find(".");
	string base = resultFilenames[0].substr(0,pos);

	labelImage->Cast<unsigned short>();		//Cannot Save as int type to tiff
	if(!labelImage->SaveAs(projectPath, base, "tiff"))
		std::cerr << "FAILED TO SAVE LABEL IMAGE" << std::endl;

	//Added by Yousef on 1/18/2009: save results into a format readable by the IDL farsight
	if(NucleusSeg)
		NucleusSeg->saveIntoIDLFormat(PrependProjectPath(projectName));

	return true;
}

bool NuclearSegmentation::SaveAll()
{
	SaveLabel();
	WriteToXML();
	WriteToMETA();
	editsNotSaved = false;
	return true;
}

ftk::Object::Box NuclearSegmentation::ExtremaBox(vector<int> ids)
{
	ftk::Object::Box extreme = myObjects.at( GetObjectIndex(ids.at(0),"nucleus") ).GetBounds().at(0);
	for(int i=1; i<ids.size(); ++i)
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
	for(int i=0; i<ids.size(); ++i)
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
	labFilter->ComputeAdvancedOn();
	labFilter->Update();

	myObjects.push_back( GetNewObject(newID, labFilter->GetFeatures(newID) ) );
	
	ftk::Object::EditRecord record;
	record.date = TimeStamp();
	std::string msg = "MERGED TO FROM: ";
	msg.append(NumToString(ids.at(0)));
	for(int i=1; i<ids.size(); ++i)
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

	for(int i=0; i<ids.size(); ++i)
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
				int pix = labelImage->GetPixel<int>(0,0,z,r,c);
				if( pix == fromId )
					labelImage->SetPixel<int>(0,0,z,r,c,toId);
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
				int pix = labelImage->GetPixel<int>(0,0,z,r,c);
				for(int i = 0; i < fromIds.size(); ++i)
				{
					if( pix == fromIds.at(i) )
						labelImage->SetPixel<int>(0,0,z,r,c,toId);
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
	for( int i=0; i<myObjects.size(); ++i )
	{
		ftk::Object obj = myObjects[i];
		if (obj.GetId() == objectID && obj.GetType() == type)
		{
			return i;
		}
	}
	return -1;
}



//Check to see if the file will filename fname exists in 
// the project path.
bool NuclearSegmentation::FileExists(const char* fname)
{
	FILE * pFile = fopen (PrependProjectPath(fname).c_str(),"r");
	if (pFile==NULL)
	{
		return false;
	}
	fclose (pFile);
	return true;
}

//********************************************************************************************
//********************************************************************************************
//********************************************************************************************
//LEGACY FUNCTIONS TO INTERACT WITH MODULE WIDGET:
void NuclearSegmentation::initConstants()
{
	//Initialize package and module names and values
	moduleNames.resize(0);
	moduleNames.push_back("Read From .dat file");
	moduleNames.push_back("Binarization");
	moduleNames.push_back("Initial Segmentation");
	moduleNames.push_back("Alpha Expansion");
}
//***************************************************************************
// Call this function first to setup the module with appropiate data
//***************************************************************************
void NuclearSegmentation::setup(string imagefilename, string paramfilename)
{
	if(dataImage)
	{
		delete dataImage;
		dataImage = NULL;
		dataFilenames.clear();
	}
	if(NucleusSeg)
	{
		delete NucleusSeg;
		NucleusSeg = NULL;
	}

	dataImage = new ftk::Image();
	string fname = PrependProjectPath(imagefilename);
	dataImage->LoadFile(fname, true);		//Scale input to 8 bits.
	dataFilenames.push_back(imagefilename);

	Image::Info *info = dataImage->GetImageInfo();
	int numStacks = info->numZSlices;
	int numRows = info->numRows;				//y-direction
	int numColumns = info->numColumns; 			//x-direction

	unsigned char *dataImagePtr = dataImage->GetSlicePtr<unsigned char>(0,0,0);	//Expects grayscale image
	char *f = (char*)imagefilename.c_str();

	NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(PrependProjectPath(paramfilename).c_str());
	NucleusSeg->setDataImage( dataImagePtr, numColumns, numRows, numStacks, f ); 
}

//***************************************************************************
// Execute the specified module
// Assumes all parameters and source data have been supplied
//***************************************************************************
void NuclearSegmentation::executeModule(int moduleNum)
{
	//Check to make sure segmentation has been initialized
	if(!NucleusSeg)
		return;

	vector<int> size = NucleusSeg->getImageSize();

	switch(moduleNum)
	{
	case 0:
		NucleusSeg->readFromIDLFormat(PrependProjectPath(dataFilenames[0]));
		createFTKLabelImg( NucleusSeg->getSegImage(), size[2],size[1],size[0] );
		break;
	case 1:
		NucleusSeg->runBinarization();
		createFTKLabelImg( NucleusSeg->getBinImage(), size[2], size[1], size[0] );
		break;
	case 2:
		NucleusSeg->runSeedDetection();
		NucleusSeg->runClustering();
		createFTKLabelImg( NucleusSeg->getClustImage(), size[2], size[1], size[0] );
		break;
	case 3:
		NucleusSeg->runAlphaExpansion3D();
		createFTKLabelImg( NucleusSeg->getSegImage(), size[2], size[1], size[0] );
		break;
	}
}

//*******************************************************************************
// INTERNAL FUNCTIONS

//creates the ftkImage labelImage that pointer in base class requires
void NuclearSegmentation::createFTKLabelImg(int* data, int numColumns, int numRows, int numStacks)
{
	if(!data)
		return;

	if(!labelImage)
	{
		labelImage = new ftk::Image();
	}
	else
	{
		delete labelImage;
		labelImage = new ftk::Image();
	}
	labelImage->ImageFromData3D((void*)data, INT, sizeof(int), numColumns, numRows, numStacks);
}

} //end namespace ftk