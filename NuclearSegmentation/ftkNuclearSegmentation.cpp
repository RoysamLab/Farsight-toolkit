/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkNuclearSegmentation.h"
#include <itkImageRegionConstIteratorWithIndex.h>
#include <ctime>

#include "yousef_core/graphColLearn_3D/sequential_coloring.cpp"

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

	if(associationFile.size())
		LoadAssociationsFromFile(associationFile);
	if(classFile.size())
		LoadClassInfoFromFile(classFile);

	return true;
}

//The function will read the given file and load the association measures into the object info;
void NuclearSegmentation::LoadAssociationsFromFile(std::string fName)
{
	if(myObjects.size() == 0)
		return;

	if(!FileExists(fName.c_str()))
		return;


	TiXmlDocument doc;
	if ( !doc.LoadFile( fName.c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ObjectAssociationRules" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return;
	}

	std::string source = rootElement->Attribute("SegmentationSource");
	int numMeasures = atoi( rootElement->Attribute("NumberOfAssociativeMeasures") );
	int numObjects = atoi( rootElement->Attribute("NumberOfObjects") );

	if( numObjects != myObjects.size() )
		return;
	
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
			return;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
}

void NuclearSegmentation::LoadClassInfoFromFile( std::string fName )
{
	if(myObjects.size() == 0)
		return;

	if(!FileExists(fName.c_str()))
		return;

	classFile = fName;

	ifstream inFile; 
	inFile.open( PrependProjectPath(fName).c_str() );
	if ( !inFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << fName << std::endl;
		return;
	}

	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//std::map< int, int > classNumber; 
	std::vector<int> classNumber;
	inFile.getline(line, MAXLINESIZE);
	int id;
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
		return;

	std::set<int> classList;
	std::set<int>::iterator it;
	for(int i=0; i<objects->size(); ++i)
	{
		int c = classNumber.at(i);

		it=classList.find(c);
		if(it==classList.end())
			classList.insert(c);

		objects->at(i).SetClass(c);
	}

	classes.clear();
	for(it=classList.begin(); it!=classList.end(); ++it)
	{
		classes.push_back(*it);
	}

}

bool NuclearSegmentation::LoadFromMETA(std::string META_file, std::string header_file, std::string data_file, std::string label_file)
{
	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	//Save the filenames & Load the Images
	dataFilenames.push_back(data_file);
	resultFilenames.push_back(label_file);
	if(dataImage) delete dataImage;
	dataImage = new ftk::Image();
	dataImage->LoadFile(data_file, true);	//Assume there is just one data file and one result file
	if(labelImage) delete labelImage;
	labelImage = new ftk::Image();
	labelImage->LoadFile(label_file, false);

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
			object.SetClass( meta.at(metaRow).at(classColumn) );
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
		std::string fname = projectPath + "/" + projectName;
		fname.append("_class");
		fname.append(NumToString(classes.at(i)));
		fname.append(".tiff");
		writer->SetFileName( fname );
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
	//WriteToMETA();
	WriteToLibSVM();
	editsNotSaved = false;
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
				for(int i = 0; i < (int)fromIds.size(); ++i)
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

//This function runs graph coloring
int NuclearSegmentation::RunGraphColoring(const char* filename)
{
	//get the label image (if not already done)
	if(!labelImage)
	{
		std::cout<<"Loading Label Image ... ";
		labelImage = new ftk::Image();
		labelImage->LoadFile(PrependProjectPath(resultFilenames[0]));
		std::cout<<"done!"<<endl;
	}
    int*** labs_im;
    int max_lab,ncolors;
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
	FILE *fp = fopen(filename,"w");
	
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
		std::string fname = projectPath + "/" + projectName;
		fname.append("_class");
		fname.append(NumToString(classes.at(i)));
		fname.append(".tiff");
		writer->SetFileName( fname );
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
} //end namespace ftk

