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
#include "ftkNuclearSegmentation.h"
#include <itkImageRegionConstIteratorWithIndex.h>
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
	//typedef ftk::LabelImageToFeatures< unsigned char, unsigned short, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( dataImage->GetItkPtr<IPixelT>(0,channelNumber), labelImage->GetItkPtr<LPixelT>(0,0) );
	labFilter->SetLevel(1);
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
		b.min.z = (int)features->BoundingBox[4];
		b.max.z = (int)features->BoundingBox[5];
		b.min.t = 0;
		b.max.t = 0;

		bBoxMap[(int)id] = b;
		centerMap[(int)id] = c;
	}
	return true;
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
	this->addObjectsToMaps(add_ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	this->removeObjectFromMaps(objID);
	EditsNotSaved = true;

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
	this->addObjectsToMaps(ids, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	this->removeObjectFromMaps(objID);
	EditsNotSaved = true;

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
	this->addObjectToMaps(newID, region.min.x, region.min.y, region.min.z, region.max.x, region.max.y, region.max.z);
	for(unsigned int i=0; i<ids.size(); ++i)
	{
		removeObjectFromMaps(ids.at(i));
	}
	EditsNotSaved = true;

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
	this->addObjectToMaps(newID, x1, y1, z1, x2, y2, z2);
	EditsNotSaved = true;

	return newID;
}

bool NuclearSegmentation::Delete(std::vector<int> ids)
{
	if(!labelImage) return false;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		ReassignLabel(ids.at(i),0);				//Turn each label in list to zero
		removeObjectFromMaps(ids.at(i));
	}
	EditsNotSaved = true;

	return true;
}

bool NuclearSegmentation::Exclude(int xy, int z)
{
	if(!labelImage) return false;

	//Find the bounds to use for exclusion margin:
	const ftk::Image::Info *info = labelImage->GetImageInfo();
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;
	int zSlices = (*info).numZSlices;

	int min_x = 0 + xy;
	int min_y = 0 + xy;
	int min_z = 0 + z;
	int max_x = totalWidth - xy - 1;
	int max_y = totalHeight - xy - 1;
	int max_z = zSlices - z - 1;

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
		removeObjectFromMaps(ids.at(i));			
	}
	EditsNotSaved = true;

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
void NuclearSegmentation::removeObjectFromMaps(int ID)
{
	centerMap.erase( ID );
	bBoxMap.erase( ID );
}

//Calculate the features within a specific region of the image for a specific ID, and update the table
bool NuclearSegmentation::addObjectToMaps(int ID, int x1, int y1, int z1, int x2, int y2, int z2)
{
	std::set<int> ids;
	ids.insert(ID);
	return addObjectsToMaps(ids,x1,y1,z1,x2,y2,z2);
}

bool NuclearSegmentation::addObjectsToMaps(std::set<int> IDs, int x1, int y1, int z1, int x2, int y2, int z2)
{
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
	return true;
}

//Calculate the features within a specific region of the image and return the filter:
//bool NuclearSegmentation::computeFeatures(int x1, int y1, int z1, int x2, int y2, int z2)
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
} //END NAMESPACE FTK
