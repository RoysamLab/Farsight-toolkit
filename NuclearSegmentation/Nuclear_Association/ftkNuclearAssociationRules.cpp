/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkNuclearAssociationRules.h"

//#include "itkImageRegionConstIterator.h"
//#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace ftk 
{	
NuclearAssociationRules::NuclearAssociationRules(string AssocFName, int numOfRules):ObjectAssociation(AssocFName, numOfRules)
{
	labImage= NULL;
	labGeometryFilter= NULL;	
	x_Size = y_Size = z_Size = 0;
	imDim=3;	
	objectType = "Nucleus";
}

/* This is the main function for computing associative features */
void NuclearAssociationRules::Compute()
{
	cout<<"Starting Associative Features Computation\n";

	//1. read the label image	
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (GetSegImgName());
	reader->Update();	
	labImage = LabImageType::New();
	labImage = reader->GetOutput();
	x_Size=labImage->GetLargestPossibleRegion().GetSize()[0];
	y_Size=labImage->GetLargestPossibleRegion().GetSize()[1];
	z_Size=labImage->GetLargestPossibleRegion().GetSize()[2];
	if(z_Size <= 1)
		imDim = 2;
	
	//2. Use the labelGeometryImageFilter to get the bounding box for each cell	
	labGeometryFilter = LabelGeometryType::New();
	labGeometryFilter->SetInput( labImage );
	//labGeometryFilter->SetIntensityInput( labImage );
	labGeometryFilter->Update();
	//get the number of objects	
	numOfLabels  = labGeometryFilter->GetNumberOfObjects();
	numOfLabels--; //because the first one is the background

	//allocate memory for the features list
	assocMeasurementsList = new float*[GetNumofAssocRules()];
	//3. then, for each type of the associations get the requested reigion based on the type and the value of the inside and outside distances.
	for(unsigned int i=0; i<GetNumofAssocRules(); i++)
	{
		assocMeasurementsList[i] = new float[numOfLabels];
		//read the ith target image (the image from which we need to compute the ith assoc. rule
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName (assocRulesList[i].GetTargetFileNmae());
		reader2->Update();		
		//cout<<"Computing Features For Association Rule "<<i+1<<": ";
		for(unsigned int j=0; j<numOfLabels; j++)
		{
			//cout<<j+1;
			cout<<"\rComputing Features For Association Rule "<<i+1<<": "<<j+1<<"/"<<numOfLabels;
			assocMeasurementsList[i][j] = ComputeOneAssocMeasurement(reader2->GetOutput(), i, j+1);
			/*std::ostringstream o;
			o<<j+1;
			for(int k=0; k<o.str().size(); k++)
				cout<<"\b";*/
		}		
		cout<<"\tdone"<<endl;
	}						
}

/* use this function to compute one associative measurement */
float NuclearAssociationRules::ComputeOneAssocMeasurement(itk::SmartPointer<TargImageType> trgIm, int ruleID, int objID)
{	
	//fisrt, get the bounding box around the object
	//The bounding box is defined by the area around the object such that the object+outerdistance are included	
	int imBounds[6] = {0,x_Size,0,y_Size,0,z_Size};
	LabelGeometryType::BoundingBoxType bbox = labGeometryFilter->GetBoundingBox(objID);
	std::vector< int > retBbox(0);
	int dist, val;
	dist = assocRulesList[ruleID].GetOutDistance();
	for(unsigned int dim=0; dim < imDim*2; ++dim)
	{  
		if(dim%2 == 0) //even
		{			
			val = int(bbox[dim])-dist-1;
			if(val<imBounds[dim])
				val=imBounds[dim];
		}
		else //odd
		{
			val = int(bbox[dim])+dist+1;
			if(val>=imBounds[dim])
				val=imBounds[dim]-1;
		}
		retBbox.push_back( val );
	}

	//the bounding box defines the region of interest we need (from both segmentation and target images)
	//so, use it to get sub images
	DistImageType::Pointer subSegImg = DistImageType::New();
	
	LabImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
	
    LabImageType::SizeType  size;
    size[0]  = retBbox[1]-retBbox[0]+1;  // size along X
    size[1]  = retBbox[3]-retBbox[2]+1;  // size along Y
	size[2]  = retBbox[5]-retBbox[4]+1;  // size along Z
  
    LabImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
	subSegImg->SetRegions( region );	
    subSegImg->Allocate();
    subSegImg->FillBuffer(0);
	subSegImg->Update();	
	

	LabImageType::IndexType start2;
    start2[0] =   retBbox[0];  // first index on X
    start2[1] =   retBbox[2];  // first index on Y    
	start2[2] =   retBbox[4];  // first index on Z    
	LabImageType::RegionType region2;
    region2.SetSize( size );
    region2.SetIndex( start2 );
	labImage->SetRequestedRegion(region2);
	trgIm->SetRequestedRegion(region2);

	typedef itk::ImageRegionIteratorWithIndex< LabImageType > IteratorType;	
	IteratorType iterator1(labImage, labImage->GetRequestedRegion());
	typedef itk::ImageRegionIteratorWithIndex< DistImageType > IteratorType2;	
	IteratorType2 iterator2(subSegImg, subSegImg->GetRequestedRegion());

	//in the sub-segmentation image, we need to mask out any pixel from another object	
	while ( ! iterator1.IsAtEnd())
	{
		int V = iterator1.Get();
		if(V == objID)
			iterator2.Set(255.0);
		else
			iterator2.Set(0.0);
		++iterator1;
		++iterator2;		
	}	
	//Compute the distance transform in the sub-segmentation image region	
	DTFilter::Pointer dt_obj= DTFilter::New() ;
	dt_obj->SetInput(subSegImg) ;
	dt_obj->SetInsideValue(255.0);
	dt_obj->SetOutsideValue(0.0);	
	try{
		dt_obj->Update() ;
	}
	catch( itk::ExceptionObject & err ){		
		return 0;
	}
	
	//now, mask out all the pixels (in the sub-seg image) that are not in the region of interest as defined by the association rule and get the intensities of the needed pixels from the target image. The intensities are saved into an std vector
	IteratorType2 iterator3(dt_obj->GetOutput(), dt_obj->GetOutput()->GetRequestedRegion());
	IteratorType iterator4(trgIm, trgIm->GetRequestedRegion());
	vector<int> trgInt;
	while ( ! iterator3.IsAtEnd())
	{
		int V = iterator3.Get();
		//if it is outside with distance less than outDistance away
		if(V>0 && V<=assocRulesList[ruleID].GetOutDistance())
			trgInt.push_back(iterator4.Get());
		//if it is inside and the whole cell is used
		else if(V<=0 && assocRulesList[ruleID].IsUseWholeObject())
			trgInt.push_back(iterator4.Get());
		//if it is inside with distance less than in Distance
		else if(V<=0 && abs(V)<=assocRulesList[ruleID].GetInDistance())
			trgInt.push_back(iterator4.Get());

		++iterator3;
		++iterator4;
	}
	if(!trgInt.size())
		return 0;
	
	//Finally, given the list of intensities compute the associative measrement based on the type defined by the association rule
	switch(assocRulesList[ruleID].GetAssocType())
	{
		case ASSOC_MIN:
			return FindMin(trgInt);
			break;
		case ASSOC_MAX:
			return FindMax(trgInt);
			break;
		case ASSOC_TOTAL:
			return ComputeTotal(trgInt);
			break;
		case ASSOC_AVERAGE:
			return ComputeAverage(trgInt);
			break;
		default:
			return ComputeAverage(trgInt);
	}	
}

/* the next functions will compute the measurements for the different cases */
//find maximum
float NuclearAssociationRules::FindMax(vector<int> LST)
{
	float mx = LST[0];
	for(int i=1; i<LST.size(); i++)
	{
		if(mx<LST[i])
			mx = LST[i];
	}
	return mx;
}

//find minimum
float NuclearAssociationRules::FindMin(vector<int> LST)
{
	float mn = LST[0];
	for(int i=1; i<LST.size(); i++)
	{
		if(mn>LST[i])
			mn = LST[i];
	}
	return mn;
}

//compute total
float NuclearAssociationRules::ComputeTotal(vector<int> LST)
{
	float tl = LST[0];
	for(int i=1; i<LST.size(); i++)
	{		
		tl += LST[i];
	}
	return tl;
}

//compute average
float NuclearAssociationRules::ComputeAverage(vector<int> LST)
{
	float av = LST[0];
	for(int i=1; i<LST.size(); i++)
	{		
		av += LST[i];
	}
	av/=LST.size();	
	return av;
}

} //end namespace ftk

