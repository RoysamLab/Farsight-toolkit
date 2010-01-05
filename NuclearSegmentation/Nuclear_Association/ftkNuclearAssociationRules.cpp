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
#include "ftkNuclearAssociationRules.h"

//#include "itkImageRegionConstIterator.h"
//#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>
#include <fstream>
#include <iomanip>
//just for testing
#include "itkImageFileWriter.h"

namespace ftk 
{	

AssociativeFeatureCalculator::AssociativeFeatureCalculator()
{
	inFilename = "";
	fPrefix = "";
	lab_im_set = false;
}

void AssociativeFeatureCalculator::SetInputFile(std::string filename)
{
	inFilename = filename;
}
	
void AssociativeFeatureCalculator::SetFeaturePrefix(std::string prefix)
{
	fPrefix = prefix;
}

//Compute features turned ON in doFeat and puts them in a new table
vtkSmartPointer<vtkTable> AssociativeFeatureCalculator::Compute(void)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc = new ftk::NuclearAssociationRules("",0);
	assoc->ReadRulesFromXML(inFilename);
	assoc->PrintSelf();
	assoc->Compute();

	//Init the table (headers):
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	table->AddColumn(column);
	for (int i=0; i < assoc->GetNumofAssocRules(); ++i)
	{
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (fPrefix+assoc->GetAssociationRules().at(i).GetRuleName()).c_str() );
		table->AddColumn(column);
	}

	//Now populate the table:
	std::vector<unsigned short> labels = assoc->GetLabels();
	float** vals = assoc->GetAssocFeaturesList();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		unsigned short id = labels.at(i);
		if(id == 0) continue;

		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant(id) );
		for (int r=0; r<assoc->GetNumofAssocRules(); ++r)
		{
			row->InsertNextValue( vtkVariant(vals[r][i]) );
		}
		table->InsertNextRow(row);
	}
	return table;
}

//Update the features in this table whose names match (sets doFeat)
void AssociativeFeatureCalculator::Update(vtkSmartPointer<vtkTable> table)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc = new ftk::NuclearAssociationRules("",0);
	assoc->ReadRulesFromXML(inFilename);
	assoc->PrintSelf();
	assoc->Compute();

	//Now update the table:
	std::vector<unsigned short> labels = assoc->GetLabels();
	float** vals = assoc->GetAssocFeaturesList();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		unsigned short id = labels.at(i);
		if(id == 0) continue;

		int row = -1;
		for(int r=0; r<table->GetNumberOfRows(); ++r)
		{
			if( table->GetValue(r,0) == id )
			{
				row = r;
				break;
			}
		}

		if(row == -1) continue;

		for (int f=0; f<assoc->GetNumofAssocRules(); ++f)
		{
			table->SetValueByName(row,(fPrefix+assoc->GetAssociationRules().at(f).GetRuleName()).c_str(), vtkVariant(vals[f][i]));
		}
	}
}

void AssociativeFeatureCalculator::SetInputImage(ftk::Image::Pointer input_labeled_image, int channel_number){

	typedef itk::Image< unsigned short, 3 > UShortImageType3D;
	typedef itk::ExtractImageFilter< UShortImageType3D, UShortImageType3D > LabelExtractType;
	typedef itk::CastImageFilter< UShortImageType3D, UShortImageType3D > LabelCastType;

	UShortImageType3D::Pointer lab_image = input_labeled_image->GetItkPtr<unsigned short>(0,channel_number);

	int z = input_labeled_image->GetImageInfo()->numZSlices;

	LabelExtractType::Pointer leFilter = LabelExtractType::New();
	UShortImageType3D::RegionType lRegion = lab_image->GetLargestPossibleRegion();
	lRegion.SetSize(2,z);
	leFilter->SetExtractionRegion(lRegion);
	leFilter->SetInput( lab_image );

	LabelCastType::Pointer lFilter = LabelCastType::New();
	lFilter->SetInput( leFilter->GetOutput() );
	try
	{
		lFilter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
	}

	lab_im = lFilter->GetOutput();

	lab_im_set=true;
}

//Update the features in this table whose names match (sets doFeat)
void AssociativeFeatureCalculator::Append(vtkSmartPointer<vtkTable> table)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc;
	if( IsLabeledImageSet() == true )
		assoc = new ftk::NuclearAssociationRules("",0,lab_im);
	else
		assoc = new ftk::NuclearAssociationRules("",0);
	assoc->ReadRulesFromXML(inFilename);
	assoc->PrintSelf();
	assoc->Compute();

	//Init the table (headers):
	for (int i=0; i < assoc->GetNumofAssocRules(); ++i)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (fPrefix+assoc->GetAssociationRules().at(i).GetRuleName()).c_str() );
		column->SetNumberOfValues( table->GetNumberOfRows() );
		table->AddColumn(column);
	}

	//Now update the table:
	std::vector<unsigned short> labels = assoc->GetLabels();
	float** vals = assoc->GetAssocFeaturesList();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		unsigned short id = labels.at(i);
		if(id == 0) continue;

		int row = -1;
		for(int r=0; r<table->GetNumberOfRows(); ++r)
		{
			if( table->GetValue(r,0) == id )
			{
				row = r;
				break;
			}
		}

		if(row == -1) continue;

		for (int f=0; f<assoc->GetNumofAssocRules(); ++f)
		{
			table->SetValueByName(row,(fPrefix+assoc->GetAssociationRules().at(f).GetRuleName()).c_str(), vtkVariant(vals[f][i]));
		}
	}
}

//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************

NuclearAssociationRules::NuclearAssociationRules(std::string AssocFName, int numOfRules):ObjectAssociation(AssocFName, numOfRules)
{
	labImage= NULL;
	labGeometryFilter= NULL;	
	x_Size = y_Size = z_Size = 0;
	imDim=3;	
	objectType = "Nucleus";
	lab_im_set = false;
}

NuclearAssociationRules::NuclearAssociationRules(std::string AssocFName, int numOfRules,LabImageType::Pointer lab_im):ObjectAssociation(AssocFName, numOfRules, lab_im){
	labImage= lab_im;
	labGeometryFilter= NULL;	
	x_Size = y_Size = z_Size = 0;
	imDim=3;	
	objectType = "Nucleus";
	lab_im_set = true;
}

/* This is the main function for computing associative features */
void NuclearAssociationRules::Compute()
{
	std::cout<<"Starting Associative Features Computation\n";

	//1. read the label image
	if( lab_im_set == false ){
		ReaderType::Pointer reader = ReaderType::New();
		std::string fname = GetSegImgName();
		reader->SetFileName(fname);
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << err << std::endl;
			return;
		}
		labImage = LabImageType::New();
		labImage = reader->GetOutput();
	}

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
	//get the number of objects	(Most probably it gets the maximum object label, so be careful!) 
	//numOfLabels  = labGeometryFilter->GetNumberOfLabels();//GetNumberOfObjects();
	//numOfLabels--; //because the first one is the background

	//Get the list of labels
	labelsList = labGeometryFilter->GetLabels();
	numOfLabels = (int)labelsList.size();
	//int maxLable = labelsList[numOfLabels-1];

	//allocate memory for the features list
	assocMeasurementsList = new float*[GetNumofAssocRules()];
	
	//3. then, for each type of the associations get the requested reigion based on the type and the value of the inside and outside distances.
	for(int i=0; i<GetNumofAssocRules(); i++)
	{
		assocMeasurementsList[i] = new float[numOfLabels];
		//read the ith target image (the image from which we need to compute the ith assoc. rule
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName(assocRulesList[i].GetTargetFileNmae());
		reader2->Update();
		if( assocRulesList[i].IsUseBackgroundSubtraction() ){
			if( assocRulesList[i].IsUseMultiLevelThresholding() )
				if( assocRulesList[i].GetNumberOfThresholds()>=assocRulesList[i].GetNumberIncludedInForeground() )
					thresh=returnthresh( reader2->GetOutput(), assocRulesList[i].GetNumberOfThresholds(), assocRulesList[i].GetNumberIncludedInForeground() );
				else
					thresh=returnthresh( reader2->GetOutput(), assocRulesList[i].GetNumberOfThresholds(), assocRulesList[i].GetNumberOfThresholds() );
			else
				thresh=returnthresh( reader2->GetOutput(), 1, 1 );
		}
		else
			thresh = 0;
		//****************************************************Add rule for ECs
		//cout<<"Computing Features For Association Rule "<<i+1<<": ";
		for(int j=0; j<numOfLabels; j++)
		{
			//cout<<j+1;
			int lbl = labelsList[j];
			if(lbl == 0) continue;
			cout<<"\rComputing Features For Association Rule "<<i+1<<": "<<j<<"/"<<numOfLabels-1;
			assocMeasurementsList[i][j] = ComputeOneAssocMeasurement(reader2->GetOutput(), i, lbl);						
		}		
		std::cout<<"\tdone"<<std::endl;
	}	
	
	//Flag invalid objects
	//allocate memory for the invalid objects list
	/*invalidObjects = new unsigned short[numOfLabels];
	for(int j=0; j<numOfLabels; j++)
	{		
		unsigned short all_neg = 1;
		for(int i=0; i<GetNumofAssocRules(); i++)
		{
			if(assocMeasurementsList[i][j] != -1.0)
				all_neg = 0;
		}
		invalidObjects[j] = all_neg;
	}*/
}

/* use this function to compute one associative measurement */
float NuclearAssociationRules::ComputeOneAssocMeasurement(itk::SmartPointer<TargImageType> trgIm, int ruleID, int objID)
{	
	//fisrt, get the bounding box around the object
	//The bounding box is defined by the area around the object such that the object+outerdistance are included	
	int imBounds[6] = {0,x_Size,0,y_Size,0,z_Size};
	LabelGeometryType::BoundingBoxType bbox = labGeometryFilter->GetBoundingBox(objID);
	//make sure the object exists
	int valid = 0;
	for(int dim=0; dim < imDim*2; ++dim)
	{
		if(bbox[dim] > 0)
			valid = 1;
	}
	if(valid == 0)
		return -1;

	std::vector< int > retBbox(0);
	int dist, val;
	dist = assocRulesList[ruleID].GetOutDistance();
	for(int dim=0; dim < imDim*2; ++dim)
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
	if(imDim == 3)
		size[2]  = retBbox[5]-retBbox[4]+1;  // size along Z
	else
		size[2] = 1;
  
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
	if(imDim == 3)
		start2[2] =   retBbox[4];  // first index on Z   
	else
		start2[2] = 0;
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
	int counter = 0;
	while ( ! iterator1.IsAtEnd())
	{		
		int V = iterator1.Get();
		if(V == objID)
		{
			iterator2.Set(255.0);
			counter++; //just for debuging purposes,, will be removed
		}
		else
			iterator2.Set(0.0);
		++iterator1;
		++iterator2;		
	}	
	//Let's try this (for debugging): save the binary mask
	/*	typedef itk::Image< unsigned short, 3 > OutputImageType;
		typedef itk::CastImageFilter< DistImageType, OutputImageType > CastType; 
		CastType::Pointer cast = CastType::New(); 
		cast->SetInput( subSegImg ); 
		cast->Update(); 

		typedef itk::ImageFileWriter< OutputImageType > WriterType; 
		WriterType::Pointer writer = WriterType::New( ); 
		writer->SetInput( cast->GetOutput() ); 
		writer->SetFileName( "c:/bin_mask.tif" ); 
		writer->Update(); */

	//Compute the distance transform in the sub-segmentation image region	
	DTFilter::Pointer dt_obj= DTFilter::New() ;
	dt_obj->SetInput(subSegImg) ;
	//dt_obj->SetInsideValue(255.0);
	//dt_obj->SetOutsideValue(0.0);	
	try{
		dt_obj->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cout << "Error in Distance Transform: " << err << std::endl; 
		return 0;	
	}

	//Let's try this (for debugging): save the distance map
		////typedef itk::Image< unsigned short, 3 > OutputImageType;
		////typedef itk::CastImageFilter< DistImageType, OutputImageType > CastType; 
		////CastType::Pointer cast = CastType::New(); 
		//cast->SetInput( dt_obj->GetOutput() ); 
		//cast->Update(); 

		////typedef itk::ImageFileWriter< OutputImageType > WriterType; 
		////WriterType::Pointer writer = WriterType::New( ); 
		//writer->SetInput( cast->GetOutput() ); 
		//writer->SetFileName( "c:/dist_map.tif" ); 
		//writer->Update(); 
	
	//now, mask out all the pixels (in the sub-seg image) that are not in the region of interest as defined by the association rule and get the intensities of the needed pixels from the target image. The intensities are saved into an std vector
	IteratorType2 iterator3(dt_obj->GetOutput(), dt_obj->GetOutput()->GetRequestedRegion());
	IteratorType iterator4(trgIm, trgIm->GetRequestedRegion());
	std::vector<int> trgInt;
	int counter_in = 0;
	int counter_at = 0;
	int counter_ot = 0;
	while ( ! iterator3.IsAtEnd())
	{
		int V = (int)iterator3.Get();
		if(V<0)
			counter_in++;
		if(V==0)
			counter_at++;
		if(V>0)
			counter_ot++;
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
  //we will never go here, just silencing a compiler warning
  return 0.0;
}

/* the next functions will compute the measurements for the different cases */
//find maximum
float NuclearAssociationRules::FindMax(std::vector<int> LST)
{
	float mx = LST[0];
	for(unsigned int i=1; i<LST.size(); i++)
	{
		if(mx<LST[i])
			mx = LST[i];
	}
	return mx;
}

//find minimum
float NuclearAssociationRules::FindMin(std::vector<int> LST)
{
	float mn = LST[0];
	for(unsigned int i=1; i<LST.size(); i++)
	{
		if(mn>LST[i])
			mn = LST[i];
	}
	return mn;
}

//compute total
float NuclearAssociationRules::ComputeTotal(std::vector<int> LST)
{
	float tl = 0;
	for(unsigned int i=0; i<LST.size(); i++)
	{
		if( LST[i] >= thresh )
			tl += LST[i];
	}
	return tl;
}

//compute average
float NuclearAssociationRules::ComputeAverage(std::vector<int> LST)
{
	float av = 0;
	float lst_sz = (float)LST.size();
	for(unsigned int i=0; i<LST.size(); i++)
	{
		if( LST[i] >= thresh )
			av += LST[i];
		else
			--lst_sz;
	}
	if( lst_sz )
		av/=lst_sz;
	return av;
}

} //end namespace ftk