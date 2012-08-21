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
//#include "itkImageFileWriter.h"

namespace ftk 
{	

AssociativeFeatureCalculator::AssociativeFeatureCalculator()
{
	inFilename = "";
	fPrefix = "";
	inputs_set = false;	
}

AssociativeFeatureCalculator::~AssociativeFeatureCalculator()
{
}

void AssociativeFeatureCalculator::SetInputFile(std::string filename)
{
	inFilename = filename;
}

void AssociativeFeatureCalculator::SetInputs(ftk::Image::Pointer inp_labeled_image, int inp_channel_number, ftk::Image::Pointer seg_labeled_image, int seg_channel_number, ftk::AssociationRule *associationrule){
	input_association = associationrule;
	if( seg_channel_number == -1 || inp_channel_number == -1 )
		return;

	typedef itk::ExtractImageFilter< LabImageType, LabImageType > LabelExtractType;
	typedef itk::ExtractImageFilter< TargImageType, TargImageType > InputExtractType;
	typedef itk::CastImageFilter< LabImageType, LabImageType > LabelCastType;
	typedef itk::CastImageFilter< TargImageType, TargImageType > InputCastType;

	lab_im = seg_labeled_image->GetItkPtr<unsigned short>(0,seg_channel_number,ftk::Image::DEEP_COPY);
	inp_im = inp_labeled_image->GetItkPtr<unsigned short>(0,inp_channel_number,ftk::Image::DEEP_COPY);

	/*int z = seg_labeled_image->GetImageInfo()->numZSlices;
	int z1 = inp_labeled_image->GetImageInfo()->numZSlices;

	LabelExtractType::Pointer leFilter = LabelExtractType::New();
	InputExtractType::Pointer leFilter1 = InputExtractType::New();
	LabImageType::RegionType lRegion = lab_image->GetLargestPossibleRegion();
	TargImageType::RegionType lRegion1 = inp_image->GetLargestPossibleRegion();
	lRegion.SetSize(2,z);
	lRegion1.SetSize(2,z1);
	leFilter->SetExtractionRegion(lRegion);
	leFilter1->SetExtractionRegion(lRegion1);
	leFilter->SetInput( lab_image );
	leFilter1->SetInput( inp_image );

	LabelCastType::Pointer lFilter = LabelCastType::New();
	InputCastType::Pointer lFilter1 = InputCastType::New();
	lFilter->SetInput( leFilter->GetOutput() );
	lFilter1->SetInput( leFilter1->GetOutput() );
	try
	{
		lFilter->Update();
		lFilter1->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
	}

	lab_im = lFilter->GetOutput();
	inp_im = lFilter1->GetOutput();*/

	/*typedef itk::ImageFileWriter< TargImageType > writertype;
	writertype::Pointer writer = writertype::New();
	writer->SetInput( inp_im );
	writer->SetFileName( "cast_image.tif" );
	writer->Update();*/

	inputs_set = true;
}

void AssociativeFeatureCalculator::SetFeaturePrefix(std::string prefix)
{
	fPrefix = prefix;
}

//Compute features turned ON in doFeat and puts them in a new table
vtkSmartPointer<vtkTable> AssociativeFeatureCalculator::Compute(void)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc;
	if( inputs_set ){
		assoc = new ftk::NuclearAssociationRules("",0,lab_im, inp_im);
		assoc->AddAssociation( input_association->GetRuleName(), "", input_association->GetOutDistance(), input_association->GetInDistance(),	input_association->IsUseWholeObject(), input_association->IsUseBackgroundSubtraction(), input_association->IsUseMultiLevelThresholding(), input_association->GetNumberOfThresholds(), input_association->GetNumberIncludedInForeground(), input_association->GetAssocType(), input_association->get_path() );
	}
	else{
		assoc = new ftk::NuclearAssociationRules("",0);
		assoc->ReadRulesFromXML(inFilename);
	}
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
	delete assoc;
	return table;
}

//Update the features in this table whose names match (sets doFeat)
void AssociativeFeatureCalculator::Update(vtkSmartPointer<vtkTable> table)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc;
	if( inputs_set ){
		assoc = new ftk::NuclearAssociationRules("",0,lab_im, inp_im);
		assoc->AddAssociation( input_association->GetRuleName(), "", input_association->GetOutDistance(), input_association->GetInDistance(),	input_association->IsUseWholeObject(), input_association->IsUseBackgroundSubtraction(), input_association->IsUseMultiLevelThresholding(), input_association->GetNumberOfThresholds(), input_association->GetNumberIncludedInForeground(), input_association->GetAssocType(), input_association->get_path() );
	}
	else{
		assoc = new ftk::NuclearAssociationRules("",0);
		assoc->ReadRulesFromXML(inFilename);
	}
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
	delete assoc;
}

//Update the features in this table whose names match (sets doFeat)
void AssociativeFeatureCalculator::Append(vtkSmartPointer<vtkTable> table)
{
	//Compute features:
	ftk::NuclearAssociationRules *assoc;
	if( inputs_set ){
		assoc = new ftk::NuclearAssociationRules("",0,lab_im, inp_im);
		assoc->AddAssociation( input_association->GetRuleName(), "", input_association->GetOutDistance(), input_association->GetInDistance(),	input_association->IsUseWholeObject(), input_association->IsUseBackgroundSubtraction(), input_association->IsUseMultiLevelThresholding(), input_association->GetNumberOfThresholds(), input_association->GetNumberIncludedInForeground(), input_association->GetAssocType(), input_association->get_path() );
	}
	else{
		assoc = new ftk::NuclearAssociationRules("",0);
		assoc->ReadRulesFromXML(inFilename);
	}
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
	//#pragma omp parallel for num_threads(4)
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
	delete assoc;
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
	inputs_set = false;
	num_rois = 8;
	assocMeasurementsList = NULL;
}

NuclearAssociationRules::NuclearAssociationRules(std::string AssocFName, int numOfRules, LabImageType::Pointer lImage, TargImageType::Pointer iImage):ObjectAssociation(AssocFName, numOfRules){
	labImage= lImage;
	inpImage= iImage;
	labGeometryFilter= NULL;	
	x_Size = y_Size = z_Size = 0;
	imDim=3;	
	objectType = "Nucleus";
	inputs_set = true;
	num_rois = 8;
	assocMeasurementsList = NULL;
}

NuclearAssociationRules::~NuclearAssociationRules()
{
	if( assocMeasurementsList!= NULL )
	{
		for(int i=0; i<GetNumofAssocRules(); i++)
			delete []assocMeasurementsList[i];
		delete [] assocMeasurementsList;
	}
}

/* This is the main function for computing associative features */
void NuclearAssociationRules::Compute(vtkSmartPointer<vtkTable> tbl)
{
	std::cout<<"Starting Associative Features Computation\n";

	/*f(assocRulesList[ruleID].GetAssocType() == ASSOC_DIST_DEVICE)
	{*/
	for(int i=0; i<GetNumofAssocRules(); i++)
	{
		if(assocRulesList[i].GetAssocType() == ASSOC_DIST_OBJECT)
		{
			std::string object_file = assocRulesList[i].GetTargetFileNmae();
			VolumeOfInterest * VOIType = new VolumeOfInterest();
			VOIType->ReadVTPVOI(object_file);
			assocMeasurementsList[i] = VOIType->CalculateCentroidDistanceToVOI(tbl);			
		}
	}

	//1. read the label image
	if( inputs_set == false ){
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
	for( int i=0; i<(int)labelsList.size(); ++i )
		if( labelsList[i] == 0 ){
			--numOfLabels;
			labelsList.erase(labelsList.begin()+i);
		}
	//int maxLable = labelsList[numOfLabels-1];

	//allocate memory for the features list
	assocMeasurementsList = new float*[GetNumofAssocRules()];
	
	//3. then, for each type of the associations get the requested reigion based on the type and the value of the inside and outside distances.
	for(int i=0; i<GetNumofAssocRules(); i++)
	{
		if(assocRulesList[i].GetAssocType() == ASSOC_DIST_OBJECT )
			continue;

		assocMeasurementsList[i] = new float[numOfLabels];
		//read the ith target image (the image from which we need to compute the ith assoc. rule
		if( inputs_set == false ){
			ReaderType::Pointer reader2 = ReaderType::New();
			reader2->SetFileName(assocRulesList[i].GetTargetFileNmae());
			reader2->Update();

			//Scale so that there is always uniform scaling of association computation when comparing across uchar and ushort images
			typedef itk::RescaleIntensityImageFilter< LabImageType, LabImageType > RescaleUsUsType;
			RescaleUsUsType::Pointer rescaleususfilter = RescaleUsUsType::New();
			rescaleususfilter->SetInput( reader2->GetOutput() );
			rescaleususfilter->SetOutputMaximum( itk::NumericTraits<unsigned char>::max() );
			rescaleususfilter->SetOutputMinimum( 0 );
			rescaleususfilter->Update();

			inpImage = rescaleususfilter->GetOutput();
		}

		if( assocRulesList[i].IsUseBackgroundSubtraction() ){
			if( assocRulesList[i].IsUseMultiLevelThresholding() )
				if( assocRulesList[i].GetNumberOfThresholds()>=assocRulesList[i].GetNumberIncludedInForeground() )
					thresh=returnthresh( inpImage, assocRulesList[i].GetNumberOfThresholds(), assocRulesList[i].GetNumberIncludedInForeground() );
				else
					thresh=returnthresh( inpImage, assocRulesList[i].GetNumberOfThresholds(), assocRulesList[i].GetNumberOfThresholds() );
			else
				thresh=returnthresh( inpImage, 1, 1 );
			//Write Binary Mask
			/*std::string out_filename;
			out_filename = assocRulesList[i].get_path()+assocRulesList[i].GetRuleName();
			if( assocRulesList[i].GetAssocType() == ASSOC_SURROUNDEDNESS )
				out_filename = out_filename + "binary_surroundedness.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_MIN )
				out_filename = out_filename + "binary_min.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_MAX )
				out_filename = out_filename + "binary_max.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_TOTAL )
				out_filename = out_filename + "binary_total.tif";
			if( assocRulesList[i].GetAssocType() == ASSOC_AVERAGE )
				out_filename = out_filename + "binary_average.tif";
			BinaryThresholdType::Pointer threshfilt = BinaryThresholdType::New();
			threshfilt->SetInput( inpImage );
			threshfilt->SetLowerThreshold(thresh);
			threshfilt->SetUpperThreshold(USHRT_MAX);
			threshfilt->SetInsideValue( USHRT_MAX );
			threshfilt->SetOutsideValue( 0 );
			threshfilt->Update();
			WriterType::Pointer writer1 = WriterType::New();
			writer1->SetInput( threshfilt->GetOutput() );
			writer1->SetFileName( out_filename.c_str() );
			writer1->Update();*/
		}
		else
			thresh = 0;
		//cout<<"Computing Features For Association Rule "<<i+1<<": ";
		if( assocRulesList[i].GetAssocType() == ASSOC_SURROUNDEDNESS ){
		std::vector<float> ec_feat_vals;
			if( numOfLabels )
				 ec_feat_vals = compute_ec_features( inpImage, labImage, num_rois, thresh, assocRulesList[i].GetOutDistance(), assocRulesList[i].GetInDistance()  );
			for(int j=0; j<numOfLabels; ++j)
				assocMeasurementsList[i][j] = ec_feat_vals[j];
		} else {

			int counterLabels = 0;
//#ifdef _OPENMP
//			std::cout << std::endl << "ASSOCIATED FEATURES WILL BE COMPUTED USING OPENMP";
//		omp_set_nested(1);
//#endif
//#pragma omp parallel for 
			for(int j=0; j<numOfLabels; ++j)
			{
				//cout<<j+1;
				int lbl = labelsList[j];
				if(lbl == 0) continue;
				//cout<<"\rComputing Features For Association Rule "<<i+1<<": "<<j<<"/"<<numOfLabels-1;
				assocMeasurementsList[i][j] = ComputeOneAssocMeasurement(inpImage, i, lbl);	
//#pragma omp critical
				//{
				//	counterLabels++;
				//	std::cout << std::endl << "Fea, Rule " << i+1 << ": " << counterLabels << " of " << numOfLabels-1 << " DONE";
				//}
			}
//#ifdef _OPENMP
			//omp_set_nested(0);
//#endif
		}
		//std::cout<<"\tdone"<<std::endl;
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
	double av = 0;
	double lst_sz = (double)LST.size();
	for(unsigned int i=0; i<LST.size(); i++)
	{
		if( LST[i] >= thresh )
			av += LST[i];
		else
			--lst_sz;
	}
	if( lst_sz > 0 && av > 0 ){
		av/=lst_sz;
		return (float)av;
	} else return 0;
}

} //end namespace ftk
