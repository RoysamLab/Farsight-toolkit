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
#include "ftkLabelImageToFeatures.h"

namespace ftk
{

IntrinsicFeatureCalculator::IntrinsicFeatureCalculator()
{
	intensityImage = NULL;
	intensityChannel = 0;
	labelImage = NULL;
	labelChannel = 0;
	fPrefix = "";
	SetFeaturesOn();
	this->SetFeatureOn(IntrinsicFeatures::SKEW, false);		//Default is to have these 4 off because histogram computation too costly
	this->SetFeatureOn(IntrinsicFeatures::ENERGY, false);
	this->SetFeatureOn(IntrinsicFeatures::ENTROPY, false);
	this->SetFeatureOn(IntrinsicFeatures::MEDIAN, false);
	useRegion = false;
	useIDs = false;
	IDs.clear();
}

bool IntrinsicFeatureCalculator::SetInputImages(ftk::Image::Pointer intImg, ftk::Image::Pointer labImg, int intChannel, int labChannel, bool CytoImage)
{
	cyto_image = CytoImage;
	if(intImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
		return false;
	if(labImg->GetImageInfo()->dataType != itk::ImageIOBase::USHORT)
		return false;

	if(intImg->GetImageInfo()->numChannels <= intChannel)
		return false;
	if(labImg->GetImageInfo()->numChannels <= labChannel)
		return false;

	intensityImage = intImg;
	intensityChannel = intChannel;
	labelImage = labImg;
	labelChannel = labChannel;
	return true;
}

//Turn on all features:
void IntrinsicFeatureCalculator::SetFeaturesOn(void)
{
	for(int i=0; i<IntrinsicFeatures::N; ++i)
	{
		doFeat[i] = true;
	}
}

void IntrinsicFeatureCalculator::SetFeaturesOn(std::set<int> onFeats)
{
	std::set<int>::iterator it;
	for( int i=0; i<IntrinsicFeatures::N; ++i )
	{
		if( onFeats.find(i) != onFeats.end() )
			doFeat[i] = true;
		else
			doFeat[i] = false;
	}
}

void IntrinsicFeatureCalculator::SetFeatureOn(int feat, bool v)
{
	if(feat<IntrinsicFeatures::N)
		doFeat[feat] = v;
}
	
void IntrinsicFeatureCalculator::SetFeaturePrefix(std::string prefix)
{
	fPrefix = prefix;
}

void IntrinsicFeatureCalculator::SetRegion(int x1, int y1, int z1, int x2, int y2, int z2)
{
	useRegion = true;

	regionIndex[0] = x1;
	regionIndex[1] = y1;
	regionIndex[2] = z1;

	regionSize[0] = x2 - x1 + 1;
	regionSize[1] = y2 - y1 + 1;
	regionSize[2] = z2 - z1 + 1;
}

void IntrinsicFeatureCalculator::SetIDs(std::set<LPixelT> ids)
{
	useIDs = true;
	IDs = ids;
}

int IntrinsicFeatureCalculator::getMaxFeatureTurnedOn(void)
{
	int max = 0;
	for(int i=0; i<IntrinsicFeatures::N; ++i)
	{
		if(doFeat[i] && i > max) 
			max = i;
	}
	return max;
}

bool IntrinsicFeatureCalculator::needTextures(void)
{
	if( getMaxFeatureTurnedOn() >= IntrinsicFeatures::T_ENERGY )
		return true;
	else
		return false;
}

bool IntrinsicFeatureCalculator::needHistogram(void)
{
	if( doFeat[IntrinsicFeatures::MEDIAN] || doFeat[IntrinsicFeatures::SKEW] || doFeat[IntrinsicFeatures::ENERGY] || doFeat[IntrinsicFeatures::ENTROPY] )
		return true;
	else
		return false;
}
int IntrinsicFeatureCalculator::needLevel(void)
{
	int max = getMaxFeatureTurnedOn();
	if(max <= IntrinsicFeatures::BBOX_VOLUME)
		return 1;
	else if(max <= IntrinsicFeatures::VARIANCE)
		return 2;
	else
		return 3;
}

//Compute features turned ON in doFeat and puts them in a new table
vtkSmartPointer<vtkTable> IntrinsicFeatureCalculator::Compute(void)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

	if(!intensityImage || !labelImage)
	{
		return table;
	}

	//Compute features:
	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	if(useRegion)
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel), regionIndex, regionSize );
	}
	else
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel) );
	}
	//labFilter->SetLevel( needLevel() );
	labFilter->SetLevel( 3); //////////////////////////modified by Yanbin from 3 to 2;
	labFilter->ComputeSurfaceOn();
	if( needHistogram() )
		labFilter->ComputeHistogramOn();
	if( needTextures() )
		labFilter->ComputeTexturesOn();
	labFilter->Update();

	//Init the table (headers):
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	table->AddColumn(column);
	
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_x" );
	table->AddColumn(column);

	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_y" );
	table->AddColumn(column);

	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_z" );
	table->AddColumn(column);

//	column = vtkSmartPointer<vtkDoubleArray>::New();
//	column->SetName( "num_z_slices" );
//	table->AddColumn(column);
//
	for (int i=0; i < IntrinsicFeatures::N; ++i)
	{
		if(doFeat[i])
		{
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (fPrefix+IntrinsicFeatures::Info[i].name).c_str() );
			table->AddColumn(column);
		}
	}

	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();

	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it
		std::vector< std::vector<double> > testVec = labFilter->GetZernikeMoments(id);
		for(int i=0; i<(int)testVec.size(); ++i)
		{
			std::stringstream ss1;
			ss1 << i;
			for(int j=0; j<(int)testVec.at(i).size(); ++j)
			{
				column = vtkSmartPointer<vtkDoubleArray>::New();
				std::stringstream ss2;
				ss2 << ((i%2)+(2*j));
				column->SetName( ("Zern_"+ss1.str()+"_"+ss2.str()).c_str() );
				table->AddColumn(column);
				std::cout<<"computing zernike"<<std::endl;
			}
		}
	}

	//Now populate the table:
	//std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant(id) );
		row->InsertNextValue( vtkVariant((int)features->Centroid[0]) );
		row->InsertNextValue( vtkVariant((int)features->Centroid[1]) );
		if( intensityImage->GetImageInfo()->numZSlices > 2 )
			row->InsertNextValue( vtkVariant((int)features->Centroid[2]) );
		else
			row->InsertNextValue( vtkVariant(0));
//		row->InsertNextValue( vtkVariant((int)intensityImage->GetImageInfo()->numZSlices) );
		for (int i=0; i<IntrinsicFeatures::N; ++i)
		{
			if(doFeat[i])
				row->InsertNextValue( vtkVariant(features->ScalarFeatures[i]) );
		}
		std::vector< std::vector<double> > zernVec = labFilter->GetZernikeMoments(id);
		for(int i=0; i<(int)zernVec.size(); ++i)
		{
			for(int j=0; j<(int)zernVec.at(i).size(); ++j)
			{
				row->InsertNextValue( vtkVariant(zernVec.at(i).at(j)) );
			}
		}
		table->InsertNextRow(row);
	}
	return table;
}

//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
void IntrinsicFeatureCalculator::GetObjectCentroids(vtkSmartPointer<vtkTable> table, int time)
{
	if(!intensityImage || !labelImage)
	{
		return;
	}

	//Compute features:
	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	if(useRegion)
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(time,intensityChannel), labelImage->GetItkPtr<LPixelT>(time,labelChannel), regionIndex, regionSize );
	}
	else
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(time,intensityChannel), labelImage->GetItkPtr<LPixelT>(time,labelChannel) );
	}
	labFilter->Update();

	//Init the table (headers):
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_x" );
	column->SetNumberOfValues(table->GetNumberOfRows());
	table->AddColumn(column);

	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_y" );
	column->SetNumberOfValues(table->GetNumberOfRows());
	table->AddColumn(column);

	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "centroid_z" );
	column->SetNumberOfValues(table->GetNumberOfRows());
	table->AddColumn(column);
	
	for (int i=0; i<(int)table->GetNumberOfRows(); ++i)
	{
		FeatureCalcType::LabelPixelType id = table->GetValue(i,0).ToInt();
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
		table->SetValueByName( i, "centroid_x", vtkVariant((int)features->Centroid[0]) );
		table->SetValueByName( i, "centroid_y", vtkVariant((int)features->Centroid[1]) );
		if( labelImage->GetImageInfo()->numZSlices > 2 )
			table->SetValueByName( i, "centroid_z", vtkVariant((int)features->Centroid[2]) );
		else
			table->SetValueByName( i, "centroid_z", vtkVariant(0) );

		
	}
	return;
}

//**************************************************************************************
//**************************************************************************************
//**************************************************************************************

//Update the features in this table whose names match (sets doFeat)
//Will create new rows if necessary:
void IntrinsicFeatureCalculator::Update(vtkSmartPointer<vtkTable> table, std::map<int, ftk::Object::Point> * cc, std::map<int, ftk::Object::Box> * bbox, vtkSmartPointer<vtkTable> NucAdjTable, int currtime)
{
	if(!intensityImage || !labelImage)
		return;

	if(table)
	{
		//Determine needed features:
		for(int i=0; i<IntrinsicFeatures::N; ++i)
		{
			std::string name = fPrefix+IntrinsicFeatures::Info[i].name;
			vtkAbstractArray * arr = table->GetColumnByName( name.c_str() );
			if(arr == NULL || arr == 0)
				doFeat[i] = false;
			else
				doFeat[i] = true;
		}
	}

	//Compute features:
	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	if(useRegion)
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel), regionIndex, regionSize );
	}
	else
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel) );
	}
	
	if(table)
	{
		labFilter->SetLevel( needLevel() );
		if( needHistogram() )
			labFilter->ComputeHistogramOn();
		if( needTextures() )
			labFilter->ComputeTexturesOn();
	}
	else
	{
		labFilter->SetLevel(3);  //modified by yanbin from 2 to 3 ;
		labFilter->ComputeHistogramOff();
		labFilter->ComputeTexturesOff();
	}

	labFilter->Update();

	//Update the Nuclear Adjacency Table
	if(NucAdjTable)
	{
		FeatureCalcType::Pointer AdjFilter = FeatureCalcType::New();
		AdjFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel) );
		AdjFilter->GetAdjacency();
		std::set<unsigned short>::iterator it;
		for(it=IDs.begin(); it!=IDs.end(); ++it)
		{
			std::vector<unsigned short> conIDs = AdjFilter->GetContactNeighbors(*it);
			if(conIDs.size()>0)
			{
				for(unsigned int i=0 ; i<conIDs.size() ; ++i)
				{
					if(conIDs[i] == 0) continue;
					vtkSmartPointer<vtkVariantArray> nextrow = vtkSmartPointer<vtkVariantArray>::New();
					nextrow->InsertNextValue(*it);
					nextrow->InsertNextValue(conIDs[i]);
					NucAdjTable->InsertNextRow(nextrow);
				}
			}
		}
	}
	//Now update the table:
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);

		if(table)
		{
			int row = -1;
			for(int r=0; r<table->GetNumberOfRows(); ++r)
			{
				if( table->GetValue(r,0) == id )
				{
					row = r;
					break;
				}
			}

			//Must Create a new row:
			if(row == -1)
			{
				vtkSmartPointer<vtkVariantArray> nrow = vtkSmartPointer<vtkVariantArray>::New();
				for( unsigned ii=0; ii<table->GetNumberOfColumns(); ++ii )
					nrow->InsertNextValue( vtkVariant(0.0) );
				table->InsertNextRow(nrow);
				row = table->GetNumberOfRows() - 1;
				table->SetValue(row, 0, vtkVariant(id));
			}

			//Update table:
			table->SetValueByName(row,"centroid_x", vtkVariant((int)features->Centroid[0]));
			table->SetValueByName(row,"centroid_y", vtkVariant((int)features->Centroid[1]));
			if( labelImage->GetImageInfo()->numZSlices > 2 )
				table->SetValueByName(row,"centroid_z", vtkVariant((int)features->Centroid[2]));
			else
				table->SetValueByName(row,"centroid_z", vtkVariant(0));

//			table->SetValueByName(row,"num_z_slices", vtkVariant((int)intensityImage->GetImageInfo()->numZSlices));
			int col_count = 0;
			for (int f=0; f<IntrinsicFeatures::N; ++f)
			{
				if(doFeat[f])
				{
					table->SetValueByName(row,(fPrefix+IntrinsicFeatures::Info[f].name).c_str(), vtkVariant(features->ScalarFeatures[f]));
					++col_count;
				}
			}
			std::vector< std::vector<double> > zernVec = labFilter->GetZernikeMoments(id);
			for(int i=0; i<(int)zernVec.size(); ++i)
			{
				for(int j=0; j<(int)zernVec.at(i).size(); ++j)
				{
					table->SetValue(row, ++col_count ,vtkVariant(zernVec.at(i).at(j)) );
				}
			}
		}

		//Update centroids:
		if(cc)
		{
			Object::Point c;
			c.x = (int)features->Centroid[0];
			c.y = (int)features->Centroid[1];
			if( intensityImage->GetImageInfo()->numZSlices > 2 )
				c.z = (int)features->Centroid[2];
			else
				c.z = 0;
			c.t = currtime;
			(*cc)[(int)id] = c;
		}

		//Update bounding boxes:
		if(bbox)
		{
			Object::Box b;
			b.min.x = (int)features->BoundingBox[0];
			b.max.x = (int)features->BoundingBox[1];
			b.min.y = (int)features->BoundingBox[2];
			b.max.y = (int)features->BoundingBox[3];
			if( intensityImage->GetImageInfo()->numZSlices > 2 )
			{
				b.min.z = (int)features->BoundingBox[4];
				b.max.z = (int)features->BoundingBox[5];
			}
			else
			{
				b.min.z = 0;
				b.max.z = 0;
			}
			b.min.t = currtime;
			b.max.t = currtime;
			(*bbox)[(int)id] = b;
		}
	}
}

//Update the features in this table whose names match (sets doFeat)
//Will creat new rows if necessary:
void IntrinsicFeatureCalculator::Append(vtkSmartPointer<vtkTable> table)
{
	if(!intensityImage || !labelImage)
		return;

	//Compute features:
	typedef ftk::LabelImageToFeatures< IPixelT, LPixelT, 3 > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	if(useRegion)
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel), regionIndex, regionSize, cyto_image );
	}
	else
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel), cyto_image );
	}
	labFilter->SetLevel( needLevel() );
	if( needHistogram() )
		labFilter->ComputeHistogramOn();
	if( needTextures() )
		labFilter->ComputeTexturesOn();
	labFilter->Update();

	//Add new columns to the table (headers):
	for (int i=0; i < IntrinsicFeatures::N; ++i)
	{
		if(doFeat[i])
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (fPrefix+IntrinsicFeatures::Info[i].name).c_str() );
			column->SetNumberOfValues( table->GetNumberOfRows() );
			table->AddColumn(column);
		}
	}

	//Now update the table:
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		int row = -1;
		for(int r=0; r<table->GetNumberOfRows(); ++r)
		{
			if( table->GetValue(r,0) == id )
			{
				row = r;
				break;
			}
		}

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);

		//Must Create a new row:
		if(row == -1)
		{
			vtkSmartPointer<vtkVariantArray> nrow = vtkSmartPointer<vtkVariantArray>::New();
			for( unsigned ii=0; ii<table->GetNumberOfColumns(); ++ii )
				nrow->InsertNextValue( vtkVariant(0.0) );
			table->InsertNextRow(nrow);
			row = table->GetNumberOfRows() - 1;
			table->SetValue(row, 0, vtkVariant(id));
		}

		for (int f=0; f<IntrinsicFeatures::N; ++f)
		{
			if(doFeat[f])
				table->SetValueByName(row,(fPrefix+IntrinsicFeatures::Info[f].name).c_str(), vtkVariant(features->ScalarFeatures[f]));
		}
	}
}

void IntrinsicFeatureCalculator::GetDistanceToSurfaceMeasures(vtkSmartPointer<vtkTable> table, std::vector< ftk::Object::Point > surfacePoints)
{
	typedef itk::Image< IPixelT, 3 > InputImageType;
	typedef itk::Image< LPixelT, 3 > OutputImageType;
	typedef itk::DanielssonDistanceMapImageFilter< OutputImageType, OutputImageType > DanielssonFilterType;
	typedef itk::RescaleIntensityImageFilter< OutputImageType, OutputImageType > RescalerType;
	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< InputImageType >  InputWriterType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	typedef itk::LineIterator< OutputImageType > LineIteratorType;

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "Dist_To_Surface" );
	column->SetNumberOfValues( table->GetNumberOfRows() );
	table->AddColumn(column);	

	OutputImageType::Pointer im;
	im = OutputImageType::New();
	OutputImageType::PointType origin;
    origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    
    im->SetOrigin( origin );

    OutputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    
	OutputImageType::Pointer temp_image = labelImage->GetItkPtr< LPixelT >(0,0);
	itk::Size<3> im_size = temp_image->GetBufferedRegion().GetSize();
	im_size[2] = 1;
  
    InputImageType::RegionType region;
    region.SetSize( im_size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	for(int p = 1; p < (int)surfacePoints.size(); ++p)
	{
		itk::Index<3> indx,indy;
		indx[0] = surfacePoints[p-1].x;
		indx[1] = surfacePoints[p-1].y;
		indx[2] = 0;
		indy[0] = surfacePoints[p].x;
		indy[1] = surfacePoints[p].y;
		indy[2] = 0;
		LineIteratorType it( im, indx, indy );
		//it.GoToBegin();
		//while(!it.IsAtEnd())
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			it.Set(255);			
		}
	}

	DanielssonFilterType::Pointer danielssonFilter = DanielssonFilterType::New();	
	WriterType::Pointer writer = WriterType::New();	
	danielssonFilter->SetInput( im );
	writer->SetFileName( "DistanceMap.tif" );
	danielssonFilter->InputIsBinaryOn();
	danielssonFilter->Update();
	OutputImageType::Pointer distMap = danielssonFilter->GetOutput();
	writer->SetInput( distMap );
	writer->Update();
	
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		OutputImageType::IndexType indx;
		indx[0] = table->GetValue(row, 1).ToInt();
		indx[1] = table->GetValue(row, 2).ToInt();
		indx[2] = 0;
		int dist = distMap->GetPixel(indx);
		table->SetValueByName(row, "Dist_To_Surface", vtkVariant(dist));
	}



}

}  // end namespace ftk
