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

bool IntrinsicFeatureCalculator::SetInputImages(ftk::Image::Pointer intImg, ftk::Image::Pointer labImg, int intChannel, int labChannel)
{
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
	labFilter->SetLevel( needLevel() );
	if( needHistogram() )
		labFilter->ComputeHistogramOn();
	if( needTextures() )
		labFilter->ComputeTexturesOn();
	labFilter->Update();

	//Init the table (headers):
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "ID" );
	table->AddColumn(column);
	for (int i=0; i < IntrinsicFeatures::N; ++i)
	{
		if(doFeat[i])
		{
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (fPrefix+IntrinsicFeatures::Info[i].name).c_str() );
			table->AddColumn(column);
		}
	}

	//Now populate the table:
	std::vector< FeatureCalcType::LabelPixelType > labels = labFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
	{
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) continue;
		if(useIDs)
			if(IDs.find(id) == IDs.end()) continue;	//Don't care about this id, so skip it

		ftk::IntrinsicFeatures * features = labFilter->GetFeatures(id);
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant(id) );
		for (int i=0; i<IntrinsicFeatures::N; ++i)
		{
			if(doFeat[i])
				row->InsertNextValue( vtkVariant(features->ScalarFeatures[i]) );
		}
		table->InsertNextRow(row);
	}
	return table;
}

//Update the features in this table whose names match (sets doFeat)
//Will create new rows if necessary:
void IntrinsicFeatureCalculator::Update(vtkSmartPointer<vtkTable> table, std::map<int, ftk::Object::Point> * cc, std::map<int, ftk::Object::Box> * bbox)
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
		labFilter->SetLevel(1);
		labFilter->ComputeHistogramOff();
		labFilter->ComputeTexturesOff();
	}

	labFilter->Update();

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
				nrow->SetNumberOfValues( table->GetNumberOfColumns() );
				table->InsertNextRow(nrow);
				row = table->GetNumberOfRows() - 1;
				table->SetValue(row, 0, vtkVariant(id));
			}

			//Update table:
			for (int f=0; f<IntrinsicFeatures::N; ++f)
			{
				if(doFeat[f])
					table->SetValueByName(row,(fPrefix+IntrinsicFeatures::Info[f].name).c_str(), vtkVariant(features->ScalarFeatures[f]));
			}
		}

		//Update centroids:
		if(cc)
		{
			Object::Point c;
			c.x = (int)features->Centroid[0];
			c.y = (int)features->Centroid[1];
			c.z = (int)features->Centroid[2];
			c.t = 0;
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
			b.min.z = (int)features->BoundingBox[4];
			b.max.z = (int)features->BoundingBox[5];
			b.min.t = 0;
			b.max.t = 0;
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
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel), regionIndex, regionSize );
	}
	else
	{
		labFilter->SetImageInputs( intensityImage->GetItkPtr<IPixelT>(0,intensityChannel), labelImage->GetItkPtr<LPixelT>(0,labelChannel) );
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
			nrow->SetNumberOfValues( table->GetNumberOfColumns() );
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

}  // end namespace ftk