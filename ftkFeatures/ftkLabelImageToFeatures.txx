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

#ifndef _ftkLabelImageToFeatures_txx
#define _ftkLabelImageToFeatures_txx
#include "ftkLabelImageToFeatures.h"

#include <itkConstNeighborhoodIterator.h>
#include <itkConstantBoundaryCondition.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkExtractImageFilter.h>
//#include<conio.h>


#include <math.h>

//******************************************************************************
// LabelImageToFeatures.h/cpp is a class similar to an ITK Filter.
// It inherits from itk::LightObject in order to take advantage of smart pointers.
//
// 2 Required Inputs: intensity image and corresponding label image
//
// Uses the itkLabelGeometryFilter by Dirk Padfield, as well as other itk filters:
// itkLabelStatisticsImageFilter, itkGradientMagnitudeImageFilter
// to compute features of label regions in N-D images (N=2,3 common).
// Several other features are calculated based on the information from these itk
// filters.
//*****************************************************************************

//*****************************************************************************
// USAGE:
//
// 1. Create a pointer to the Feature Calculator:
//		LabelFeaturesType::Pointer features = LabelFeaturesType::New();
// 2. Set the inputs:
//		features->SetImageInputs( intensityImage, labelImage );
// 3. Set Desired Level and optional computations (see below):
//		features->SetLevel(2);
//		features->ComputeHistogramOn();
// 4. Update:		
//		features->Update();
// 5. Use these functions to retrieve information:
//		GetLabels, GetPercentSharedBoundary, GetContactNeighbors, GetFeatures
//
// OPTIONS:
// The Filter is broken into 3 levels of computation with these features:
// LEVEL 1:
//	volume, integrated intensity, centroid, weighted centroid, axis length,
//  eccentricity, elongation, orientation, bounding box, bounding box volume
// LEVEL 2:
//  LEVEL 1 + sum, mean, min, max, sigma, variance
// LEVEL 3:
//  LEVEL 2 + surface gradient, interior gradient, surface intensity, interior intensity,
//  intensity ratio, average distance, radius variation, distance variation, 
//  surface area, shape, percent shared boundary.
//
// ComputeHistogramOn(): (forces level >=2)
//  median, skew, energy, entropy
// ComputeTexturesOn():	 (forces level >=1)
//  texture(s)
//
// LEVEL 3 allows for boundary sharing features to be calculated, but this
// information is not stored in the feature structure and should be retrieved
// manually
//
// This class defaults to LEVEL 2 with all options set to OFF.
//*****************************************************************************

namespace ftk 
{

//Constructor
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::LabelImageToFeatures()
{
	intensityImage = NULL;
	labelImage = NULL;
	zernikeOrder = 5;
	
	boundaryPix.clear();
	interiorPix.clear();
	sharePix.clear();
	labels.clear();
	this->featureVals.clear();
	
	//Defaults:
	this->computationLevel = 2;
	this->computeHistogram = false;		
	this->computeTextures = false;
	this->surfacecomputation = false;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, bool CytoImage )
{
	cyto_image = CytoImage;

	typename LabelImageType::RegionType lblRegion = lblImgIn->GetRequestedRegion();
	typename IntensityImageType::RegionType intRegion = intImgIn->GetRequestedRegion();

	//Need to check size
	if( lblRegion != intRegion )
		return false;
		
	//Need to check regions:
	if( intRegion != intImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requested Region
		typedef itk::ExtractImageFilter< IntensityImageType, IntensityImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(intImgIn);
		cropFilter->SetExtractionRegion(intRegion);
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter1: " << err << std::endl;
		}

		intensityImage = cropFilter->GetOutput();
	}
	else
	{
		intensityImage = intImgIn;		//Use Full Image
	} 
	
	
	if( lblRegion != lblImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requestedd Region
		typedef itk::ExtractImageFilter< LabelImageType, LabelImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(lblImgIn);
		cropFilter->SetExtractionRegion(lblRegion);
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter2: " << err << std::endl;
		}

		labelImage = cropFilter->GetOutput();
	}
	else
	{
		labelImage = lblImgIn;
	}
	
	return true;
}




template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetCompleteImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, bool CytoImage )
{
	cyto_image = CytoImage;
		
	typename LabelImageType::RegionType lblRegion = lblImgIn->GetLargestPossibleRegion();
	typename IntensityImageType::RegionType intRegion = intImgIn->GetLargestPossibleRegion();

	//Need to check size
	if( lblRegion != intRegion )
		return false;
		
	//Need to check regions:
	if( intRegion != intImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requested Region
		typedef itk::ExtractImageFilter< IntensityImageType, IntensityImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(intImgIn);
		cropFilter->SetExtractionRegion(intRegion);
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter3: " << err << std::endl;
		}

		intensityImage = cropFilter->GetOutput();
	}
	else
	{
		intensityImage = intImgIn;		//Use Full Image
	} 
	
	
	if( lblRegion != lblImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requestedd Region
		typedef itk::ExtractImageFilter< LabelImageType, LabelImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(lblImgIn);
		cropFilter->SetExtractionRegion(lblRegion);
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter4: " << err << std::endl;
		}
		labelImage = cropFilter->GetOutput();
	}
	else
	{
		labelImage = lblImgIn;
	}
	
	return true;
}







template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, TLPixel index[VImageDimension], TLPixel size[VImageDimension], bool CytoImage)
{
	cyto_image = CytoImage;
	
	typename IntensityImageType::RegionType intRegion;
	typename IntensityImageType::IndexType intIndex;
	typename IntensityImageType::SizeType intSize;
	
	typename LabelImageType::RegionType labRegion;
	typename LabelImageType::IndexType labIndex;
	typename LabelImageType::SizeType labSize;

	for(unsigned int i=0; i<VImageDimension; ++i)
	{
		intIndex[i] = index[i];
		intSize[i] = size[i];
		labIndex[i] = index[i];
		labSize[i] = size[i];
	}
	//For the case when 2D images are passed as 3D
	if( intImgIn->GetLargestPossibleRegion().GetSize()[2] == 1 )
	{
		intSize[2] = 1;
		labIndex[2] = 0;
		intIndex[2] = 0;
		labSize[2] = 1;
	}

	intRegion.SetIndex(intIndex);
	intRegion.SetSize(intSize);
	labRegion.SetIndex(labIndex);
	labRegion.SetSize(labSize);
		
	//Need to check regions:
	if( intRegion != intImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requested Region
		typedef itk::ExtractImageFilter< IntensityImageType, IntensityImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(intImgIn);
		cropFilter->SetExtractionRegion(intRegion);
		cropFilter->SetDirectionCollapseToIdentity();
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter5: " << err << std::endl;
		}
		intensityImage = cropFilter->GetOutput();
	}
	else
	{
		intensityImage = intImgIn;		//Use Full Image
	} 
	
	if( labRegion != lblImgIn->GetBufferedRegion() )
	{
		//Crop Image to Requestedd Region
		typedef itk::ExtractImageFilter< LabelImageType, LabelImageType > CropFilterType;
		typename CropFilterType::Pointer cropFilter = CropFilterType::New();
		cropFilter->SetInput(lblImgIn);
		cropFilter->SetExtractionRegion(labRegion);
		cropFilter->SetDirectionCollapseToIdentity();
		try
		{
			cropFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Error in cropFilter6: " << err << std::endl;
		}
		labelImage = cropFilter->GetOutput();
	}
	else
	{
		labelImage = lblImgIn;
	}

	return true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ComputeHistogramOn()
{
	if(computationLevel <= 1)
		this->SetLevel(2);
	this->computeHistogram = true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ComputeSurfaceOn()
{
	this->surfacecomputation = true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ComputeTexturesOn()
{
	this->computeTextures = true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetLevel(short int newLevel)
{
	if(newLevel >= 3)
	{
		computationLevel = 3;
	}
	else if (newLevel <= 1)
	{
		computationLevel = 1;
	}
	else
	{
		computationLevel = 2;
	}
}


template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::Update()
{	

	if(!intensityImage || !labelImage)
	{
		std::cerr << " Please Set Input Images " << std::endl;
		return;
	}

	//LEVEL 1:
	if(computationLevel >= 1)
	{
		if(!RunLabelGeometryFilter()) return;	//Should throw exception
	}
	
	//LEVEL 2:
	if(computationLevel >= 2)
	{
		if(!RunLabelStatisticsFilter()) return; //Should throw exception
	}
	
	//LEVEL 3:
	if(computationLevel >= 3)
	{
		LabelImageScan();
		CalculateScanFeatures();
		RunZernikeFilter();
	}
	
	//TEXTURE CALCULATOR:
	if (this->computeTextures)
    {
		//if(!this->RunTextureFilter()) return;		//Should throw exception
    }
    
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetAdjacency()
{
	RunLabelGeometryFilter();
	LabelImageScan();
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
TLPixel LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetMaxLabel()
{
	LabelPixelType max = 0;
	
	for (int i = 0; i < (int)labels.size(); ++i)
	{
		if ( labels.at(i) > max )
			max = labels.at(i);
	}
	return max;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
IntrinsicFeatures * LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::GetFeatures( TLPixel label )
{
	typename FeatureMapType::iterator it;
	it = this->featureVals.find( label );
	if ( it == this->featureVals.end() )
    {
		// label does not exist, return a NULL value
		return NULL;
    }
	else
    {
		return &featureVals[ label ];
    }
}


template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::RunLabelGeometryFilter()
{
	typedef itk::ExtractImageFilter< IntensityImageType, Intensity2DType > DataExtractIntensity2DType;
	typedef itk::ExtractImageFilter< LabelImageType, Label2DType > DataExtractLabel2DType;

	if(!intensityImage || !labelImage)
		return false;

	typename Intensity2DType::Pointer intensity2D;
	typename Label2DType::Pointer label2D;

	//Check If Image is 2D
	if( labelImage->GetLargestPossibleRegion().GetSize()[2] == 1 )
	{
		typename DataExtractIntensity2DType::Pointer deFilter = DataExtractIntensity2DType::New();
		typename IntensityImageType::RegionType dRegion = intensityImage->GetLargestPossibleRegion();
		dRegion.SetSize(2,0);
		deFilter->SetExtractionRegion( dRegion );
		deFilter->SetDirectionCollapseToIdentity();
		deFilter->SetInput( intensityImage );

		typename DataExtractLabel2DType::Pointer deFilter1 = DataExtractLabel2DType::New();
		typename LabelImageType::RegionType dRegion1 = labelImage->GetLargestPossibleRegion();
		dRegion1.SetSize(2,0);
		deFilter1->SetExtractionRegion( dRegion1 );
		deFilter1->SetDirectionCollapseToIdentity();
		deFilter1->SetInput( labelImage );

		try
		{
			deFilter1->Update();
			deFilter->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << excp << std::endl;
			return false;
		}
		intensity2D = deFilter->GetOutput();
		label2D = deFilter1->GetOutput();

		typedef itk::LabelGeometryImageFilter< Label2DType, Intensity2DType > LabelGeometryType;
		typedef typename LabelGeometryType::Pointer LabelGeometryPointer;
		LabelGeometryPointer labelGeometryFilter = LabelGeometryType::New();		//Dirk's Filter;

		labelGeometryFilter->SetInput( label2D );
		labelGeometryFilter->SetIntensityInput( intensity2D );
		
		//SET ADVANCED (OPTIONAL) ITEMS FOR THIS FILTER:
		labelGeometryFilter->CalculatePixelIndicesOff();
		labelGeometryFilter->CalculateOrientedBoundingBoxOff();
		labelGeometryFilter->CalculateOrientedLabelRegionsOff();
		labelGeometryFilter->CalculateOrientedIntensityRegionsOff();

		//UPDATE THE FILTER	
		try
		{
			labelGeometryFilter->Update();
		}
		catch (itk::ExceptionObject & e) 
		{
			std::cerr << "Exception in Dirk's Label Geometry Filter: " << e << std::endl;
			return false;
		}

		//Now populate vector of labels & fill in the the map
		labels.clear();
		std::vector< typename LabelGeometryType::LabelPixelType > ls = labelGeometryFilter->GetLabels();
		for (int i = 0; i < (int)ls.size(); ++i)
		{
			TLPixel l = ls.at(i);
			labels.push_back( (unsigned short)l );
			LtoIMap[ l ] = 	i;
		}

		//Now extract the features:
		for (int i = 0; i < (int)labels.size(); ++i)
		{
			TLPixel label = labels.at(i);
			
			featureVals[label].ScalarFeatures[IntrinsicFeatures::VOLUME] = (float)labelGeometryFilter->GetVolume( label );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::INTEGRATED_INTENSITY] = (float)labelGeometryFilter->GetIntegratedIntensity( label );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ECCENTRICITY] = float( labelGeometryFilter->GetEccentricity( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ELONGATION] = float( labelGeometryFilter->GetElongation( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ORIENTATION] = float( labelGeometryFilter->GetOrientation( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::BBOX_VOLUME] = float( labelGeometryFilter->GetBoundingBoxVolume( label ) );
			
			typename LabelGeometryType::LabelPointType c = labelGeometryFilter->GetCentroid( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].Centroid[i] = float( c[i] );
			}

			c = labelGeometryFilter->GetWeightedCentroid( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].WeightedCentroid[i] = float( c[i] );
			}

			typename LabelGeometryType::AxesLengthType aL = labelGeometryFilter->GetAxesLength( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].AxisLength[i] = float( aL[i] );
			}

			typename LabelGeometryType::BoundingBoxType bbox = labelGeometryFilter->GetBoundingBox( label );
			for (unsigned int i = 0; i < VImageDimension*2; ++i)
			{
				featureVals[label].BoundingBox[i] = float( bbox[i] );
			}

		}
	}
	else
	{
		typedef itk::LabelGeometryImageFilter< LabelImageType, IntensityImageType > LabelGeometryType;
		typedef typename LabelGeometryType::Pointer LabelGeometryPointer;
		LabelGeometryPointer labelGeometryFilter = LabelGeometryType::New();		//Dirk's Filter;

		labelGeometryFilter->SetInput( labelImage );
		labelGeometryFilter->SetIntensityInput( intensityImage );
		
		//SET ADVANCED (OPTIONAL) ITEMS FOR THIS FILTER:
		labelGeometryFilter->CalculatePixelIndicesOff();
		labelGeometryFilter->CalculateOrientedBoundingBoxOff();
		labelGeometryFilter->CalculateOrientedLabelRegionsOff();
		labelGeometryFilter->CalculateOrientedIntensityRegionsOff();

		//UPDATE THE FILTER	
		try
		{
			labelGeometryFilter->Update();
		}
		catch (itk::ExceptionObject & e) 
		{
			std::cerr << "Exception in Dirk's Label Geometry Filter: " << e << std::endl;
			return false;
		}

		//Now populate vector of labels & fill in the the map
		labels.clear();
		std::vector< typename LabelGeometryType::LabelPixelType > ls = labelGeometryFilter->GetLabels();
		for (int i = 0; i < (int)ls.size(); ++i)
		{
			TLPixel l = ls.at(i);
			labels.push_back( (unsigned short)l );
			LtoIMap[ l ] = 	i;
		}

		//Now extract the features:
		//Surface feature Extraction
		if(this->surfacecomputation)
		{
			//this->RunSurfaceFeature(labelGeometryFilter);
		}
		//Geometry feature Extraction
		for (int i = 0; i < (int)labels.size(); ++i)
		{
			TLPixel label = labels.at(i);
			
			featureVals[label].ScalarFeatures[IntrinsicFeatures::VOLUME] = (float)labelGeometryFilter->GetVolume( label );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::INTEGRATED_INTENSITY] = (float)labelGeometryFilter->GetIntegratedIntensity( label );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ECCENTRICITY] = float( labelGeometryFilter->GetEccentricity( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ELONGATION] = float( labelGeometryFilter->GetElongation( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ORIENTATION] = float( labelGeometryFilter->GetOrientation( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::BBOX_VOLUME] = float( labelGeometryFilter->GetBoundingBoxVolume( label ) ); 
			
			typename LabelGeometryType::LabelPointType c = labelGeometryFilter->GetCentroid( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].Centroid[i] = float( c[i] );
			}

			c = labelGeometryFilter->GetWeightedCentroid( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].WeightedCentroid[i] = float( c[i] );
			}

			typename LabelGeometryType::AxesLengthType aL = labelGeometryFilter->GetAxesLength( label );
			for (unsigned int i = 0; i < VImageDimension; ++i)
			{
				featureVals[label].AxisLength[i] = float( aL[i] );
			}

			typename LabelGeometryType::BoundingBoxType bbox = labelGeometryFilter->GetBoundingBox( label );
			for (unsigned int i = 0; i < VImageDimension*2; ++i)
			{
				featureVals[label].BoundingBox[i] = float( bbox[i] );
			}

		}

	}

	return true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::RunLabelStatisticsFilter()
{
	if(!intensityImage || !labelImage) return false;
	
	// Run Label Statistics Image Filter (part of ITK)
	typedef itk::LabelStatisticsImageFilter< IntensityImageType , LabelImageType > LabelStatisticsType;
	typedef typename LabelStatisticsType::Pointer LabelStatisticsPointer;
	
	LabelStatisticsPointer labelStatisticsFilter = LabelStatisticsType::New();	//ITK label statistics filter
	labelStatisticsFilter->SetInput( intensityImage );
	labelStatisticsFilter->SetLabelInput( labelImage );
	
	//SET ADVANCED (OPTIONAL) ITEMS FOR THIS FILTER:
	if(computeHistogram)
	{
		labelStatisticsFilter->UseHistogramsOn();
		int numBins, lowerBound, upperBound;
		this->SetHistogramParameters(&numBins,&lowerBound, &upperBound);
		labelStatisticsFilter->SetHistogramParameters(numBins, lowerBound, upperBound);
	}
	else
	{
		labelStatisticsFilter->UseHistogramsOff();
	}
	
	//UPDATE THE FILTER	
	try
	{
		labelStatisticsFilter->Update();
	}
	catch (itk::ExceptionObject & e) 
	{
		std::cerr << "Exception in ITK Label Statistics Filter: " << e << std::endl;
		return false;
	}
	
	//Now populate the features information
	for (int i = 0; i < (int)labels.size(); ++i)
	{
		TLPixel label = labels.at(i);
		if (label <= 0) continue;
		
		featureVals[label].ScalarFeatures[IntrinsicFeatures::SUM] = (float)labelStatisticsFilter->GetSum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MEAN] = (float)labelStatisticsFilter->GetMean( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MEDIAN] = (float)labelStatisticsFilter->GetMedian( label ); //Should be 0 when no Histogram
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MINIMUM] = (float)labelStatisticsFilter->GetMinimum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MAXIMUM] = (float)labelStatisticsFilter->GetMaximum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::SIGMA] = (float)labelStatisticsFilter->GetSigma( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::VARIANCE] = (float)labelStatisticsFilter->GetVariance( label );
	
		if(computeHistogram)
		{
			typename LabelStatisticsType::HistogramPointer histo;
			double vol, mean, sigma, t_skew, t_energy, t_entropy;
			double log, diff, cube, prob;

			histo = labelStatisticsFilter->GetHistogram( label );
			vol = this->featureVals[label].ScalarFeatures[IntrinsicFeatures::VOLUME];
			mean = this->featureVals[label].ScalarFeatures[IntrinsicFeatures::MEAN];
			sigma = this->featureVals[label].ScalarFeatures[IntrinsicFeatures::SIGMA];

			t_skew = 0; t_energy = 0; t_entropy = 0;
			
			for (int i=0; i<(int)histo->Size(); ++i)
			{
				double v = histo->GetMeasurement(i,0);					//if not 1 bin for each value
				diff = v-mean;											//for skew
				cube = diff*diff*diff;									//for skew
				prob = histo->GetFrequency(i) / vol;					//for skew,energy,entropy
				t_skew = t_skew + (cube * prob);						//for skew
				t_energy = t_energy + (prob*prob);						//for energy
				if (prob == 0) log = 0;									//for entropy
				else log = log10(prob) / log10(double(2));				//for entropy
				t_entropy = t_entropy + (prob*log);						//for entropy
			}	
			if(sigma == 0)
				featureVals[label].ScalarFeatures[IntrinsicFeatures::SKEW] =  float(0);
			else
				featureVals[label].ScalarFeatures[IntrinsicFeatures::SKEW] =  float( t_skew / ( sigma * sigma * sigma ) );
				
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ENERGY] = float( t_energy );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::ENTROPY] = float( -1 * t_entropy );
		}
	}
	
	return true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetHistogramParameters(int* numBins, int* lowerBound, int* upperBound)
{
	//This function should use TIPixel (the type of pixels in the intensity image)
	// to set the histogram bins
	// FOR NOW ASSUME UNSIGNED VALUES OF UCHAR, USHORT
	*lowerBound = 0;
	
	if( typeid(TIPixel) == typeid(unsigned char) )
	{
		*upperBound = UCHAR_MAX;
	}
	else if ( typeid(TIPixel) == typeid(unsigned short) )
	{
		*upperBound = USHRT_MAX;
	}
	else
	{
		*upperBound = 1;
	}
	//*numBins = *upperBound - *lowerBound + 1;
	*numBins = 256;
}

//****************************************************************************************************
// void LabelImageScan(void)
// 
// This function scans through the label image and extracts the following information
// for each label found:
//  1. interior pixels
//  2. boundary pixels
//  3. Neighbor IDs (number of edges shared between label pairs)
//****************************************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::LabelImageScan()
{
	if(!labelImage) return;
	
	typedef itk::ConstantBoundaryCondition<LabelImageType> boundaryConditionType;
	typedef itk::ConstNeighborhoodIterator<LabelImageType, boundaryConditionType > NeighborhoodIteratorType;

	// The offsets for the neighboring pixels for 4-connectivity
	unsigned int dim;
	if ( cyto_image == true)
		dim = 4*VImageDimension;
	else
		dim = 2*VImageDimension;
		
	std::vector< typename NeighborhoodIteratorType::OffsetType > offsets(dim);
	for (unsigned int i=0; i<(dim); ++i)
	{
		offsets[i].Fill(0);
	}
	unsigned int p = 0, o = 0;
	while ( p < VImageDimension )
	{
		if ( cyto_image == true)
			offsets[o++][p] = -2;
		offsets[o++][p] = -1;
		offsets[o++][p] = 1;
		if ( cyto_image == true)
			offsets[o++][p] = 2;
		p++;
	}
	
	typename NeighborhoodIteratorType::RadiusType radius;
	if ( cyto_image == true)
		radius.Fill(2);
	else
		radius.Fill(1);	

	boundaryPix.clear();
	interiorPix.clear();
	int numLabels = (int)labels.size();
	boundaryPix.resize( numLabels );
	interiorPix.resize( numLabels );
	
	sharePix.clear();
	sharePix.resize( numLabels );

	NeighborhoodIteratorType it( radius, labelImage, labelImage->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) 
	{
		TLPixel v = it.GetCenterPixel();  // in the mask

		if ( v <= 0 ) continue;

		bool allSame = true;
		for (unsigned int i=0; i<dim; ++i)
		{
			TLPixel p = it.GetPixel( offsets[i] );
			
			if ( v != p )
			{
				allSame = false;
				int index = LtoIMap[v];
				typename std::map<TLPixel,int>::iterator loc = sharePix.at( index ).find(p);
				if( loc == sharePix.at( index ).end() )
					sharePix.at( index )[p] = 1;
				else
					++( (*loc).second );

			}
		}
	
		typename LabelImageType::IndexType index = it.GetIndex();
			
		if(allSame == false)
		{
			//I have a boundary point!!!
			boundaryPix[ LtoIMap[v] ].push_back(index);
		}
		else 
		{
			//interior point
			interiorPix[ LtoIMap[v] ].push_back(index);
		}
	}

}

//****************************************************************************************
// Calculate secondary features that depend on previous calculation:
//
// Features that depend on the Gradient mag (computed here):
// (also requires filters and labelScan to be complete)
// 1. Surface Gradient - average of surface gradients
// 2. interior gradient - average of interior gradients
//
// Features that use surface pixels list
// 1. Surface Area - no computation needed - just get number of surface voxels
// 2. Shape - (surface area)^3 / ( 36*pi*volume^2);
// 3. Eccentricity - ratio of max to min distance of the centroid to the object surface
// 4. Radius Variation - standard deviation of distance from surface voxels to centroid
//
// Features that use all object pixels:
// 1. equivalent distance - the average distance from all voxels to the centroid
// 2. distance variation - standard deviation of distance from all voxels to centroid
//****************************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::CalculateScanFeatures()
{
	if( !boundaryPix.size() || !labels.size() || !intensityImage)
		return;
		
	//FIRST CREATE the Gradient Magnitude Image:
	typedef float FloatPixelType;
	typedef itk::Image< FloatPixelType, VImageDimension > FloatImageType;
	typedef typename FloatImageType::Pointer FloatImagePointer;
	typedef itk::GradientMagnitudeImageFilter< IntensityImageType, FloatImageType > GradMagFilterType;
	typename GradMagFilterType::Pointer gmFilter = GradMagFilterType::New();
	gmFilter->SetInput( intensityImage );
	try
	{
		gmFilter->Update();
	}
	catch (itk::ExceptionObject & e) 
	{
		std::cerr << "Exception in Gradient Magnitude Filter: " << e << std::endl;
		return;
	}
	FloatImagePointer gmImage = gmFilter->GetOutput();				//Calculated gradient magnitude image

	//Local Variables:
	double sum_interior_grad, sum_surface_grad;				//gradient sums
	double sum_interior_intensity, sum_surface_intensity;	//intensity sums
	double max_bound_dist, min_bound_dist;					//from centroid to boundary 
	std::vector< double > bounddistances(0);				//from centroid to boundary
	std::vector< double > interiordistances(0);				//from centroid to interior

	float * centroid;
	typename LabelImageType::IndexType point;
	TLPixel currentLabel;
	
	for (int lab=0; lab<(int)labels.size(); ++lab)
	{
		currentLabel = labels.at(lab);
		if ((int)currentLabel <= 0) continue;

		centroid = this->featureVals[currentLabel].Centroid;

		max_bound_dist = 0.0;
		min_bound_dist = 100.0;
		bounddistances.clear();
		interiordistances.clear();
		sum_surface_grad = 0;
		sum_interior_grad = 0;
		sum_surface_intensity = 0;
		sum_interior_intensity = 0;

		//Get boundary pixel information for this label
		for (int i=0; i<(int)boundaryPix.at(LtoIMap[currentLabel]).size(); ++i)
		{
			point = boundaryPix.at(LtoIMap[currentLabel])[i];
			sum_surface_grad += gmImage->GetPixel(point);				//Add up surface gradients
			sum_surface_intensity += intensityImage->GetPixel(point);	//Add up intensity values

			//Calculate the distance to the centroid (city block)
			double dist = 0;
			for(unsigned int n = 0; n < VImageDimension; n++)
			{
				dist += fabs( (double)(point[n]) - (double)(centroid[n]) );
			}
			dist = sqrt(dist);
			bounddistances.push_back(dist);

			if      (dist < max_bound_dist)  min_bound_dist = dist;
			else if (dist > max_bound_dist)  max_bound_dist = dist;
		}



	//Added by Raghav : Convex Hull calculation.
		/*int ctr =0;
		int flag = 0;
		
		for (int i=0; i<(int)boundaryPix.at(LtoIMap[currentLabel]).size(); ++i)
		{
			++ctr;
			flag =1;
		}
		
		int myDimension = VImageDimension;
		
		// Need to correct for dimension for QHull to work correctly.
		if(labelImage->GetLargestPossibleRegion().GetSize()[2]==1)
		{
			myDimension = 2;
		}
		
		//int dim= myDimension;             /* dimension of points */
		//int numpoints;            /* number of points */
		//boolT ismalloc= True;    /* True if qhull should free points in qh_freeqhull() or reallocation */
		//char flags[250] = "qhull file.txt";          /* option flags for qhull, see qh_opt.htm */
		//FILE *outfile= NULL;    /* output from qh_produce_output() use NULL to skip qh_produce_output() */
		//FILE *errfile= fopen("qhull file.txt","w+");    /* error messages from qhull code */
		//int exitcode;             /* 0 if no error from qhull */
		//facetT *facet;            /* set by FORALLfacets */
		//int curlong, totlong;     /* memory remaining after qh_memfreeshort */
		//int i;
		/*	
		if(flag==1)
		{
		coordT *points = new coordT[(dim+1)*ctr];
		coordT *pointz;	
		ctr = 0;
		
		for (int i=0; i<(int)boundaryPix.at(LtoIMap[currentLabel]).size(); ++i)
		{		
				point = boundaryPix.at(LtoIMap[currentLabel])[i];
				pointz = points + ctr*myDimension;
				for(int i=myDimension; i-- ; )
				{	
					pointz[i] = static_cast<double>(point[i]);
				}
				ctr++;
			
		} 
		
		int cvHull = qh_new_qhull(myDimension, ctr, points, ismalloc,flags, outfile, errfile);
		qh_memfreeshort (&curlong, &totlong);
		if(cvHull!=0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CONVEXITY] = (featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::VOLUME])/cvHull;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CONVEXITY] = 0.0; // Changed By Amin, it gives inf as an output.
		
		delete [] points;
		}
		
	fclose(errfile);*/


		//Get interior pixel information for this label
		for (int i=0; i<(int)interiorPix.at(LtoIMap[currentLabel]).size(); ++i)
		{
			point = interiorPix.at(LtoIMap[currentLabel])[i];
			sum_interior_grad += gmImage->GetPixel(point);				//Add up interior gradients
			sum_interior_intensity += intensityImage->GetPixel(point);	//Add up intensity values
			
			//Calculate the distance to the centroid (city block)
			double dist = 0;
			for(unsigned int n = 0; n < VImageDimension; n++)
			{
				dist += fabs( (double)(point[n]) - (double)(centroid[n]) );
			}
			dist = sqrt(dist);
			interiordistances.push_back(dist);
		}

		//Compute Gradient information for this label
		if( boundaryPix[ LtoIMap[currentLabel] ].size() == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_GRADIENT] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_GRADIENT] = float( sum_surface_grad / boundaryPix[ LtoIMap[currentLabel] ].size() );
			
		if( interiorPix[ LtoIMap[currentLabel] ].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_GRADIENT] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_GRADIENT] = float( sum_interior_grad / interiorPix[ LtoIMap[currentLabel] ].size() );

		//Compute Intensities
		if( boundaryPix[ LtoIMap[currentLabel] ].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY]	= float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY] = float( sum_surface_intensity / boundaryPix[ LtoIMap[currentLabel] ].size() );
		
		if( interiorPix[ LtoIMap[currentLabel] ].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] = float( sum_interior_intensity / interiorPix[ LtoIMap[currentLabel] ].size() );
		
		if (featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTENSITY_RATIO] = 0;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTENSITY_RATIO] \
				= this->featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY] \
				/ this->featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY];

		//Find mean and std-dev of distances
		double interior_sum = 0;
		double surface_sum = 0;
		for (int i=0; i<(int)bounddistances.size(); ++i)
		{
			surface_sum += bounddistances[i];
		}
		for (int i=0; i<(int)interiordistances.size(); ++i)
		{
			interior_sum += interiordistances[i];
		}
		double surface_mean;
		if(bounddistances.size() == 0)
			surface_mean = 0;
		else
			surface_mean = surface_sum / bounddistances.size();
		
		double interior_mean;
		if(interiordistances.size())
			interior_mean = 0;
		else
			interior_mean = interior_sum / interiordistances.size();
		
		//allFeatures[currentLabel].averagedistance = (surface_mean + interior_mean) / 2;
		
		double surface_sq_sum = 0;
		for (int i=0; i<(int)bounddistances.size(); ++i)
		{
			double diff = bounddistances[i] - surface_mean;
			double sq = diff*diff;
			surface_sq_sum += sq;
		}
		if(bounddistances.size() == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::RADIUS_VARIATION] = 0;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::RADIUS_VARIATION] = float( sqrt( surface_sq_sum / bounddistances.size() ) );
		
		//double interior_sq_sum = 0;
		//for (int i=0; i<(int)interiordistances.size(); ++i)
		//{
		//	double diff = interiordistances[i] - interior_mean;
		//	double sq = diff*diff;
		//	interior_sq_sum += sq;
		//}
		//float interiorvariation = float( sqrt( interior_sq_sum / interiordistances.size() ) );
		//allFeatures[currentLabel].distancevariation = (allFeatures[currentLabel].radiusvariation + interiorvariation) / 2;
		
		//surface area:
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_AREA] = float( boundaryPix[ LtoIMap[currentLabel] ].size() );
		
		//shape:
		double sa = this->featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_AREA];
		double pi = 3.1415;
		double v = this->featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::VOLUME];
		if(v == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHAPE] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHAPE] = float( sa*sa*sa / ( 36*pi*v*v) );
		
		//percent shared boundary:
		int index = LtoIMap[currentLabel];
		int zeroBound = sharePix.at( index )[0];
		int nonzeroBound = 0;
		typename std::map<TLPixel, int>::iterator it;
		for (it=sharePix.at( index ).begin(); it!=sharePix.at( index ).end(); ++it)
		{
			if( (*it).first != 0 )
				nonzeroBound += (*it).second;
		}
		if(zeroBound == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHARED_BOUNDARY] = 0;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHARED_BOUNDARY] = float(nonzeroBound) / float(nonzeroBound+zeroBound);
	}
	
	//CLEAR THESE BECAUSE I DON'T USE THEM AGAIN
	boundaryPix.clear();
	interiorPix.clear();
}

//**************************************************************************
// RUN THE ITK TEXTURE FILTER
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::RunTextureFilter()
{
	if(!intensityImage || !labelImage) return false;
	if(!labels.size()) return false;
	if(	labelImage->GetLargestPossibleRegion().GetSize()[0]==1 ||
		labelImage->GetLargestPossibleRegion().GetSize()[1]==1 ||
		labelImage->GetLargestPossibleRegion().GetSize()[2]==1	)
  {
    return this->RunTextureFilter2D(); 
  }
    
	LabelImagePointer tempIntensityImage;

	typedef itk::RescaleIntensityImageFilter< IntensityImageType, LabelImageType > RescaleFilterType;
	typedef typename RescaleFilterType::Pointer RescaleFilterPointer;
	RescaleFilterPointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->SetInput(intensityImage);
	rescaleFilter->Update();
	tempIntensityImage = rescaleFilter->GetOutput();
	
	typedef itk::Statistics::ScalarImageToTextureFeaturesFilter< LabelImageType > TextureCalcType;
	typedef typename TextureCalcType::Pointer TextureCalcPointer;
	TextureCalcPointer textureCalculator = TextureCalcType::New();
	textureCalculator->SetInput(tempIntensityImage);
	textureCalculator->SetMaskImage(labelImage);
	textureCalculator->SetPixelValueMinMax(0,255);
	textureCalculator->SetFastCalculations(true);
	
	/*typedef itk::ImageFileWriter<LabelImageType> LabelImageFileWriterType;
	typedef LabelImageFileWriterType::Pointer LabelImageFileWriter;
	LabelImageFileWriter label_image_file_writer = LabelImageFileWriterType::New();
	
	label_image_file_writer->SetFileName("labelImage.tif");
	label_image_file_writer->SetInput(labelImage);

	try
	{
		label_image_file_writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "label_image_file_writer error: " << err << std::endl;
		//throw std::exception;
	}

	label_image_file_writer->SetFileName("tempIntensityImage.tif");
	label_image_file_writer->SetInput(tempIntensityImage);

	try
	{
		label_image_file_writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "label_image_file_writer error: " << err << std::endl;
		//throw std::exception;
	}*/
	
	TLPixel currentLabel;
	for (int lab=0; lab<(int)labels.size(); ++lab)
	{
		currentLabel = labels.at(lab);
		if ((int)currentLabel <= 0) continue;
		
		typename LabelImageType::RegionType region;
		typename LabelImageType::SizeType size;
		typename LabelImageType::IndexType index;

		int loc = 0;
		for (unsigned int dim=0; dim<VImageDimension; dim++)
		{
			index[dim] = (long int)this->featureVals[currentLabel].BoundingBox[loc];	//bbox min
			size[dim] = (long unsigned int)this->featureVals[currentLabel].BoundingBox[loc+1]- index[dim] + 1;	//bbox max - min + 1
			loc += 2;
		}

        region.SetSize(size);
        region.SetIndex(index);
        labelImage->SetRequestedRegion(region);
        tempIntensityImage->SetRequestedRegion(region);
		textureCalculator->SetInsidePixelValue( currentLabel );
		
		try
		{
			textureCalculator->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cerr << "textureCalculator error: " << err << std::endl;
		}
		/*
		FeatureValueVector * 	GetFeatureMeans ()
		FeatureValueVector * 	GetFeatureStandardDeviations ()
		const FeatureNameVector * 	GetRequestedFeatures ()
		SetRequestedFeatures (const FeatureNameVector *_arg)
		*/
		
		typename TextureCalcType::FeatureValueVector *tex = textureCalculator->GetFeatureMeans();
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENERGY] = float( tex->ElementAt(0) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENTROPY] = float( tex->ElementAt(1) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INVERSE_DIFFERENCE_MOMENT] = float( tex->ElementAt(2) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INERTIA] = float( tex->ElementAt(3) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_SHADE] = float( tex->ElementAt(4) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_PROMINENCE] = float( tex->ElementAt(5) );
	}
	
	tempIntensityImage = 0;
	return true;
}

//**************************************************************************
// RUN THE ITK TEXTURE FILTER FOR A 2D IMAGE
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::RunTextureFilter2D()
{
	
  //define our 2D image types
  typedef itk::Image< IntensityPixelType, 2 > IntensityImageType2D;
	typedef typename IntensityImageType2D::Pointer IntensityImagePointer2D;
	typedef itk::Image< LabelPixelType, 2 > LabelImageType2D;
	typedef typename LabelImageType2D::Pointer LabelImagePointer2D;
	
  //recast labelImage into two dimensions
  typename LabelImageType::IndexType labelStart;
  labelStart.Fill(0);
  typename LabelImageType::SizeType labelSize; 
  labelSize = labelImage->GetLargestPossibleRegion().GetSize();
  labelSize[2] = 0;  //this tells the filter to collapse the Z-dimension
  
  typename LabelImageType::RegionType labelRegion;
  labelRegion.SetSize(labelSize);
  labelRegion.SetIndex(labelStart);

  typedef itk::ExtractImageFilter< LabelImageType, LabelImageType2D >
    ExtractLabelFilterType;
  typename ExtractLabelFilterType::Pointer extractLabel = ExtractLabelFilterType::New();
  extractLabel->SetExtractionRegion(labelRegion);
  extractLabel->SetInput(labelImage);
  extractLabel->SetDirectionCollapseToIdentity(); // This is required.
  extractLabel->Update();
  LabelImagePointer2D labelImage2D = extractLabel->GetOutput();

  //do the same for intensityImage
  typename IntensityImageType::IndexType intensityStart;
  intensityStart.Fill(0);
  typename IntensityImageType::SizeType intensitySize;
  intensitySize = intensityImage->GetLargestPossibleRegion().GetSize();
  intensitySize[2] = 0;  //this tells the filter to collapse the Z-dimension
  
  typename IntensityImageType::RegionType intensityRegion;
  intensityRegion.SetSize(intensitySize);
  intensityRegion.SetIndex(intensityStart);

  typedef itk::ExtractImageFilter< IntensityImageType, IntensityImageType2D >
    ExtractIntensityFilterType;
  typename ExtractIntensityFilterType::Pointer extractIntensity = ExtractIntensityFilterType::New();
  extractIntensity->SetExtractionRegion(intensityRegion);
  extractIntensity->SetInput(intensityImage);
  extractIntensity->SetDirectionCollapseToIdentity(); // This is required.
  extractIntensity->Update();
  IntensityImagePointer2D intensityImage2D = extractIntensity->GetOutput();

	typedef itk::RescaleIntensityImageFilter< IntensityImageType2D, LabelImageType2D > RescaleFilterType;
	typedef typename RescaleFilterType::Pointer RescaleFilterPointer;
	RescaleFilterPointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->SetInput(intensityImage2D);
	rescaleFilter->Update();
	LabelImagePointer2D tempIntensityImage2D = rescaleFilter->GetOutput();
	
	typedef itk::Statistics::ScalarImageToTextureFeaturesFilter< LabelImageType2D > TextureCalcType;
	typedef typename TextureCalcType::Pointer TextureCalcPointer;
	TextureCalcPointer textureCalculator = TextureCalcType::New();
	textureCalculator->SetInput(tempIntensityImage2D);
	textureCalculator->SetMaskImage(labelImage2D);
	textureCalculator->SetPixelValueMinMax(0,255);
	textureCalculator->SetFastCalculations(true);
	
	TLPixel currentLabel;
	for (int lab=0; lab<(int)labels.size(); ++lab)
	{
		currentLabel = labels.at(lab);
		if ((int)currentLabel <= 0) continue;
		
		typename LabelImageType2D::RegionType region;
		typename LabelImageType2D::SizeType size;
		typename LabelImageType2D::IndexType index;

		int loc = 0;
		for (unsigned int dim=0; dim<2; dim++)
		{
			index[dim] = (long int)this->featureVals[currentLabel].BoundingBox[loc];	//bbox min
			size[dim] = (long unsigned int)this->featureVals[currentLabel].BoundingBox[loc+1]- index[dim] + 1;	//bbox max - min + 1
			loc += 2;
		}

        region.SetSize(size);
        region.SetIndex(index);
        labelImage2D->SetRequestedRegion(region);
        tempIntensityImage2D->SetRequestedRegion(region);
		textureCalculator->SetInsidePixelValue( currentLabel );
		
		try
		{
			textureCalculator->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cerr << "textureCalculator error: " << err << std::endl;
		}
		
		typename TextureCalcType::FeatureValueVector *tex = textureCalculator->GetFeatureMeans();
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENERGY] = float( tex->ElementAt(0) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENTROPY] = float( tex->ElementAt(1) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INVERSE_DIFFERENCE_MOMENT] = float( tex->ElementAt(2) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INERTIA] = float( tex->ElementAt(3) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_SHADE] = float( tex->ElementAt(4) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_PROMINENCE] = float( tex->ElementAt(5) );
	}
	
	tempIntensityImage2D = 0;
	return true;
}


//**************************************************************************
// RUN THE ZERNIKE FILTER
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::RunZernikeFilter()
{
	#ifdef ZERNIKE
	{
		typedef zernike::ImageType zernikeImageType;
		typedef itk::ImageRegionIterator< zernikeImageType> zerIteratorType;
		
		if(!intensityImage || !labelImage) return;
		if(!labels.size()) return;
		
		IntensityImagePointer nucleusImage;
		for(int i = 0; i<(int)labels.size(); ++i)
		{
			TLPixel label = labels.at(i);
			int bb[6];
			for(int co = 0; co < 6; co++)
			{
				bb[co] = int(featureVals[label].BoundingBox[co]+0.5);
			}
			
			LabelImageType::IndexType labIndex;
			LabelImageType::SizeType labSize;
			LabelImageType::RegionType labRegion;
			labIndex[0] = bb[0];
			labIndex[1] = bb[2];
			labIndex[2] = bb[4];
			labSize[0] = bb[1]-bb[0]+1;
			labSize[1] = bb[3]-bb[2]+1;
			labSize[2] = bb[5]-bb[4]+1;
			labRegion.SetIndex(labIndex);
			labRegion.SetSize(labSize);
			labIteratorType labIter(labelImage,labRegion);
			
			IntensityImageType::IndexType intIndex;
			IntensityImageType::SizeType intSize;
			IntensityImageType::RegionType intRegion;
			intIndex[0] = bb[0];
			intIndex[1] = bb[2];
			intIndex[2] = bb[4];
			intSize[0] = bb[1]-bb[0]+1;
			intSize[1] = bb[3]-bb[2]+1;
			intSize[2] = bb[5]-bb[4]+1;
			intRegion.SetIndex(intIndex);
			intRegion.SetSize(intSize);
			intIteratorType intIter(intensityImage,intRegion);
			
			zernikeImageType::Pointer zerImg = zernikeImageType::New();
			zernikeImageType::IndexType zerIndex;
			zernikeImageType::SizeType zerSize;
			zernikeImageType::RegionType zerRegion;
			zerIndex.Fill(0);
			zerSize[0] = bb[1]-bb[0]+1;
			zerSize[1] = bb[3]-bb[2]+1;
			zerSize[2] = bb[5]-bb[4]+1;
			zerRegion.SetIndex(zerIndex);
			zerRegion.SetSize(zerSize);
			zerImg->SetRegions(zerRegion);
			zerImg->Allocate();
			zerImg->FillBuffer(0);
			zerIteratorType zerIter(zerImg,zerImg->GetLargestPossibleRegion());
			
			for(;!zerIter.IsAtEnd(); ++zerIter,++intIter,++labIter)
			{
				if(labIter.Get() == label)
				{
					zerIter.Set(intIter.Get());
				}
				else
				{
					zerIter.Set(0);
				}
			}
			
			zernike* myzernike = new zernike( zerImg, zernikeOrder);
			IDtoZernikeMap[label] = myzernike->GetZernike();		
	
		}
	}
	#endif
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::RunSurfaceFeature(LabelGeometryPointer labelGeometryFilter)
{
	// Iterator over the label image.
	typedef itk::ImageRegionConstIterator< LabelImageType> labelItType;	
	typedef itk::ImageRegionConstIteratorWithIndex< LabelImageType > IndexItType;

	labelItType labelIt( labelGeometryFilter->GetInput(), labelGeometryFilter->GetInput()->GetBufferedRegion() );
	labelItType labelItxp;  
	labelItType labelItxn;
	labelItType labelItyp;
	labelItType labelItyn;
	labelItType labelItzp;
	labelItType labelItzn;

	typename LabelGeometryType::LabelPixelType label;
	typename LabelGeometryType::LabelPixelType labelxp;
	typename LabelGeometryType::LabelPixelType labelxn;
	typename LabelGeometryType::LabelPixelType labelyp;
	typename LabelGeometryType::LabelPixelType labelyn;
	typename LabelGeometryType::LabelPixelType labelzp;
	typename LabelGeometryType::LabelPixelType labelzn;

	typename IndexItType::IndexType index;
	typename IndexItType::IndexType indexxp;
	typename IndexItType::IndexType indexxn;
	typename IndexItType::IndexType indexyp;
	typename IndexItType::IndexType indexyn;
	typename IndexItType::IndexType indexzp;
	typename IndexItType::IndexType indexzn;

	typename LabelImageType::SizeType labelSize; 
	labelSize = labelImage->GetLargestPossibleRegion().GetSize();
	itk::SizeValueType sizex = labelSize[0];
	itk::SizeValueType sizey = labelSize[1];
	itk::SizeValueType sizez = labelSize[2];

	int counterx = 0;
	int countery = 0;
	int counterz = 0;

	// Do the work
	std::vector< std::vector<double[3] > > surfacecoordinates;
	while ( !labelIt.IsAtEnd() )
    {
		label = labelIt.Get();
		index = labelIt.GetIndex();
		//if(label != 0)
		//{		
		//	bool surface = false;
		//	if(counterx != 0)
		//	{
		//		labelItxn = labelIt - 1;
		//		labelxn = labelItxn.Get();
		//		indexxn = labelItxn.GetIndex();
		//		if(indexxn[0] != index[0])
		//			surface = true;
		//	}
		//	if(counterx != sizex - 1)
		//	{
		//		labelItxp = labelIt + 1;
		//		labelxp = labelItxp.Get();
		//		indexxp = labelItxp.GetIndex();
		//		if(indexxp[0] != index[0])
		//			surface = true;
		//	}
		//	
		//	if(countery != 0)
		//	{
		//		labelItyn = labelIt - sizex;
		//		labelyn = labelItyn.Get();
		//		indexyn = labelItyn.GetIndex();
		//		if(indexyn[1] != index[1])
		//			surface = true;
		//	}
		//	if(countery != sizey - 1)
		//	{
		//		labelItyp = labelIt + sizex;
		//		labelyp = labelItyp.Get();
		//		indexyp = labelItyp.GetIndex();
		//		if(indexyp[1] != index[1])
		//			surface = true;
		//	}

		//	if(counterz != 0)
		//	{
		//		labelItzn = labelIt - sizex*sizey;
		//		labelzn = labelItzn.Get();
		//		indexzn = labelItzn.GetIndex();
		//		if(indexzn[2] != index[2])
		//			surface = true;
		//	}
		//	if(counterz != sizez - 1)
		//	{
		//		labelItzp = labelIt + sizex*sizey;
		//		labelzp = labelItzp.Get();
		//		indexzp = labelItzp.GetIndex();
		//		if(indexzp[2] != index[2])
		//			surface = true;
		//	}
		//	if(surface == true)
		//	{
		//		if(surfacecoordinates.size() < label)
		//			surfacecoordinates.resize (label)
		//		double coordinate[3];
		//		coordinate[0] = index[0];
		//		coordinate[1] = index[1];
		//		coordinate[2] = index[2];
		//		surfacecoordinates[label].push_back(coordinate);
		//	}
		//}

		counterx++;
		if(counterx == sizex - 1)
		{
			counterx = 0;
			countery++;
		}
		if(countery == sizey - 1)
		{
			countery = 0;
			counterz++;
		}
		labelIt++;
	}

	///////////////////////////////////////////////write to files;
	const char* filename = "Surface_Coordinate.txt";
	FILE *fp = fopen(filename,"w");
	for(int i = 0; i < 20; i++)
	{
		fprintf(fp,"%s", "lable ");
		fprintf(fp,"\t");
		fprintf(fp,"%d", label);
		fprintf(fp,"\n");
		for(int j = 0; j < surfacecoordinates[i].size (); j++)
		{
			fprintf(fp,"%d", surfacecoordinates[i][j][0]);
			fprintf(fp,"\t");
			fprintf(fp,"%d", surfacecoordinates[i][j][1]);
			fprintf(fp,"\t");
			fprintf(fp,"%d", surfacecoordinates[i][j][2]);
			fprintf(fp,"\n");
		}
	}	
	fclose(fp);
}

//**************************************************************************
// This function will compute the percent of the object's boundary with label 
// 'focusLabel' is shared with the 'neighborLabel' object
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
float LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetPercentSharedBoundary(TLPixel focusLabel, TLPixel neighborLabel)
{
	if(sharePix.size() <= 0)
		return 0;
	
	//Check to be sure focusLabel and neighborLabel exist:
	if(	LtoIMap.find(focusLabel) == LtoIMap.end() || \
		LtoIMap.find(neighborLabel) == LtoIMap.end() )
	{
		return 0;
	}
	 
	int totalBound = 0;
	int sharedBound = 0;
	int index = LtoIMap[focusLabel];
	typename std::map<TLPixel, int>::iterator it;
	for (it=sharePix.at(index).begin(); it!=sharePix.at(index).end(); ++it)
	{
		totalBound += (*it).second;
		if( (*it).first == neighborLabel )
			sharedBound = (*it).second;
	}

	return ( float(sharedBound) / float(totalBound) );
}

//**************************************************************************
// This function returns the object labels that 'label' touches
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
std::vector<TLPixel> LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetContactNeighbors(TLPixel label)
{
	std::vector<TLPixel> nbs(0);
	
	if(sharePix.size() <= 0)
	{
		return nbs;
	}

	//Check to be sure label exists:
	if(	LtoIMap.find(label) == LtoIMap.end() )
	{
		return nbs;
	}
	
	int index = LtoIMap[label];
	typename std::map<TLPixel, int>::iterator it;
	for (it=sharePix.at(index).begin(); it!=sharePix.at(index).end(); ++it)
	{
		nbs.push_back( (*it).first );
	}
	return nbs;
}

} //end namespace ftk
#endif
