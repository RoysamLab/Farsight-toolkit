#ifndef _ftkLabelImageToFeatures_txx
#define _ftkLabelImageToFeatures_txx
#include "ftkLabelImageToFeatures.h"

#include <itkConstNeighborhoodIterator.h>
#include <itkConstantBoundaryCondition.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkExtractImageFilter.h>

//#include <math.h>

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
// ComputeHistogramOn(): (forces level 2)
//  median, skew, energy, entropy
// ComputeAdvancedOn():	 (forces level 3)
//  solidity, texture(s)
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
	gmImage = NULL;
	
	labelGeometryFilter = NULL;
	labelStatisticsFilter = NULL;
	
	boundaryPix.clear();
	interiorPix.clear();
	sharePix.clear();
	labels.clear();
	featureVals.clear();
	
	//Defaults:
	computationLevel = 2;
	computeHistogram = false;					
	computeAdvanced = false;
						
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn )
{
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
		cropFilter->Update();
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
		cropFilter->Update();
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
	computeHistogram = true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ComputeAdvancedOn()
{
	this->SetLevel(3);
	computeAdvanced = true;
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
		RunLabelGeometryFilter();
		ReadLabelGeometryFeatures();
	}
	
	//LEVEL 2:
	if(computationLevel >= 2)
	{
		RunLabelStatisticsFilter();
		ReadLabelStatisticsFeatures();
	}
	
	//LEVEL 3:
	if(computationLevel >= 3)
	{
		CreateGradientMagnitudeImage();
		LabelImageScan();
		CalculateScanFeatures();
	}
	
	//OPTIONS:
	if(computeHistogram)
		CalculateHistogramFeatures();
	
	if(computeAdvanced)
	{
		RunTextureFilter();
	}
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
	it = featureVals.find( label );
	if ( it == featureVals.end() )
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
	if(!intensityImage || !labelImage) return false;
	
	//Run Dirk's Label Image Geometry Filter
	if(labelGeometryFilter) labelGeometryFilter->Delete();
	labelGeometryFilter = LabelGeometryType::New();
	labelGeometryFilter->SetInput( labelImage );
	labelGeometryFilter->SetIntensityInput( intensityImage );
	
	//SET ADVANCED (OPTIONAL) ITEMS FOR THIS FILTER:
	if(computeAdvanced)
	{
		labelGeometryFilter->CalculateOrientedBoundingBoxOn();
	}
	else
	{
		labelGeometryFilter->CalculatePixelIndicesOff();
		labelGeometryFilter->CalculateOrientedBoundingBoxOff();
	}
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
		labelGeometryFilter->Delete();
		labelGeometryFilter = NULL;
		return false;
	}
	
	//Now populate vector of labels
	this->labels.clear();
	std::vector< typename LabelGeometryType::LabelPixelType > ls = labelGeometryFilter->GetLabels();
	for (int i = 0; i < (int)ls.size(); ++i)
	{
		this->labels.push_back( (unsigned short)ls.at(i) );
	}
	
	return true;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::RunLabelStatisticsFilter()
{
	if(!intensityImage || !labelImage) return false;
	// Run Label Statistics Image Filter (part of ITK)
	if(labelStatisticsFilter) labelStatisticsFilter->Delete();
	labelStatisticsFilter  = LabelStatisticsType::New();
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
	*numBins = *upperBound - *lowerBound - 1;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::CreateGradientMagnitudeImage()
{
	if(!intensityImage) return false;		//If no intensity image do nothing
	
	if(gmImage) gmImage->Delete();			//If gmImage exists, delete it
	
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
		return false;
	}
	
	gmImage = gmFilter->GetOutput();
	return true;
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
	std::vector< typename NeighborhoodIteratorType::OffsetType > offsets(2*VImageDimension);
	for (int i=0; i<(2*VImageDimension); ++i)
	{
		offsets[i].Fill(0);
	}
	int p = 0, o = 0;
	while ( p < VImageDimension )
	{
		offsets[o++][p] = -1;
		offsets[o++][p] = 1;
		p++;
	}
	
	typename NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);

	boundaryPix.clear();
	interiorPix.clear();
	int maxLabel = (int)GetMaxLabel();
	boundaryPix.resize( maxLabel + 1 );
	interiorPix.resize( maxLabel + 1 );
	
	sharePix.clear();
	sharePix.resize( maxLabel + 1 );
	for (int i=0; i<(int)sharePix.size(); ++i)
		sharePix[i].resize(i);

	NeighborhoodIteratorType it( radius, labelImage, labelImage->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) 
	{
		TLPixel v = it.GetCenterPixel();  // in the mask

		if ( v <= 0 ) continue;

		bool allSame = true;
		for (int i=0; i<2*VImageDimension; ++i)
		{
			TLPixel p = it.GetPixel( offsets[i] );
			
			if ( v != p )
			{
				allSame = false;
				if ( v > p )
					++sharePix[v][p];
				else if ( v < p )
					++sharePix[p][v];
			}
		}
	
		typename LabelImageType::IndexType index = it.GetIndex();
			
		if(allSame == false)
		{
			//I have a boundary point!!!
			boundaryPix[v].push_back(index);
		}
		else 
		{
			//interior point
			interiorPix[v].push_back(index);
		}
	}

}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ReadLabelGeometryFeatures()
{
	if(!labelGeometryFilter) return;
	
	//Now populate the features information
	for (int i = 0; i < (int)labels.size(); ++i)
	{
		TLPixel label = labels.at(i);
		
		featureVals[label].ScalarFeatures[IntrinsicFeatures::VOLUME] = (float)labelGeometryFilter->GetVolume( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::INTEGRATED_INTENSITY] = (float)labelGeometryFilter->GetIntegratedIntensity( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::ECCENTRICITY] = float( labelGeometryFilter->GetEccentricity( label ) );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::ELONGATION] = float( labelGeometryFilter->GetElongation( label ) );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::ORIENTATION] = float( labelGeometryFilter->GetOrientation( label ) );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::BBOX_VOLUME] = float( labelGeometryFilter->GetBoundingBoxVolume( label ) );
		
		string n = IntrinsicFeatures::Info[IntrinsicFeatures::VOLUME].name;

		if(computeAdvanced)
		{
			double objVol = double( labelGeometryFilter->GetVolume( label ) );
			double boxVol = double( labelGeometryFilter->GetOrientedBoundingBoxVolume( label ) );
			featureVals[label].ScalarFeatures[IntrinsicFeatures::SOLIDITY] = float( objVol / boxVol );
		}
		
		this->GetCentroid( label );
		this->GetWeightedCentroid( label );
		this->GetAxisLength( label );
		this->GetBoundingBox( label );
	}
	
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::ReadLabelStatisticsFeatures()
{
	if(!labelStatisticsFilter) return;
	
	//Now populate the features information
	for (int i = 0; i < (int)labels.size(); ++i)
	{
		TLPixel label = labels.at(i);
		
		featureVals[label].ScalarFeatures[IntrinsicFeatures::SUM] = (float)labelStatisticsFilter->GetSum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MEAN] = (float)labelStatisticsFilter->GetMean( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MEDIAN] = (float)labelStatisticsFilter->GetMedian( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MINIMUM] = (float)labelStatisticsFilter->GetMinimum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::MAXIMUM] = (float)labelStatisticsFilter->GetMaximum( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::SIGMA] = (float)labelStatisticsFilter->GetSigma( label );
		featureVals[label].ScalarFeatures[IntrinsicFeatures::VARIANCE] = (float)labelStatisticsFilter->GetVariance( label );
	}
}

//****************************************************************************************
// Calculate secondary features that depend on previous calculation:
//
// Features that depend on the Gradiant mag:
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
	if( !boundaryPix.size() || !labels.size() )
		return;

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

		centroid = featureVals[currentLabel].Centroid;

		max_bound_dist = 0.0;
		min_bound_dist = 100.0;
		bounddistances.clear();
		interiordistances.clear();
		sum_surface_grad = 0;
		sum_interior_grad = 0;
		sum_surface_intensity = 0;
		sum_interior_intensity = 0;

		//Get boundary pixel information for this label
		for (int i=0; i<(int)boundaryPix[currentLabel].size(); ++i)
		{
			point = boundaryPix[currentLabel][i];
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

		//Get interior pixel information for this label
		for (int i=0; i<(int)interiorPix[currentLabel].size(); ++i)
		{
			point = interiorPix[currentLabel][i];
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
		if( boundaryPix[currentLabel].size() == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_GRADIENT] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_GRADIENT] = float( sum_surface_grad / boundaryPix[currentLabel].size() );
			
		if( interiorPix[currentLabel].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_GRADIENT] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_GRADIENT] = float( sum_interior_grad / interiorPix[currentLabel].size() );

		//Compute Intensities
		if( boundaryPix[currentLabel].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY]	= float( 0 );
		else
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY] = float( sum_surface_intensity / boundaryPix[currentLabel].size() );
		
		if( interiorPix[currentLabel].size() == 0 )
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] = float( sum_interior_intensity / interiorPix[currentLabel].size() );
		
		if (featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY] == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTENSITY_RATIO] = 0;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTENSITY_RATIO] \
				= featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_INTENSITY] \
				/ featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INTERIOR_INTENSITY];

		/*
		//Now use min/max to calculate the eccentricity
		if (min_dist != 0)
			eccentricityF14[currentLabel] = max_bound_dist / min_bound_dist;
		else
			eccentricityF14[currentLabel] = 100;
		*/

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
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_AREA] = float( boundaryPix[currentLabel].size() );
		
		//shape:
		double sa = featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SURFACE_AREA];
		double pi = 3.1415;
		double v = featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::VOLUME];
		if(v == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHAPE] = float( 0 );
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHAPE] = float( sa*sa*sa / ( 36*pi*v*v) );
		
		//percent shared boundary:
		int zeroBound = sharePix[currentLabel][0];
		int nonzeroBound = 0;
		for (TLPixel i=1; i<currentLabel; ++i)
		{
			nonzeroBound += sharePix[currentLabel][i];
		}
		for (TLPixel i=currentLabel+1; i<=GetMaxLabel(); ++i)
		{
			nonzeroBound = nonzeroBound + sharePix[i][currentLabel];
		}
		if(zeroBound == 0)
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHARED_BOUNDARY] = 0;
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SHARED_BOUNDARY] = float(nonzeroBound) / float(nonzeroBound+zeroBound);
	}
}

//**************************************************************************
// Calculate the values that rely on the histogram:
// (all three from Umbaugh, Wei et al. 1997)
//
// 1. Skew of the normalized intensity histogram 
// 2. Energy of the normalized intensity histogram
// 3. Entropy of the intensity histogram
//**************************************************************************
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::CalculateHistogramFeatures()
{
	if (!labelStatisticsFilter) return;

	typename LabelStatisticsType::HistogramPointer histo;
	TLPixel currentLabel;
	double vol, mean, sigma, t_skew, t_energy, t_entropy;
	double log, diff, cube, prob;
	
	for (int lab=0; lab<(int)labels.size(); ++lab)
	{
		currentLabel = labels[lab];
		if (currentLabel <= 0) continue;

		histo = labelStatisticsFilter->GetHistogram( currentLabel );
		vol = featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::VOLUME];
		mean = featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::MEAN];
		sigma = featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SIGMA];

		t_skew = 0;
		t_energy = 0;
		t_entropy = 0;
		
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
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SKEW] =  float(0);
		else
			featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::SKEW] =  float( t_skew / ( sigma * sigma * sigma ) );
			
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::ENERGY] = float( t_energy );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::ENTROPY] = float( -1 * t_entropy );
	}
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
	
	LabelImagePointer tempIntensityImage;

	typedef itk::RescaleIntensityImageFilter< IntensityImageType, LabelImageType > RescaleFilterType;
	typedef typename RescaleFilterType::Pointer RescaleFilterPointer;
	RescaleFilterPointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->SetInput(intensityImage);
	rescaleFilter->Update();
	tempIntensityImage = rescaleFilter->GetOutput();
	
	typedef itk::Statistics::ScalarImageTextureCalculator< LabelImageType > TextureCalcType;
	typedef typename TextureCalcType::Pointer TextureCalcPointer;
	TextureCalcPointer textureCalculator = TextureCalcType::New();
	textureCalculator->SetInput(tempIntensityImage);
	textureCalculator->SetImageMask(labelImage);
	textureCalculator->SetPixelValueMinMax(0,255);
	textureCalculator->SetFastCalculations(true);
	
	TLPixel currentLabel;
	for (int lab=0; lab<(int)labels.size(); ++lab)
	{
		currentLabel = labels.at(lab);
		if ((int)currentLabel <= 0) continue;
		
		LabelImageType::RegionType region;
		LabelImageType::SizeType size;
		LabelImageType::IndexType index;

		int loc = 0;
		for (int dim=0; dim<VImageDimension; dim++)
		{
			index[dim] = featureVals[currentLabel].BoundingBox[loc];	//bbox min
			size[dim] = featureVals[currentLabel].BoundingBox[loc+1]- index[dim] + 1;	//bbox max - min + 1
			loc += 2;
		}

        region.SetSize(size);
        region.SetIndex(index);
        labelImage->SetRequestedRegion(region);
        tempIntensityImage->SetRequestedRegion(region);
		textureCalculator->SetInsidePixelValue( currentLabel );
		textureCalculator->Compute();
		
		/*
		FeatureValueVector * 	GetFeatureMeans ()
		FeatureValueVector * 	GetFeatureStandardDeviations ()
		const FeatureNameVector * 	GetRequestedFeatures ()
		SetRequestedFeatures (const FeatureNameVector *_arg)
		*/
		
		TextureCalcType::FeatureValueVector *tex = textureCalculator->GetFeatureMeans();
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENERGY] = float( tex->ElementAt(0) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::T_ENTROPY] = float( tex->ElementAt(1) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INVERSE_DIFFERENCE_MOMENT] = float( tex->ElementAt(2) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::INERTIA] = float( tex->ElementAt(3) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_SHADE] = float( tex->ElementAt(4) );
		featureVals[currentLabel].ScalarFeatures[IntrinsicFeatures::CLUSTER_PROMINENCE] = float( tex->ElementAt(5) );
	}
	
	tempIntensityImage = 0;
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
		
	int totalBound = 0;
	for (TLPixel i=0; i<focusLabel; ++i)
	{
		totalBound = totalBound + sharePix[focusLabel][i];
	}
	for (TLPixel i=focusLabel+1; i<=GetMaxLabel(); ++i)
	{
		totalBound = totalBound + sharePix[i][focusLabel];
	}

	int sharedBound;
	if (focusLabel > neighborLabel)
		sharedBound = sharePix[focusLabel][neighborLabel];
	else if (focusLabel < neighborLabel)
		sharedBound = sharePix[neighborLabel][focusLabel];

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

	for (TLPixel i=1; i<label; ++i)
	{
		if(sharePix[label][i] > 0)
			nbs.push_back(i);
	}
	for (TLPixel i=label+1; i<=GetMaxLabel(); ++i)
	{
		if(sharePix[i][label] > 0)
			nbs.push_back(i);
	}
	return nbs;
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetCentroid(TLPixel label)
{
	typename LabelGeometryType::LabelPointType c = labelGeometryFilter->GetCentroid( label );
	for (unsigned int i = 0; i < VImageDimension; ++i)
	{
		featureVals[label].Centroid[i] = float( c[i] );
	}
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetWeightedCentroid(TLPixel label)
{
	typename LabelGeometryType::LabelPointType c = labelGeometryFilter->GetWeightedCentroid( label );
	for (unsigned int i = 0; i < VImageDimension; ++i)
	{
		featureVals[label].WeightedCentroid[i] = float( c[i] );
	}
}
	
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension>
::GetAxisLength(TLPixel label)
{
	typename LabelGeometryType::AxesLengthType aL = labelGeometryFilter->GetAxesLength( label );
	for (unsigned int i = 0; i < VImageDimension; ++i)
	{
		featureVals[label].AxisLength[i] = float( aL[i] );
	}
}
	
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension >
void LabelImageToFeatures< TIPixel, TLPixel, VImageDimension >
::GetBoundingBox(TLPixel label)
{
	typename LabelGeometryType::BoundingBoxType bbox = labelGeometryFilter->GetBoundingBox( label );
	for (unsigned int i = 0; i < VImageDimension*2; ++i)
	{
		featureVals[label].BoundingBox[i] = float( bbox[i] );
	}
}


} //end namespace ftk
#endif
