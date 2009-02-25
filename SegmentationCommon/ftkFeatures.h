/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkFeatures_h
#define __ftkFeatures_h

#include <string>
#include <vector>

namespace ftk
{

typedef struct {
	std::string name;
	std::string units;
	std::string description;
} FeatureInfoType;

class IntrinsicFeatures
{
public:

	static const enum
	{
		VOLUME, INTEGRATED_INTENSITY, ECCENTRICITY, ELONGATION, ORIENTATION, BBOX_VOLUME, \
		SUM, MEAN, MEDIAN, MINIMUM, MAXIMUM, SIGMA, VARIANCE, \
		SURFACE_GRADIENT, INTERIOR_GRADIENT, SURFACE_INTENSITY, INTERIOR_INTENSITY, \
		INTENSITY_RATIO, RADIUS_VARIATION, SURFACE_AREA, SHAPE, SHARED_BOUNDARY, \
		SOLIDITY, SKEW, ENERGY, ENTROPY 
	};	//FEATURES WILL GET ASSIGNED INT 0,1,...N-1

	static const int N = ENTROPY + 1;	//This is the number of scalar intrinsic features

	int Dimensions;

	float ScalarFeatures[N];

	std::vector<float> Centroid;
	std::vector<float> WeightedCentroid;
	std::vector<float> AxisLength;
	std::vector<float> BoundingBox;

	static FeatureInfoType Info[N];
};


//*****************************************************************************************
// THIS MAP HOLDS ALL FEATURE INFORMATION FOR ONE OBJECT
//*****************************************************************************************
//typedef std::map< IntrinsicFeaturesType, FeatureInfoType > FeatureInfoMapType;

}  // end namespace ftk

#endif