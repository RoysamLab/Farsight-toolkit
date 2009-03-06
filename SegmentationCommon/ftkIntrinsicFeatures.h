/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkIntrinsicFeatures_h
#define __ftkIntrinsicFeatures_h

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

	enum ScalarTypes
	{
		VOLUME, INTEGRATED_INTENSITY, ECCENTRICITY, ELONGATION, ORIENTATION, BBOX_VOLUME, \
		SUM, MEAN, MEDIAN, MINIMUM, MAXIMUM, SIGMA, VARIANCE, \
		SURFACE_GRADIENT, INTERIOR_GRADIENT, SURFACE_INTENSITY, INTERIOR_INTENSITY, \
		INTENSITY_RATIO, RADIUS_VARIATION, SURFACE_AREA, SHAPE, SHARED_BOUNDARY, \
		SKEW, ENERGY, ENTROPY,
		T_ENERGY, T_ENTROPY, INVERSE_DIFFERENCE_MOMENT, INERTIA, CLUSTER_SHADE, CLUSTER_PROMINENCE
	};	//FEATURES WILL GET ASSIGNED INT 0,1,...N-1

	static const int N = CLUSTER_PROMINENCE + 1;	//This is the number of scalar intrinsic features

	int Dimensions;

	float ScalarFeatures[N];

	float Centroid[3];
	float WeightedCentroid[3];
	float AxisLength[3];
	float BoundingBox[6];		//min(X), max(X), min(Y), max(Y), min(Z), max(Z),...]

	static FeatureInfoType Info[N];

	void Print();

	//Tracking Stuff:
	int num;
	int time;
	int tag;
};


//*****************************************************************************************
// THIS MAP HOLDS ALL FEATURE INFORMATION FOR ONE OBJECT
//*****************************************************************************************
//typedef std::map< IntrinsicFeaturesType, FeatureInfoType > FeatureInfoMapType;

}  // end namespace ftk

#endif