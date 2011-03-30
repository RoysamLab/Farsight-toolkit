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
		INTENSITY_RATIO, CONVEXITY,RADIUS_VARIATION, SURFACE_AREA, SHAPE, SHARED_BOUNDARY, \
		SKEW, ENERGY, ENTROPY,
		T_ENERGY, T_ENTROPY, INVERSE_DIFFERENCE_MOMENT, INERTIA, CLUSTER_SHADE,CLUSTER_PROMINENCE
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

