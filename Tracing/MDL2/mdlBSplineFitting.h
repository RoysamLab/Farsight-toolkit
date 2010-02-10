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
/**************************************************************************  
// BSplineFitting
 *  Adapted Jan. 2010 by Isaac Abbott 
 *        
 *************************************************************************/
#ifndef __mdlBSplineFitting_h
#define __mdlBSplineFitting_h

#include "mdlTypes.h"

#include <vector>
#include <iostream>
#include <map>
#include <stdio.h>

#include "itkImage.h"
#include "itkPointSet.h"
#include "itkReview/itkBSplineScatteredDataPointSetToImageFilter.h"

namespace mdl
{

class BSplineFitting
{
public:
	BSplineFitting(ImageType::Pointer inImage);
	~BSplineFitting();

	//Setup:
	void SetDebug(bool inp = true){ debug = inp; };
	void SetOrder(int order){ splineOrder = order; };
	void SetLevels(int levels){ splineLevels = levels; };

	//Input:
	void SetNodes( std::vector<Point3D> * nds ){ nodes = nds; };
	void SetBBPairs( std::vector<pairE> * bbp ){ bbpairs = bbp; };

	//Processing:
	bool Update();

	//Output:

private:
	//Graph properties for each node: degree and connectivity
	struct Graphprop{ int deg; int outVert[8]; };

	bool debug;
	int splineOrder;
	int splineLevels;

	//Images & size
	ImageType::Pointer m_inputImage;
	ImageType::RegionType region;
	int sizeX;
	int sizeY;
	int sizeZ;
	long numPix;

	//Inputs:
	std::vector<Point3D> * nodes; //all nodes in skeleton
	std::vector<pairE> * bbpairs; //all pairs in backbone (initial)

	double dist2pts(Point3D p1, Point3D p2);
	int round(float number);
	int selffloor(float number);
	int BackboneBSplineFitting(Point3D * Points, int Num,
                           Point3D *SamplePoints, int Num2, 
						   int order, int levels);
};

}  // end namespace mdl

#endif
