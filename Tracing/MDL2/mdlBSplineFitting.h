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
#include "mdlUtils.h"

#include <vector>
#include <iostream>
#include <map>
#include <list>
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
	void SetOrder(unsigned int order)
	 {// now the itk B-spline only support the order 0,1,2,3,
	   splineOrder = order;
	   if (splineOrder >3) splineOrder = 3;
	   if (splineOrder <0) splineOrder = 0;
	 };
	void SetLevels(unsigned int levels)
	 { 
	   // typical value 4,5,6,7,8
	   splineLevels = levels;
	   if (splineLevels >8) splineLevels =8;
	   if (splineLevels <4) splineLevels =4;
	 };

	//Input:
	void SetNodes( std::vector<fPoint3D> * nds ){ nodes = nds; };
	void SetBBPairs( std::vector<pairE> * bbp ){ bbpairs = bbp; };
	//void SetSpinePairs( std::vector<pairE> *spp ) {spnpairs = spp };

	//Processing:
	bool Update();

	//Output:
	std::vector<fPoint3D> GetNodes(){ return nodes_out; };
	std::vector<pairE> GetBBPairs(){ return bbpairs_out; };
	//std::vector<pairE> GetSpinePairs(){ return spine_out; };


private:
	typedef std::set<int> SetType;
	typedef std::list<int> ListType;
	typedef std::map<int, ListType> ListMapType;
	typedef std::map<int, int> IntMapType;
	struct Graphprop{ int deg; SetType outVert; };

	bool debug;
	unsigned int splineOrder;
	unsigned int splineLevels;

	//Images & size
	ImageType::Pointer m_inputImage;
	ImageType::RegionType region;
	int sizeX;
	int sizeY;
	int sizeZ;
	long numPix;

	//Inputs:
	std::vector<fPoint3D> * nodes; //all nodes in skeleton
	std::vector<pairE> * bbpairs; //all pairs in backbone (initial)

	//The key is the branch number, the list contains all indexes in that branch
	ListMapType branches;			//Intermediate: holds branches info
	SetType branchPts;				//List of branch point ids (in nodes)

	//Outputs:
	std::vector<fPoint3D> nodes_out;
	std::vector<pairE> bbpairs_out;
	std::vector<pairE> spine_out;


	double dist2pts(fPoint3D p1, fPoint3D p2);
	int round(float number);
	int selffloor(float number);
	void findBranches();
	void smoothBranches();
	void detectExtraSpines();

	std::vector<fPoint3D> bbBSplineFitting(std::vector<fPoint3D> inPts, int numOut, unsigned int order, unsigned int levels);
};

}  // end namespace mdl

#endif
