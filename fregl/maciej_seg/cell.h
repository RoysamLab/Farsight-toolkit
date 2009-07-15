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

/** @file cell.h
*   @brief cell class
*   This is the class that represents the cells in the image
*
*   @author Maciej Wotjon
*/

#define _USE_MATH_DEFINES	//for math.h to define constants

#ifndef __CELL_H_
#define __CELL_H_

#include "pix_t.h"
#include <vector>

/** @brief main cell class, holds all the feature values of each cell
*/
class cell
{
public:
	int mVolume;		//volume
	double mAvgInt;	//avg intensity
	double mAvgBoundInt;	//avg bound intensity
	int mBoundPix;	//# of boundary pix
	int mCenterX;	//center coordinates
	int mCenterY;
	int mCenterZ;
	double mInsideGrad;		//avg inside gradient
	double mBoundGrad;		//avg boundary gradient
	int mTotInts;		//total intensity
	
	int mLabel;		//label
	double mBoundIntsRatio;		//boundary intensity ratio
	double mVarRad;		
	int mSharedPix;	//# of shered pix
	double mTexture;		//texture
	double mOrientation;		//orientation
	int mMajorAxis;		//major axis
	int mMinorAxis;		//minor axis
	double mEccentricity;		//eccencricity
	double mConvexity;		//convexity
	double mShapeFact;		//shape factor
	double mScore;		//score
	int mDepth;	//plane depth
	int mStartPlane;	//starting plane
	double	mBendEng;	//bending energy
	int mClass;		//class of cell
	double mVolGradVar;	//grad variance
	double mVolGrad;	//avg grad
	double mRadVar;	//radius variance
	double mPerNbr;	//percent pix nex to nbrs
	bool mIsMerged;	//is the cell merged
	bool mIsTrain;	//is the cell a trainig cell
	double mBoundGradVar;	//boudnary gradient variance
	double bendengstd;	//std of bending energy
	double bendengoutliers;	//number of bend enrg outliers, above 2*std
	int changebound;	//number of pix that are not in an ellipse
	
	std::vector<int> mvNbrs;		//cell nbrs
	std::vector<int> mvNbrsPix;		//amount of pix shared by cell with nbr
	std::vector<pix_t*> mvPix;		//all pix of cell
	std::vector<pix_t*> mvBoundPix2D;	//2d bound pix of cell
	std::vector<pix_t*> mvBoundPix3D;	//3d bound pix of cell
	std::vector<pix_t*> mvVolPix2D;		//2d vol pix of cell
	std::vector<pix_t*> mvVolPix3D;		//3d vol pix of cell

	cell(int label)
		: mVolume(-1)
		, mAvgInt(0.0)
		, mAvgBoundInt(0.0)
		, mBoundPix(-1)
		, mCenterX(-1)
		, mCenterY(-1)
		, mCenterZ(-1)
		, mInsideGrad(-1)
		, mBoundGrad(0.0)
		, mTotInts(-1)
		, mLabel(label)
		, mBoundIntsRatio(-1)
		, mVarRad(-1)
		, mSharedPix(-1)
		, mTexture(-1)
		, mOrientation(-1)
		, mMajorAxis(-1)
		, mMinorAxis(-1)
		, mEccentricity(0.0)
		, mConvexity(0.0)
		, mShapeFact(0.0)
		, mScore(0.0)
		, mDepth(-1)
		, mStartPlane(-1)
		, mBendEng(0.0)
		, mClass(-1)
		, mVolGradVar(0.0)
		, mVolGrad(0.0)
		, mRadVar(0.0)
		, mPerNbr(0.0)
		, mIsMerged(false)
		, mIsTrain(false)
		, mBoundGradVar(0.0)
		, bendengstd(0.0)
		, bendengoutliers(0.0)
		, changebound(0)
	{}
};

#endif
