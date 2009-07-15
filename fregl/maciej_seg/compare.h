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

/** @file compare.h
*   @brief helper compare class for sort function
*
*   @author Maciej Wotjon
*/

/** @brief helper compare class for sort function
*	need to be able to sort pix values from ITK image pointer, this allows for the pix class
*	,which only has the coordinates to be sorted according to image values
*/
#ifndef __COMPARE__H__
#define __COMPARE__H__

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "pix_t.h"

class compare
{
private:
	itk::Image< float, 3 >::Pointer im;
	int mImageSize[3];
	itk::Image< float, 3 >::IndexType ind;
public:
	compare(itk::Image< float, 3 >::Pointer i,int* ims)
	{
		im = itk::Image< float, 3 >::New();
		im = i;
		mImageSize[0] = ims[0];
		mImageSize[1] = ims[1];
		mImageSize[2] = ims[2];
	};
	~compare()
	{};
	bool operator()(const pix_t* p1,const pix_t* p2)
	{
		//compare pix values  by getting the iterator value at the miIndex
		itk::ImageRegionConstIterator< itk::Image< float, 3 > > it(im,im->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator< itk::Image< float, 3 > > it2(im,im->GetLargestPossibleRegion());
		ind[0] = p1->x_;
		ind[1] = p1->y_;
		ind[2] = p1->z_;
		it.SetIndex(ind);
		ind[0] = p2->x_;
		ind[1] = p2->y_;
		ind[2] = p2->z_;
		it2.SetIndex(ind);
		return it.Get()<it2.Get();
	}
};

/** @brief helper compare class for sort function for convex hull
*	sorts pix by increasing x and then y 
*/
class sort_xy
{
public:
	bool operator()(const pix_t* p1,const pix_t* p2)
	{
		if( p1->x_ < p2->x_ )
			return true;
		else if( p1->x_ == p2->x_ && p1->y_ < p2->y_ )
			return true;
		else
			return false;
	}
};
#endif
