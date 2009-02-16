/** @file wtd_lg.h
*   @brief main class for watershed algorithm
*
*   @author Maciej Wotjon
*/

#ifndef __WTD_LG_H_
#define __WTD_LG_H_

//#include "FarSight.h"
#include "cell.h"
#include "pix_t.h"
#include "compare.h"

#include <vector>
#include <functional>
#include <queue>
#include <cmath>
#include <cassert>
#include <iostream>

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


//----------------------------------------------------
//
// This program is written by Gang Lin.
// All rights reserved.
// Several routinues to realize the model-based 
// 3D watershed algorithm, application to
// nucleus cell segmentation, which are 
// called by id_imageL.
//
//----------------------------------------------------

#define MASK -2
#define WSHED 0

/**	@brief main watershed class
*/
class wtd_lg
{
private:
	
	int mImageSize[3];
	int id_image(int x,int y,int z) {return z*mImageSize[0]*mImageSize[1]+y*mImageSize[0]+x;}

	std::vector<pix_t*> mvPix;
	std::vector<unsigned char> mvDist;
	unsigned char nb6_flag;
	static const unsigned int msDim = 3;
	typedef itk::Image< float, msDim > mImageType;
	//image for labeled output
	mImageType::Pointer mImgpLabel;
	//filtered input image
	mImageType::Pointer mImgpFilt;
	typedef itk::ImageRegionConstIterator< mImageType > mCitType;
	typedef itk::ImageRegionIterator< mImageType > mItType;
	mImageType::IndexType miIndex;

public:
	wtd_lg(int ims[],std::vector<pix_t*> p,mImageType::Pointer ip);

	~wtd_lg();

	void run_watershed(short value_min,short value_max, short step_size, unsigned char set_nb6_flag);

	void get_neighbors( pix_t* p, std::vector<pix_t*>& nbs , bool tdflag=true);

	void remove_wtd(std::vector<pix_t*>& pix);

	short wtd3d(std::vector<pix_t*>& pix, short h_min, short h_max, short step_size);

	mImageType::Pointer getlabel(void);
};
#endif
