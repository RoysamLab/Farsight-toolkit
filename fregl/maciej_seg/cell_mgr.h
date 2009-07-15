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

/** @file cell_mgr.h
*   @brief cell class that calc all features
*
*   @author Maciej Wotjon
*/

#define _USE_MATH_DEFINES	//for math.h to define constants

#ifndef __CELL_MGR_H_
#define __CELL_MGR_H_

#include "pix_t.h"
#include "cell.h"
#include "compare.h"
#include "filter.h"	//for helper class to get train data

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"

/** @brief main cell class
*/
class cell_mgr: public filter
{
public:
	static const unsigned int msDim = 3;
	typedef itk::Image< float, msDim > mImageType;
	typedef itk::ImageRegionConstIterator< mImageType > mCitType;
	typedef itk::ImageRegionIterator< mImageType > mItType;
	//image for labeled output
	mImageType::Pointer mImgpIntensity;
	mImageType::Pointer mImgpGradient;
	mImageType::Pointer mImgpFilt;
	mImageType::Pointer mImgpLabel;
	mImageType::Pointer mImgpText;

	mImageType::IndexType miIndex;

	int mImageSize[3];
	std::vector<pix_t*> mvPix;	//pix of image
	std::vector<cell*> mvCell;	//cells
	std::vector<cell*> mvOrigCell; //unmerged cells from segmentation
	std::vector<cell*> mvTrain;	//training data sample 
	int mMaxLabel;

	vnl_matrix<double> mMtrxGlia;	//train data for glia
	vnl_matrix<double> mMtrxNeuron;	//train data for neuron
	vnl_matrix<double> mNoclass;

	vnl_matrix<double> mMeanN;
	vnl_matrix<double> mMeanG;

	double mProbGlia;	//priori prob
	double mProbNeuron;

	double mDetGlia;	//determinant of glia mtrx
	double mDetNeuron;	//determinant of neuron mtrx
	double mDetTrain;	//det for train cells
	double mDetAll;

	vnl_matrix<double> mMtrxCovGlia;	//covariance mtrx of glia
	vnl_matrix<double> mMtrxCovNeuron;	//covariance mtrx of neuron
	vnl_matrix<double> mMtrxCovTrain;	//cov of train data
	vnl_matrix<double> mMtrxCovAll;

	vnl_matrix<double> mMtrxInvCovGlia;	//inv cov mtrx
	vnl_matrix<double> mMtrxInvCovNeuron;
	vnl_matrix<double> mMtrxInvCovTrain;
	vnl_matrix<double> mMtrxInvCovAll;
	
	vnl_matrix<double> mMtrxTrain;		//train data matrix

	double mForeScore;
	double mBackScore;

	std::vector<int> bla;



	// Copyright 2000, softSurfer (www.softsurfer.com)
	// This code may be freely used and modified for any purpose
	// providing that this copyright notice is included with it.
	// SoftSurfer makes no warranty for this code, and cannot be held
	// liable for any real or imagined damage resulting from its use.
	// Users of this code must verify correctness for their application.
	// isLeft(): tests if a point is Left|On|Right of an infinite line.
	//    Input:  three points P0, P1, and P2
	//    Return: >0 for P2 left of the line through P0 and P1
	//            =0 for P2 on the line
	//            <0 for P2 right of the line
	//    See: the January 2001 Algorithm on Area of Triangles
	inline float
		is_left( pix_t *P0, pix_t *P1, pix_t *P2 )
	{
		return (P1->x_ - P0->x_)*(P2->y_ - P0->y_) - (P2->x_ - P0->x_)*(P1->y_ - P0->y_);
	}

public:	
	cell_mgr(){}
	cell_mgr(int ims[],mImageType::Pointer filt,mImageType::Pointer intensity,mImageType::Pointer gradient);
	~cell_mgr(void);

	void setup_cells(const short &rMin, const short &rMax, const short &rStep, std::string &file);

	mImageType::Pointer getlabelimage(void);

	void add_pix_vol_2d(pix_t* p,cell *c);

	void add_pix_bound_2d(pix_t* p,cell *c);

	void add_pix_vol_3d(pix_t* p,cell *c);

	void add_pix_bound_3d(pix_t* p,cell *c);

	void update_center(cell *c);

	cell* create_merged_cells(const std::vector<int> &v, std::vector<cell*> &mCell , cell *c2=NULL);

	void update_cells(std::vector<int> v,cell *c,bool del_cells=false);

	void relabel_cells(bool bRelabel = false);

	void build_rag_train(void);

	void update_nbrs(cell *c,const std::vector<int> &v);

	void merge_cells(void);

	bool compare_path(const std::vector<int> &p,const std::vector<int> &v);

	bool compare_path_to_tree(const std::vector<std::vector<int>* > &ptree,const std::vector<int> *pt);

	void score_cell_train(cell *c,std::vector<int> v);

	void update_convexity(cell *c);

	std::vector<pix_t*> chain_hull_2D(cell *c,std::vector<pix_t*> p);

	void update_depth(cell *c);

	void update_volume(cell *c);

	void get_data( std::string fileglia , std::string fileneuron );

	void update_avg_int(cell *c);

	std::vector<cell*> get_cells(void);

	void update_shape_fact(cell *c);

	void update_bend_eng(cell *c);

	void set_text_img(mImageType::Pointer t);

	void update_texture(cell *c);

	vnl_matrix<double> vnl_covariance_mtrx(vnl_matrix<double> m);

	void get_chain_dir(int& rx, int& ry, int d);

	void clear(void);

	void output_cells(std::string sN);

	void update_vol_grad_var(cell *c);

	void update_avg_bound_int(cell *c);

	void update_vol_grad(cell *c);

	void update_bound_grad(cell *c);

	double fast_ang(const int& crX1,const int& crY1,const int& crX2,const int& crY2);

	void merge_test_cells(void);

	void update_eccentricity( cell *c , const std::vector<int> &v);

	void update_ints_ratio( cell *c );
	
	void update_per_nbr( cell *c , const std::vector<int> &v );

	void score_background(std::string param);

	void score_foreground(void);
	
	void calc_train_data( const double &conv , const double &bend , const double &shape , const double &pnbr , const double &min_vol,const double &max_vol );

	void sep_glia_neur(const double &param_text , const double &param_int);

	void clear_train_data(void);
	
	void calc_more_train_data(const double &conv , const double &bend , const double &shape , const double &pnbr , const double &min_vol,const double &max_vol );

	bool check_shape( cell* c , const double &conv , const double &bend , const double &shape , const double &pnbr , const double &vol );

	bool check_shape_dist( cell* c , const double &conv , const double &bend , const double &shape , const double &pnbr , const double &vol );

	void classify_cells(void);

	void build_rag(const int &mode,std::string &file,const int &rDepth=5);

	void calc_ng_data(const bool &calcMtrx=true);

	void score_cell(cell *c,std::vector<int> v);

	void relabel_orig_cells(void);

	double get_perc_shared( cell *c1 , cell *c2 );

	void calc_mean_vect(const double &param=1.0);

	void check_all_cells(const double &conv , const double &bend , const double &shape , const double &pnbr , const	double &min_vol );

	void build_rag_mah(void);

	void score_cell_mah(cell *c,std::vector<int> v);

	void check_score_cells( const double &per );

	void score_foreground_mah( std::string param );

	void update_bound_grad_var(cell *c);

	inline double sigmoid( const double &a, const double &b, const double &net );

	void del_some_cells(void);

	void save_labeled(std::string &file);

	void project_cells(void);
	
	void remove_border_cells(void);

	void save_train_data(void);

	void remove_non_median_cells(void);
};

///**	@brief	helper class to gen correct train data
//*/
//class gen_train_data: public cell_mgr, public filter
//{
//public:
//	gen_train_data()
//	{
//	}
//	gen_train_data(int ims[],mImageType::Pointer filt,mImageType::Pointer intensity,mImageType::Pointer gradient)
//	{
//		//set image size
//		mImageSize[0] = ims[0];
//		mImageSize[1] = ims[1];
//		mImageSize[2] = ims[2];
//
//		//init pts
//		for( int z=0;z<mImageSize[2];++z )
//		{
//			for( int y=0;y<mImageSize[1];++y )
//			{
//				for( int x=0;x<mImageSize[0];++x )
//				{
//					pix_t* ptr=NULL;
//					ptr=new pix_t(x, y, z);
//					cell_mgr::mvPix.push_back(ptr);
//				}
//			}
//		}
//		//save filtered image
//		filter::mImgpFilt = mImageType::New();
//		filter::mImgpFilt = filt;
//
//		//save intensity image
//		mImgpIntensity = mImageType::New();
//		mImgpIntensity = intensity;
//
//		//save gradient image
//		cell_mgr::mImgpGradient = mImageType::New();
//		cell_mgr::mImgpGradient = gradient;
//
//		mMaxLabel=0;
//	}
//
//	/**	@brief	loads the labeled image from file and then sets up the segmentation
//	*/
//	void
//	setup_cells(const short &rMin, const short &rMax, const short &rStep)
//	{
//		data_rec data;
//		//get labeled image
//		mImgpLabel = data.get_image();
//		//sort labeled image to init cells
//		std::sort(cell_mgr::mvPix.begin(),cell_mgr::mvPix.end(),compare(mImgpLabel,mImageSize));
//
//		mCitType cit(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
//
//		//get nbr iterator
//		typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
//		NeighborhoodIteratorType::RadiusType radius;
//		radius[0] = 1;
//		radius[1] = 1;
//		radius[2] = 1;
//
//		NeighborhoodIteratorType it( radius, mImgpLabel , mImgpLabel->GetRequestedRegion() );
//
//		NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
//		NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
//		NeighborhoodIteratorType::OffsetType offset3 = {{0,1,0}};
//		NeighborhoodIteratorType::OffsetType offset4 = {{0,-1,0}};
//		NeighborhoodIteratorType::OffsetType offset5 = {{0,0,1}};
//		NeighborhoodIteratorType::OffsetType offset6 = {{0,0,-1}};
//
//		int istart;
//		//find first non background cell
//		for(istart=0;istart<cell_mgr::mvPix.size();istart++)
//		{
//			miIndex[0] = cell_mgr::mvPix[istart]->x_;
//			miIndex[1] = cell_mgr::mvPix[istart]->y_;
//			miIndex[2] = cell_mgr::mvPix[istart]->z_;
//			cit.SetIndex(miIndex);
//			if(cit.Get() != -1)
//				break;
//		}
//		//keep track of cell labels
//		int count = 1;
//		//label of previous pix, to know when we start pix of a new cell
//		int prev = -1;
//		cell *c=NULL;
//		for(int i=istart;i<cell_mgr::mvPix.size();i++)
//		{
//			miIndex[0] = cell_mgr::mvPix[i]->x_;
//			miIndex[1] = cell_mgr::mvPix[i]->y_;
//			miIndex[2] = cell_mgr::mvPix[i]->z_;
//			it.SetLocation(miIndex);
//			//check if new label found
//			if(prev != it.GetCenterPixel() )
//			{
//				//check to see if c has a cell and push it back on the vector
//				if(c != NULL)
//					cell_mgr::mvCell.push_back(c);
//				//create new cell
//				c = new cell(count);
//				//set old lebel as prev
//				prev = (int)it.GetCenterPixel();
//				count++;
//			}
//			//check 2d nbr labels, if all same as pix then add pix as volume otherwise as bound
//			if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
//				|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4) )
//			{
//				cell_mgr::add_pix_bound_2d(cell_mgr::mvPix[i],c);
//			}
//			else
//				cell_mgr::add_pix_vol_2d(cell_mgr::mvPix[i],c);
//
//			//check 3d nbr labels, if all same as pix then add pix as volume otherwise as bound
//			if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
//				|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4)
//				|| it.GetCenterPixel() != it.GetPixel(offset5) || it.GetCenterPixel() != it.GetPixel(offset6))
//			{
//				cell_mgr::add_pix_bound_3d(cell_mgr::mvPix[i],c);
//			}
//			else
//				cell_mgr::add_pix_vol_3d(cell_mgr::mvPix[i],c);
//		}
//		//add in last cell
//		if(c != NULL)
//			cell_mgr::mvCell.push_back(c);
//	}
//};

#endif
