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

/** @file cell.cpp
*   @brief cell class
*   This is the class that represents the cells in the image
*
*   @author Maciej Wotjon
*/



//#include "FarSight.h" // For VS6 issues
#include "cell_mgr.h"
#include "pix_t.h"
#include "wtd_lg.h"
#include "compare.h"
#include "filter.h"
#include <vnl/algo/vnl_determinant.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>	//for accumulate

//rand
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/**	@brief	constructor
//	@param ims	size of image
//	@param filt	ITK pointer to filtered image
//	@param intensity ITK pointer to intensity image
//	@param gradient	ITK pointer to gradient image
*/
cell_mgr::cell_mgr(int ims[],mImageType::Pointer filt,mImageType::Pointer intensity,mImageType::Pointer gradient):
mDetGlia(0),	
mDetNeuron(0),
mDetTrain(0)
{
	//set image size
	mImageSize[0] = ims[0];
	mImageSize[1] = ims[1];
	mImageSize[2] = ims[2];

	//init pts
	for( int z=0;z<mImageSize[2];++z )
	{
		for( int y=0;y<mImageSize[1];++y )
		{
			for( int x=0;x<mImageSize[0];++x )
			{
				pix_t* ptr=NULL;
				ptr=new pix_t(x, y, z);
				mvPix.push_back(ptr);
			}
		}
	}

	//save filtered image
	mImgpFilt = mImageType::New();
	mImgpFilt = filt;

	//save intensity image
	mImgpIntensity = mImageType::New();
	mImgpIntensity = intensity;

	//save gradient image
	mImgpGradient = mImageType::New();
	mImgpGradient = gradient;

	mMaxLabel=0;
}

/**	@brief	destructor
//	cleans up
*/
cell_mgr::~cell_mgr(void)
{
	for(unsigned int i=0;i<mvPix.size();i++)
	{
		delete mvPix[i];
	}

	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] != NULL )
		{
			if( mvCell[i]->mIsMerged )	//only delete merged cells, rest will get deleted next
				delete mvCell[i];		//merged cells are not in mvOrigCell, also we cannot del all
			mvCell[i] = NULL;			//cells in mvCell b/c then ptrs in mvOrigCell would point to deallocated space
		}

	}

	for(unsigned int i=0;i<mvOrigCell.size();i++)
	{
		delete mvOrigCell[i];	//del all cells from orig seg
		mvOrigCell[i] = NULL;
	}
}

/**	@brief	runs the watershed algorithm
//	@param	rMin watershed param
//	@param	rMax watershed param
//	@param	rStep	step size of watershed flooding
//	@param	file	name of file to load labeled image
*/
void
cell_mgr::setup_cells(const short &rMin, const short &rMax, const short &rStep, std::string &file)
{
	std::ifstream in(std::string( "_label" + file ).c_str() );
	if( in.is_open() )
	{
		std::cout << "Loading labels" << std::endl;
		mImgpLabel = mImageType::New();
		mImgpLabel->SetRegions( mImgpIntensity->GetLargestPossibleRegion() );
		mImgpLabel->CopyInformation( mImgpIntensity );
		mImgpLabel->Allocate();

		mItType it( mImgpLabel,mImgpLabel->GetLargestPossibleRegion() );
		
		for( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{	
			int input=0;
			in >> input;
			it.Set(input);
		}
	}
	else
	{
		//init watershed
		wtd_lg w(mImageSize,mvPix,mImgpFilt);
		//run watershed
		w.run_watershed(rMin,rMax,rStep,1);
		mImgpLabel = mImageType::New();
		//get labeled image
		mImgpLabel = w.getlabel();
	}
	//sort labeled image to init cells
	std::sort(mvPix.begin(),mvPix.end(),compare(mImgpLabel,mImageSize));

	mCitType cit(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	//get nbr iterator
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 1;

	NeighborhoodIteratorType it( radius, mImgpLabel , mImgpLabel->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
	NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
	NeighborhoodIteratorType::OffsetType offset3 = {{0,1,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{0,-1,0}};
	NeighborhoodIteratorType::OffsetType offset5 = {{0,0,1}};
	NeighborhoodIteratorType::OffsetType offset6 = {{0,0,-1}};

	unsigned int istart;
	//find first non background cell
	for(istart=0;istart<mvPix.size();istart++)
	{
		miIndex[0] = mvPix[istart]->x_;
		miIndex[1] = mvPix[istart]->y_;
		miIndex[2] = mvPix[istart]->z_;
		cit.SetIndex(miIndex);
		if(cit.Get() != -1)
			break;
	}
	//keep track of cell labels
	int count = 1;
	//label of previous pix, to know when we start pix of a new cell
	int prev = -1;
	cell *c=NULL;
	for(unsigned int i=istart;i<mvPix.size();i++)
	{
		miIndex[0] = mvPix[i]->x_;
		miIndex[1] = mvPix[i]->y_;
		miIndex[2] = mvPix[i]->z_;
		it.SetLocation(miIndex);
		//check if new label found
		if(prev != it.GetCenterPixel() )
		{
			//check to see if c has a cell and push it back on the vector
			if(c != NULL)
				mvOrigCell.push_back(c);
			//create new cell
			c = new cell(count);
			//set old lebel as prev
			prev = (int)it.GetCenterPixel();
			count++;
		}

		if( mvPix[i]->x_ == mImageSize[0] - 1  || mvPix[i]->x_ == 0 || mvPix[i]->y_ == 0 || mvPix[i]->y_ == mImageSize[1] - 1 )
		{
			add_pix_bound_2d(mvPix[i],c);
			add_pix_vol_3d(mvPix[i],c);
			continue;
		}
		//check 2d nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
			|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4) )
		{
			add_pix_bound_2d(mvPix[i],c);
		}
		else
			add_pix_vol_2d(mvPix[i],c);

		//check 3d nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
			|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4)
			|| it.GetCenterPixel() != it.GetPixel(offset5) || it.GetCenterPixel() != it.GetPixel(offset6))
		{
			add_pix_bound_3d(mvPix[i],c);
		}
		else
			add_pix_vol_3d(mvPix[i],c);
	}
	//add in last cell
	if(c != NULL)
		mvOrigCell.push_back(c);

	for(unsigned int i=0;i<mvOrigCell.size();i++)
	{
		cell *temp = mvOrigCell[i];
		update_center(temp);
		mvCell.push_back(temp);
	}
	if( in.is_open() )
	{
		for(unsigned int i=0;i<mvCell.size(); i++ )
		{
			in >> mvCell[i]->mClass;
		}
		in.close();
	}

}

/**	@brief	get the segmented image
*	@return	ITK image pointer to segmented image
*/
cell_mgr::mImageType::Pointer
cell_mgr::getlabelimage(void)
{
	return mImgpLabel;
}

/** @brief add 2d volume pix to cell, vol on each plane
*	@param p added pix
*	@param c cell
*/
void 
cell_mgr::add_pix_vol_2d(pix_t* p,cell *c)
{
	c->mvPix.push_back(p);
	c->mvVolPix2D.push_back(p);
}

/** @brief add 2d boundary pix to cell, boundary on each plane
*	@param p added pix
*	@param c cell
*/
void 
cell_mgr::add_pix_bound_2d(pix_t* p,cell *c)
{
	c->mvPix.push_back(p);
	c->mvBoundPix2D.push_back(p);
}

/** @brief add 3d volume pix to cell
*	@param p added pix
*	@param c cell
*/
void 
cell_mgr::add_pix_vol_3d(pix_t* p,cell *c)
{
	c->mvVolPix3D.push_back(p);
}

/** @brief add 3d boundary pix to cell
*	@param p added pix
*	@param c cell
*/
void 
cell_mgr::add_pix_bound_3d(pix_t* p,cell *c)
{
	c->mvBoundPix3D.push_back(p);
}

/**	@brief update the center values for cell
*/
void cell_mgr::update_center(cell *c)
{
	//calc mean of each coordinate component

	c->mCenterX=0;
	c->mCenterY=0;
	c->mCenterZ=0;
	for(unsigned int i=0;i<c->mvPix.size();i++)
	{
		c->mCenterX+=c->mvPix[i]->x_;
		c->mCenterY+=c->mvPix[i]->y_;
		c->mCenterZ+=c->mvPix[i]->z_;
	}
	c->mCenterX/=c->mvPix.size();
	c->mCenterY/=c->mvPix.size();
	c->mCenterZ/=c->mvPix.size();
}

/**	@brief merges cells given by vector v
*	first the labels of all the cells are set as one label and then the cell is created and then labels are reset
*	@param	v	std vector of the cell labels to be merged
*/
cell*
cell_mgr::create_merged_cells(const std::vector<int> &v, std::vector<cell*> &mCell , cell *c2)
{
	if( v.empty() )	//if vector is empty someone did not want anything merged
		return NULL;

	if( v.size() == 1)	//return same cell if only one in vector
	{
		if( c2 != NULL )
			update_nbrs( c2 , bla );
		update_nbrs(mCell[v[0]-1] , bla);
		return mCell[v[0]-1];
	}

	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	cell* c = new cell(v[0]);

	//the first cell label will be used for the merged cells
	int label = v[0];

	//go through each cell that needs to be merged to add pix to new cell
	for(unsigned int i=0;i<v.size();i++)
	{
		if(mCell[v[i]-1] == NULL)	//no merge, cell already has been merged
		{
			delete c;
			return NULL;
		}
		//add pix of cell to be merged to merged cell
		c->mvPix.insert(c->mvPix.end(), mCell[v[i]-1]->mvPix.begin(),mCell[v[i]-1]->mvPix.end());
	}

	//go through each pix of the cell and relabel the cells as one label
	for(unsigned int j=0;j<c->mvPix.size();j++)
	{
		//set miIndex
		miIndex[0] = c->mvPix[j]->x_;
		miIndex[1] = c->mvPix[j]->y_;
		miIndex[2] = c->mvPix[j]->z_;
		it.SetIndex(miIndex);
		it.Set(label);
	}

	//create nbr iterator
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 1;

	NeighborhoodIteratorType nit( radius, mImgpLabel , mImgpLabel->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{0,-1,0}};
	NeighborhoodIteratorType::OffsetType offset2 = {{0,1,0}};
	NeighborhoodIteratorType::OffsetType offset3 = {{-1,0,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{1,0,0}};
	NeighborhoodIteratorType::OffsetType offset5 = {{0,0,1}};
	NeighborhoodIteratorType::OffsetType offset6 = {{0,0,-1}};

	for(unsigned int i=0;i<c->mvPix.size();i++)
	{
		miIndex[0] = c->mvPix[i]->x_;
		miIndex[1] = c->mvPix[i]->y_;
		miIndex[2] = c->mvPix[i]->z_;
		nit.SetLocation(miIndex);
		//check nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( nit.GetCenterPixel() != nit.GetPixel(offset1) || nit.GetCenterPixel() != nit.GetPixel(offset2)
			|| nit.GetCenterPixel() != nit.GetPixel(offset3) || nit.GetCenterPixel() != nit.GetPixel(offset4) )
		{
			c->mvBoundPix2D.push_back(c->mvPix[i]);
		}
		else
			c->mvVolPix2D.push_back(c->mvPix[i]);

		//check 3d nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( nit.GetCenterPixel() != nit.GetPixel(offset1) || nit.GetCenterPixel() != nit.GetPixel(offset2)
			|| nit.GetCenterPixel() != nit.GetPixel(offset3) || nit.GetCenterPixel() != nit.GetPixel(offset4)
			|| nit.GetCenterPixel() != nit.GetPixel(offset5) || nit.GetCenterPixel() != nit.GetPixel(offset6))
		{
			c->mvBoundPix3D.push_back(c->mvPix[i]);
		}
		else
			c->mvVolPix3D.push_back(c->mvPix[i]);
	}

	update_nbrs( c , bla);
	if( c2 != NULL )
		update_nbrs( c2 , bla);


	for(unsigned int i=0;i<v.size();i++)
	{
		//go through each pix of the cell and relabel the cells to orig label
		for(unsigned int j=0;j<mCell[v[i]-1]->mvPix.size();j++)
		{
			//set miIndex
			miIndex[0] = mCell[v[i]-1]->mvPix[j]->x_;
			miIndex[1] = mCell[v[i]-1]->mvPix[j]->y_;
			miIndex[2] = mCell[v[i]-1]->mvPix[j]->z_;
			it.SetIndex(miIndex);
			it.Set(mCell[v[i]-1]->mLabel);
		}
	}

	update_center(c);
	return c;
}

/**	@brief	save the merged cell into first merged cell position and del the rest
*	@param v vector of ints of cells to be merged
*	@param c merged cell
*	@param del_cells flag if cells need to be deleted (defualt false)
*/
void
cell_mgr::update_cells(std::vector<int> v,cell *c,bool del_cells)
{
	if( c == NULL || v.empty() || v.size() == 1)	//check to see if cell is NULL, vector to merge is empty or 										
		return;										//merge only one cell
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	//save merged cell into the position
	mvCell[v[0]-1] = c;
	//relabel the image
	for(unsigned int i=0;i<c->mvPix.size();i++)
	{
		miIndex[0] = c->mvPix[i]->x_;
		miIndex[1] = c->mvPix[i]->y_;
		miIndex[2] = c->mvPix[i]->z_;
		it.SetIndex(miIndex);
		it.Set(c->mLabel);
	}
	//remove the merged cells
	for(unsigned int i=1;i<v.size();i++)
	{
		mvCell[v[i]-1] = NULL;
	}

}

/**	@brief	relabel cells and image and erase NULL cells
*/
void
cell_mgr::relabel_cells(bool bRelabel)
{
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	for(it.GoToBegin();!it.IsAtEnd();++it)
	{
		it.Set( -1 );
	}

	//relabel the cells and the image
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if(mvCell[i] == NULL)
		{
			mvCell.erase(mvCell.begin() + i);
			i--;
		}
		else if(mvCell[i]->mLabel != (int)(i+1) )	//check if label needs to be updated
		{
			if( bRelabel == true )
				mvCell[i]->mLabel = i+1;
			for(unsigned int j=0;j<mvCell[i]->mvPix.size();j++)
			{
				miIndex[0] = mvCell[i]->mvPix[j]->x_;
				miIndex[1] = mvCell[i]->mvPix[j]->y_;
				miIndex[2] = mvCell[i]->mvPix[j]->z_;
				it.SetIndex(miIndex);
				it.Set( mvCell[i]->mLabel );
			}
			update_center(mvCell[i]);
		}
	}
}

/**	@brief	merges the cells
*/
void
cell_mgr::build_rag_train(void)
{
	//merge tree for each rag
	std::vector<std::vector<int>* > ttree;
	//temp pointer
	std::vector<int> *pt = NULL;
	//front of queue
	std::vector<int>* qfront = NULL;
	std::ofstream out("_mergecellstrain.txt");

	cell *c;
	for(unsigned int i=0;i<mvCell.size();i++)
	{		
		std::vector<std::vector<int>* > ptree;	//paths for current rag
		cell * pMxScC = NULL;	//ptr to cell with max score
		std::queue<std::vector<int> *> q;
		double max_score = 1E100 ;	//max score of merged cells, init to current cell score

		std::cout << "On Cell: " << i+1 << " of " << mvCell.size() << std::endl;

		//if cell is already merged move on to next cell
		if(mvCell[i] == NULL)
		{
			out << "skipping: " << i+1 << std::endl << std::endl;
			continue;
		}

		//queue
		pt = new std::vector<int>();

		//create init path
		pt->push_back(mvCell[i]->mLabel);

		//add path to queue
		q.push(pt);
		ttree.push_back(pt);

		//update_volume(mvCell[i]);	//for max rag length

		while(!q.empty() && (*q.front()).size() < 5 )
		{
			//get the cell from top of queue
			qfront = q.front();
			//delete top cell
			q.pop();

			//get the cell at end of path
			c = mvCell[qfront->at( qfront->size() - 1) - 1];

			if( c == NULL )	//if call already merged move on
			{
				//(*qfront).swap(std::vector<int>());
				//delete qfront;
				//qfront = NULL;
				continue;
			}

			update_nbrs(c,bla);
			//std::cout << "Updated nbrs\n";

			if(c->mvNbrs.empty())	//if a loner cell score it and move on
			{
				out << "loner: " << i+1 << std::endl << std::endl;
				score_cell_train( mvCell[i],bla);
				continue;
			}
			//check each nbr
			for(unsigned int j=0;j<c->mvNbrs.size();j++)
			{
				//check to make sure that the path is not going in a loop and check if the new cell has not been merged before
				if( std::find( qfront->begin(),qfront->end(),c->mvNbrs[j]) != qfront->end() || mvCell[c->mvNbrs[j]-1] == NULL)
				{
					continue;
				}	

				cell *tcell = create_merged_cells( *qfront ,mvCell,mvCell[c->mvNbrs[j]-1]);
				if( tcell == NULL )
				{ 
					continue;
				}
				update_avg_int( tcell );
				update_avg_int( mvCell[c->mvNbrs[j]-1] );
				update_texture( tcell );
				update_texture( mvCell[c->mvNbrs[j]-1] );
				out << "Trying: ";
				for(unsigned int k=0;k<(*qfront).size();k++)
				{
					out << (*qfront).at(k) << ' ';
				}
				out << "AND " << c->mvNbrs[j] << std::endl;
				//check if int and text are fairly similar for new segmented to be added to be merged
				if( tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() > 8000 || sqrt( pow(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt,2) + pow(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture,2) ) > 35 || get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) < 0.3 )
				
				//if( tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() > 8000 || abs(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) > 15|| abs(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) > 15 || get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1])) - 1 ) < .4 )
				{
					out << "Cell Text and Int not Good Enough: " <<  mvCell[c->mvNbrs[j]-1]->mLabel << ' ' << abs(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) << ' ' << abs(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) << ' ' << tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() << ' ' << get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) << ' ' << std::endl;

					if( (*qfront).size() != 1 )	//if its a merged cell delete it
						delete tcell;
					continue;
				}
				out << "Cell Text and Int GOOD: " <<  mvCell[c->mvNbrs[j]-1]->mLabel << ' ' << abs(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) << ' ' << abs(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) << ' ' << tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() << ' ' << get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) << std::endl;

				//add the nbr to the path
				pt = new std::vector<int>(*qfront);
				pt->push_back(c->mvNbrs[j]);
				//check if the path is already in tree
				if(!compare_path_to_tree(ptree,pt))
				{
					cell *tcell = create_merged_cells( *pt ,mvCell,mvCell[c->mvNbrs[j]-1]);
					if( tcell != NULL )
					{
						if( tcell->mvPix.size() < 8000 )
						{
							//add new path
							ptree.push_back(pt);
						}
					}
					delete tcell;
				}
				ttree.push_back(pt);
				q.push(pt);
			}
		}
		if(ptree.empty())	//if there is nothing to merge move on
		{
			out << "Nothing to merge" << std::endl;
			score_cell_train( mvCell[i],bla );
			continue;
		}

		int max_miIndex = -1;
		
		out << "trying: " << (*ptree[0]).at(0) << ' ' << max_score << std::endl;

		//create new cells with each new path
		for(unsigned int j=0;j<ptree.size();j++)
		{
			out << "trying to merge ";
			for(unsigned int k=0;k<(*ptree[j]).size();k++)
			{
				out << (*ptree[j]).at(k) << ' ';
			}
			c = create_merged_cells(*ptree[j],mvCell);
			//score each cell
			score_cell_train(c,*ptree[j]);
			out << ' ' << c->mScore;

			//check if new cell has a higher score than n
			if(c->mScore < max_score)
			{
				max_score = c->mScore;	//if yes set new max score
				max_miIndex = j;	//set max index	
				delete pMxScC;	//delete prev max score cell
				pMxScC = c;	//set ptr to new max cell
				c = NULL;
			}
			else
			{
				delete c;
				c = NULL;
			}
			out << ' ' << max_score << ' ' << max_miIndex << std::endl;
		}

		//if( max_miIndex != -1 )
		//{
		//	for(int k=0;k<(*ptree[max_miIndex]).size();k++)
		//	{
		//		score_cell_train( mvCell[ (*ptree[max_miIndex]).at(k) - 1 ],bla );
		//		if( mvCell[ (*ptree[max_miIndex]).at(k) - 1 ]->mScore < max_score )
		//		{
		//			out << (*ptree[max_miIndex]).at(k) << " is better" << std::endl;
		//			max_miIndex = -1;
		//			delete pMxScC;
		//			pMxScC = NULL;
		//			break;
		//		}
		//	}
		//}

		//if there is a cell with higher score than current update cells
		if( max_miIndex != -1 )
		{
			out << "Merging ";
			for(unsigned int k=0;k<(*ptree[max_miIndex]).size();k++)
			{
				out << (*ptree[max_miIndex]).at(k) << ' ';
			}

			update_cells(*ptree[max_miIndex],pMxScC);
			pMxScC = NULL;
			mvCell[i]->mIsMerged = true;
			out << std::endl;
			out << mvCell[i]->mLabel << ' ' << mvCell[i]->mScore << std::endl;
			--i;
		}
		out << std::endl;
		
	}
	out.close();
	for(unsigned int i=0;i<ttree.size();i++)
		delete ttree[i];
}

/**	@brief	updates the cell nbrs
*	@param	c cell ptr
*	@param	
*/
void
cell_mgr::update_nbrs(cell *c , const std::vector<int> &v)
{
	c->mvNbrs.clear();
	c->mvNbrsPix.clear();

	//create nbr iterator
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;

	NeighborhoodIteratorType nit( radius, mImgpLabel , mImgpLabel->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{0,1}};
	NeighborhoodIteratorType::OffsetType offset2 = {{0,-1}};
	NeighborhoodIteratorType::OffsetType offset3 = {{-1,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{1,0}};

	//check each pix nbrhood
	for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)
	{
		miIndex[0] = c->mvBoundPix2D[i]->x_;
		miIndex[1] = c->mvBoundPix2D[i]->y_;
		miIndex[2] = c->mvBoundPix2D[i]->z_;
		nit.SetLocation(miIndex);

		std::vector<int>::iterator it = std::find(c->mvNbrs.begin(),c->mvNbrs.end(),nit.GetPixel(offset1));
		//see if each nbr is included in the rag already
		if( it == c->mvNbrs.end() && nit.GetPixel(offset1) != c->mLabel
			&& nit.GetPixel(offset1) != -1)
		{
			//update the nbr vector
			c->mvNbrs.push_back((int)nit.GetPixel(offset1));
			c->mvNbrsPix.push_back( 0 );
		}
		else if( it != c->mvNbrs.end() )
		{
			c->mvNbrsPix[ it - c->mvNbrs.begin() ]++;
		}

		it = std::find(c->mvNbrs.begin(),c->mvNbrs.end(),nit.GetPixel(offset2));
		if(it == c->mvNbrs.end() && nit.GetPixel(offset2) != c->mLabel
			&& nit.GetPixel(offset2) != -1)
		{
			//update the nbr vector
			c->mvNbrs.push_back( (int) nit.GetPixel(offset2));
			c->mvNbrsPix.push_back( 0 );
		}
		else if( it != c->mvNbrs.end() )
		{
			c->mvNbrsPix[ it - c->mvNbrs.begin() ]++;
		}

		it = std::find(c->mvNbrs.begin(),c->mvNbrs.end(),nit.GetPixel(offset3));
		if( it == c->mvNbrs.end() && nit.GetPixel(offset3) != c->mLabel
			&& nit.GetPixel(offset3) != -1)
		{
			//update the nbr vector
			c->mvNbrs.push_back( (int) nit.GetPixel(offset3));
			c->mvNbrsPix.push_back( 0 );
		}
		else if( it != c->mvNbrs.end() )
		{
			c->mvNbrsPix[ it - c->mvNbrs.begin() ]++;
		}

		it = std::find(c->mvNbrs.begin(),c->mvNbrs.end(),nit.GetPixel(offset4));
		if( it == c->mvNbrs.end() && nit.GetPixel(offset4) != c->mLabel
			&& nit.GetPixel(offset4) != -1)
		{
			//update the nbr vector
			c->mvNbrs.push_back( (int) nit.GetPixel(offset4));
			c->mvNbrsPix.push_back( 0 );
		}
		else if( it != c->mvNbrs.end() )
		{
			c->mvNbrsPix[ it - c->mvNbrs.begin() ]++;
		}
	}
}

/**	@brief	see if both vectors have the same elements but not in same order
*	@param	p vector
*	@param	v vector
*	@return	true if the vectors have some elements not in same order
*/
bool 
cell_mgr::compare_path(const std::vector<int> &p,const std::vector<int> &v)
{
	//if size different elements def not same
	if( v.size() != p.size())
		return false;
	for(unsigned int i=0;i<p.size();i++)
	{

		//if cant find one element in other vector return false
		if(std::find(v.begin(),v.end(),p[i]) == v.end() )
			return false;	
	}
	return true; //if the path is found
}

/**	@brief	checks if path is already included in whole tree
*	@param	ptree the vector with all the paths
*	@param	t path that needs to be checked
*	@return	true if path is already in tree
*/
bool
cell_mgr::compare_path_to_tree(const std::vector<std::vector<int>* > &ptree,const std::vector<int> *pt)
{
	//go throug all paths in tree
	for(unsigned int i=0;i<ptree.size();i++)
	{
		//compare each path
		if(compare_path(*ptree[i],*pt))
			return true;
	}
	return false;
}

/**	@breif scores cell
*	@param c cell to be scored
*/
void
cell_mgr::score_cell_train(cell *c,std::vector<int> v)
{
	if( c == NULL )
		return;
	vnl_matrix<double> x(mMtrxTrain.cols(),1);	//init vector
	vnl_matrix<double> t(mMtrxTrain.cols(),1);	//init vector

	update_convexity( c );
	x( 0 , 0 ) = 1 - c->mConvexity ;

	update_shape_fact( c );
	x( 1 , 0 ) = c->mShapeFact ;

	update_bend_eng( c );
	x( 2 , 0 ) = c->mBendEng ;

	x( 3 , 0 ) = c->changebound;

	//x( 3 , 0 ) = c->bendengstd ;

	//x( 3 , 0 ) = c->bendengoutliers ;

	update_per_nbr( c , v );
	x( 4 , 0 ) = c->mPerNbr;

	c->mScore = ( x.transpose()*mMtrxInvCovTrain*x ).get(0,0);
}

/**	@brief updates the chull value of cell, the ratio of volume to chull volume
*	@param c cell
*	uses a chain hull to get chull and then fills in the outline and calc the ratio 
*	of volume to the total volume of the chull
*/

void
cell_mgr::update_convexity(cell *c)
{
	c->mConvexity = 0.0;
	update_depth(c);
	//tot pix of chull
	double hull = 0.0;
	//go through each plane
	for(int d=c->mStartPlane;d<=c->mStartPlane+c->mDepth;d++)
	{
		std::vector<pix_t*> tch;
		//		std::cout << "Getting plane: " << d << std::endl ;
		//get all pix from current plane
		for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)
		{
			if(c->mvBoundPix2D[i]->z_ == d)
				tch.push_back(c->mvBoundPix2D[i]);
		}

		//convex hull pix
		std::vector<pix_t*> ch;
		//calc chull
		ch = chain_hull_2D( c , tch );

		//		std::cout << "Ran chainhull\n";

		//now we need to outline the chull
		//find min max x and y
		int minx = mImageSize[0], miny = mImageSize[1], maxx = 0, maxy = 0;

		//get bounds for temp image
		for(unsigned int i=0;i<ch.size();i++)
		{
			//compare each max min to curr pix and see if it is better match for max min
			if( ch[i]->x_ > maxx)
				maxx = ch[i]->x_;
			if( ch[i]->x_ < minx)
				minx = ch[i]->x_;
			if( ch[i]->y_ < miny)
				miny = ch[i]->y_;
			if( ch[i]->y_ > maxy)
				maxy = ch[i]->y_;
		}
		if( maxx == minx || maxy == miny )
			continue;
		//		std::cout << maxx << ' ' << maxy << ' ' << minx << ' ' << miny << ' ' << std::endl;
		//		std::cout << "Building convex image\n";
		//create a temp image
		std::vector<char> timg( (maxy-miny+1)*(maxx-minx+1) );

		for(int i=0;i<(maxy-miny+1)*(maxx-minx+1);i++)
		{
			timg[i] = '.';
		}
		//		std::cout << "initialized conv image\n";

		//now for each pair of chull pix draw a line between then
		double x = 0.0,y = 0.0,h = 0.0;
		//		std::cout << "size ch: " << ch.size() << std::endl;
		for(unsigned int i=0;i<tch.size()-1;i++)
		{
			timg[ (tch[i]->y_ - miny )*(maxx-minx+1) + tch[i]->x_-minx ] = 'x';
		}

		for(unsigned int i=0;i<ch.size()-1;i++)
		{
			//use basic trig to get the line, x and y are changes in x for every change in y
			x = ch[i+1]->x_ - ch[i]->x_;
			y = ch[i+1]->y_ - ch[i]->y_;
			//length of line
			h=sqrt( y*y + x*x );

			//std::cout << i << ' ' << x << ' ' << y << ' ' << h <<std::endl;
			//step through the line and add an x
			for(double r=0;r<int(h);r+= 0.5)
			{

				timg[int( ( ch[i]->y_ - miny + int(r*(y/h)) )*(maxx-minx+1) + ch[i]->x_-minx + int(r*(x/h)) )] = 'x';

			}
		}
		x = ch[ ch.size() - 1 ]->x_ - ch[ 0 ]->x_;
		y = ch[ ch.size() - 1 ]->y_ - ch[ 0 ]->y_;
		//length of line
		h=sqrt( y*y + x*x );

		//std::cout << "Connecte 1st and last\n";

		if( h > 0 )	//check if left most upper and lower pix is same
		{
			//step through the line and add an x
			for(double r=0; r <= h + 0.1; r += 0.5)
			{
				if( (unsigned int)( ( ch[ 0 ]->y_ - miny + int(r*(y/h)) )*(maxx-minx+1) + ch[ 0 ]->x_-minx + int(r*(x/h)) ) < timg.size() )
					timg[int( ( ch[ 0 ]->y_ - miny + int(r*(y/h)) )*(maxx-minx+1) + ch[ 0 ]->x_-minx + int(r*(x/h)) )] = 'x';
			}
		}

		//now we need to fill in the chull, for each line search from left and right until we find the left and right chull
		//then we just fill in between the two pix, if there is no left and right do not fill in the line and if left right 
		//are equal also do not fill anything
		for(int yy=0;yy<=maxy-miny;yy++)	//step though each line
		{
			bool found_left=false;
			int left=-1;
			int right=-1;
			bool found_right=false;
			for(int xx=0;xx<=maxx-minx;xx++)	//step through each column
			{
				if ( !found_left && timg[(maxx-minx+1)*yy + xx] == 'x' )	//if havent found left and found chull pix
				{
					found_left = true;	//set left
					left=xx;
				}
				if ( !found_right && timg[(maxx-minx+1)*yy + (maxx-minx-xx)] == 'x' )	//if havent found right and found chull pix
				{
					found_right = true;	//set right
					right=maxx-minx-xx;
				}
			}
			if(right == -1)		//if havent found right
				right = left;
			if(left == -1)		//if havent found left
				left=right;
			for(int xx=left;xx<=right;xx++)		//fill in between left and right
			{
				timg[(maxx-minx+1)*yy + xx] = 'x';
			}
		}
		//count all the chull volume pix
		for(int i=0;i<(maxy-miny+1)*(maxx-minx+1);i++)
		{
			if(timg[i] == 'x')
				hull++;
		}
		//		std::cout << "Filled in x\n";
	}
	//calc ratio of volume to hull
	update_volume(c);
	if( hull != 0 )
		c->mConvexity += double(c->mVolume) / hull;
}

/**	@brief	calc the chain hull
*	@param p	boundary points of cell
*	@return	the chull of the points
*	This alg first sort the input pix by first increasing x and then y, then
*	consecutive pix are added to stack and checked to make sure that they are not 
*	on the left of the previous two pix
*	This is an implementation of the Andrew's Monotone Chain Algorithm
*/

std::vector<pix_t*>
cell_mgr::chain_hull_2D(cell *c,std::vector<pix_t*> p )
{
	std::vector<pix_t*> thull;
	std::vector<pix_t*> hull;
	pix_t *prev;
	std::sort(p.begin(),p.end(),sort_xy());
	thull.push_back(p[0]);
	prev = p[0];
	//find all the min x pix
	for(unsigned int i=1;i<p.size();i++)
	{
		if( p[i]->x_ != prev->x_ )
		{
			prev = p[i];
			thull.push_back(p[i]);
		}
	}
	//Now we need to do the same for bottom part of hull, this time we get the max y at each x and need to go backwards to keep order of
	//convex hull
	thull.push_back(p[p.size() - 1]);
	prev = p[p.size() - 1];
	//find all the max x pix
	for(int i=p.size() - 2;i>=0;i--)
	{
		if( p[i]->x_ != prev->x_ )
		{
			prev = p[i];
			thull.push_back(p[i]);
		}
	}

	//hold the curr pix miIndex
	unsigned int count = 0;
	while( count < thull.size() )
	{
		//add pix to hull
		hull.push_back(thull[count]);
		//need 3 pix to determine direction
		if( hull.size() < 3 )
		{
			count++;
			continue;
		}
		//check if the pix is on left
		if( is_left( hull[ hull.size() - 3 ] , hull[ hull.size() - 2 ] , hull[ hull.size() - 1 ] ) < 0 )
		{
			//if yes delete the two top pix the last pix will be added and direction recalculated
			hull.pop_back();
			hull.pop_back();
		}
		else
		{
			//if pix is on right then continue to next pix
			count++;
		}
	}
	return hull;
}

/**	@brief	updates the starting plane and plane depth
*	@param	c	cell
*/
void
cell_mgr::update_depth(cell *c)
{
	int maxd = 0;
	int mind = mImageSize[2];
	for(unsigned int i=0;i<c->mvPix.size();i++)
	{
		if( c->mvPix[i]->z_ > maxd )
			maxd = c->mvPix[i]->z_;
		if( c->mvPix[i]->z_ < mind)
			mind = c->mvPix[i]->z_;
	}
	c->mDepth = maxd - mind;
	c->mStartPlane = mind;
}

/**	@brief update the cell volume
*	@param	c cell
*/
void
cell_mgr::update_volume(cell *c)
{
	c->mVolume = c->mvPix.size();
}

/**	@brief	load the training data
*	saves the training data into a vnl mtrx
*/
void
cell_mgr::get_data(std::string fileglia , std::string fileneuron)
{
	vnl_matrix<double> temp;
	//load train data
	//vcl_ifstream in_glia("Gang/glia.txt");
	//vcl_ifstream in_neuron("Gang/neuron.txt");
	vcl_ifstream in_glia(fileglia.c_str() );
	vcl_ifstream in_neuron(fileneuron.c_str() );

	//	vcl_ofstream out("out1.txt");

	//check if file exists
	if(!in_glia)
	{
		std::cerr << "Could not find file in_glia" << std::endl;
		throw "Could not find file in_neuron";

	}
	if(!in_neuron)
	{
		std::cerr << "Could not find file in_neuron" << std::endl;
		throw "Could not find file in_neuron";
	}

	//save priori prob
	mProbGlia = 0.5;	
	//save mtrx
	in_glia >> mMtrxGlia;
	std::cout << "Glia: " << mMtrxGlia.rows() << ' ' << mMtrxGlia.cols() << std::endl;

	//save priori prob
	mProbNeuron = 0.5;
	//save mtrx
	in_neuron >> mMtrxNeuron;
	std::cout << "Neuron: " << mMtrxNeuron.rows() << ' ' << mMtrxNeuron.cols() << std::endl;

	//calc cov mtrx
	mMtrxCovGlia = vnl_covariance_mtrx(mMtrxGlia);

	mMtrxCovNeuron = vnl_covariance_mtrx(mMtrxNeuron);	

	//calc determinant
	mDetGlia = vnl_determinant(mMtrxCovGlia);

	mDetNeuron = vnl_determinant(mMtrxCovNeuron);

	//calc inv
	mMtrxInvCovGlia = vnl_matrix_inverse<double> (mMtrxCovGlia);
	mMtrxInvCovNeuron = vnl_matrix_inverse<double> (mMtrxCovNeuron);

	//	out << "GLIA\nCov\n" << mMtrxCovGlia << "\n\n";
	//	out << "DET: " << mDetGlia << "\n\n";
	//	out << "INV: " << mMtrxInvCovGlia << "\n\n";
	//
	//	out << "NEURON\nCov\n" << mMtrxCovNeuron << "\n\n";
	//	out << "DET: " << mDetNeuron << "\n\n";
	//	out << "INV: " << mMtrxInvCovNeuron << "\n\n";
}

/**	@brief	update the average intensity of cell
*	@param	c	cell
*/
void
cell_mgr::update_avg_int(cell *c)
{
	c->mAvgInt = 0.0;
	if(c->mvVolPix2D.empty())
		return;
	mCitType cit(mImgpIntensity,mImgpIntensity->GetLargestPossibleRegion());
	for(unsigned int i=0;i<c->mvVolPix2D.size();i++)
	{
		miIndex[0] = c->mvVolPix2D[i]->x_;
		miIndex[1] = c->mvVolPix2D[i]->y_;
		miIndex[2] = c->mvVolPix2D[i]->z_;
		cit.SetIndex(miIndex);
		c->mAvgInt += cit.Get();
	}
	c->mAvgInt /= double( c->mvVolPix2D.size() );
}

/**	@brief	get cells
*	@return vector of cells
*/
std::vector<cell*>
cell_mgr::get_cells(void)
{
	return mvCell;
}

/**	@brief	update shape factor
*	@param c cell
*/
void
cell_mgr::update_shape_fact(cell *c)
{
	if(!c->mvVolPix2D.empty())
		c->mShapeFact = double( c->mvBoundPix2D.size()*c->mvBoundPix2D.size()*c->mvBoundPix2D.size() ) / (36.0 * M_PI * double(c->mvVolPix2D.size()*c->mvVolPix2D.size()) );
}

/**	@brief	update the bend energy of cell
*	@param	c cell
*	@param	v vector of the bound pts
*	//implementation of chain code alg
*/
void
cell_mgr::update_bend_eng(cell *c)
{
	std::vector<double> totbend;
	c->mBendEng = 0.0;
	update_depth(c);
	int count = 0;	//bound pix count
	//go through each plane
	for(int d=c->mStartPlane;d<=c->mStartPlane+c->mDepth;d++)
	{
		//find min max x and y
		int minx = mImageSize[0], miny = mImageSize[1], maxx = 0, maxy = 0;

		//coord for current, previous and next pix
		int startX = 0, startY = 0, currX = 0, currY = 0, prevX = 0, prevY = 0, nextX = 0, nextY = 0;

		//direction
		short dir = 5;

		double ang1 = 0.0, ang2 = 0.0;	//angles
		double be = 0.0;	//temo bend eng


		std::vector<pix_t*> vt;

		//get all pix from current plane
		for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)
		{
			if(c->mvBoundPix2D[i]->z_ == d)
				vt.push_back(c->mvBoundPix2D[i]);
		}
		//get bounds for temp image
		for(unsigned int i=0;i<vt.size();i++)
		{
			//compare each max min to curr pix and see if it is better match for max min
			if( vt[i]->x_ > maxx)
				maxx = vt[i]->x_;
			if( vt[i]->x_ < minx)
				minx = vt[i]->x_;
			if( vt[i]->y_ < miny)
				miny = vt[i]->y_;
			if( vt[i]->y_ > maxy)
				maxy = vt[i]->y_;
		}
		//create a temp image
		std::vector<char> pTimg( (maxy-miny+3)*(maxx-minx+3) );

		//init temp image with '-'
		for(int i=0;i<(maxy-miny+3)*(maxx-minx+3);i++)
		{
			pTimg[i] = '-';
		}


		//set all bound pix to '+'
		for(unsigned int i=0;i<vt.size();i++)
		{
			pTimg[(vt[i]->y_-miny+1)*(maxx-minx+3) + vt[i]->x_-minx+0+1] = '+';
		}	

		//find top most pix
		for(int y=1;y<maxy-miny+2;y++)
		{
			for(int x=1;x<maxx-minx+2;x++)
			{
				if( pTimg[(maxx-minx+3)*y + x] == '+' )
				{
					currX = startX = x;
					currY = startY = y;
					pTimg[(maxx-minx+3)*currY + currX] = 'C';
					goto jump1;
				}
			}
		}
jump1:
		//temp
		int tX = 0, tY = 0;
		//if nbr found
		bool found = false;
		//step through direction 5 - 0
		for(int i=0;i<4;i++)
		{
			//check if dir > 7
			if(dir > 7)
				dir-=8;
			tX = currX;
			tY = currY;
			//get coord of the pix with this dir
			get_chain_dir(tX,tY,dir);
			dir++;
			//if its a bound pix
			if( pTimg[(maxx-minx+3)*tY + tX] == '+' )
			{
				//update next
				nextX = tX;
				nextY = tY;
				pTimg[(maxx-minx+3)*nextY + nextX] = 'N';
				found = true;
				break;
			}
		}

		//if cant find bound pix return
		if(!found)
			continue;

		//inv dir
		dir += 4;

		found = false;
		for(int i=0;i<7;i++)
		{
			//check if dir > 7
			if(dir > 7)
				dir-=8;
			tX = nextX;
			tY = nextY;
			//get coord of the pix with this dir
			get_chain_dir(tX,tY,dir);
			dir++;
			//if its a bound pix
			if( pTimg[(maxx-minx+3)*tY + tX] == '+' )
			{
				//update
				prevX = currX;
				prevY = currY;
				pTimg[(maxx-minx+3)*prevY + prevX] = '+';
				currX = nextX;
				currY = nextY;
				pTimg[(maxx-minx+3)*currY + currX] = 'C';
				nextX = tX;
				nextY = tY;
				pTimg[(maxx-minx+3)*nextY + nextX] = 'N';
				//calc angl
				ang1 = atan( double(currY- prevY) / double(currX - prevX) );
				ang2 = atan( double(nextY - currY) / double(nextX - currX) );
				be += (ang2 - ang1) * (ang2 - ang1);
				if( is_left( &pix_t(prevX,prevY,0) , &pix_t(currX,currY,0) , &pix_t(nextX,nextY,0) ) )
					c->changebound++;
				totbend.push_back( be );
				count++;
				found = true;
				break;
			}
		}

		//if cant find bound pix return
		if(!found)
			continue;
		//while current pix isnt the staring pix
		while( currX != startX || currY != startY )
		{
			found = false;
			dir += 4;
			for(int i=0;i<7;i++)
			{
				//check if dir > 7
				if(dir > 7)
					dir-=8;
				tX = nextX;
				tY = nextY;
				//get coord of the pix with this dir
				get_chain_dir(tX,tY,dir);
				dir++;
				//if its a bound pix
				if( pTimg[(maxx-minx+3)*tY + tX] == '+' )
				{
					//update
					prevX = currX;
					prevY = currY;
					pTimg[(maxx-minx+3)*prevY + prevX] = '+';
					currX = nextX;
					currY = nextY;
					pTimg[(maxx-minx+3)*currY + currX] = 'C';
					nextX = tX;
					nextY = tY;
					pTimg[(maxx-minx+3)*nextY + nextX] = 'N';
					//calc angle
					ang1 = atan( double(currY- prevY) / double(currX - prevX) );
					ang2 = atan( double(nextY - currY) / double(nextX - currX) );
					be += (ang2 - ang1) * (ang2 - ang1);
					if( is_left( &pix_t(prevX,prevY,0) , &pix_t(currX,currY,0) , &pix_t(nextX,nextY,0) ) )
						c->changebound++;
					totbend.push_back( be );
					count++;
					found = true;
					break;
				}
			}
			//if cant find bound pix go to next plane
			if(!found)
				break;

		}
		c->mBendEng += be;
	}
	if( count > 0 )	
	{
		c->mBendEng /= double( count );	
		if( count == 1 )	//if only one pix, exit
			return;
		for(unsigned int i=0; i<totbend.size();i++)		//calc std
		{
			c->bendengstd += ( totbend[i] - c->mBendEng ) * ( totbend[i] - c->mBendEng );
		}
		c->bendengstd /= count - 1;
		c->bendengstd = sqrt( c->bendengstd );

		for(unsigned int i=0; i<totbend.size();i++)		//calc percent outliers
		{
			if( ( totbend[i] - c->mBendEng ) > ( 2*c->bendengstd ) )
				c->bendengoutliers++;
		}
	}
}


/**	@brief	save texture image
*	@param t	itk image pointer for text image
*/
void
cell_mgr::set_text_img(mImageType::Pointer t)
{
	//save texture image
	mImgpText = mImageType::New();
	mImgpText = t;
}


/**	@brief	update the texture
*	@param c cell
*/
void
cell_mgr::update_texture(cell *c)
{
	//calc avg text
	c->mTexture = 0.0;

	//if no volume save text as 0
	if( c->mvVolPix3D.empty() )
		return;

	//get nbr iterator
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;

	NeighborhoodIteratorType it( radius, mImgpIntensity , mImgpIntensity->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{-1,-1,0}};
	NeighborhoodIteratorType::OffsetType offset2 = {{1,-1,0}};
	NeighborhoodIteratorType::OffsetType offset3 = {{-1,1,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{1,1,0}};

	//calc avg text
	for(unsigned int i=0;i<c->mvVolPix2D.size();i++)
	{
		miIndex[0] = c->mvVolPix2D[i]->x_;
		miIndex[1] = c->mvVolPix2D[i]->y_;
		miIndex[2] = c->mvVolPix2D[i]->z_;
		it.SetLocation(miIndex);

		c->mTexture += sqrt( ( it.GetCenterPixel() - it.GetPixel(offset1) )*( it.GetCenterPixel() - it.GetPixel(offset1) )
			+ ( it.GetCenterPixel() - it.GetPixel(offset2) )*( it.GetCenterPixel() - it.GetPixel(offset2) )
			+ ( it.GetCenterPixel() - it.GetPixel(offset3) )*( it.GetCenterPixel() - it.GetPixel(offset3) )
			+ ( it.GetCenterPixel() - it.GetPixel(offset4) )*( it.GetCenterPixel() - it.GetPixel(offset4) ) );
	}
	c->mTexture /= double( c->mvVolPix3D.size() );
}

/**	@brief	calculates the covariance mtrx of the input mtrx
*	@param m input mtrx of type vnl_matrix<double>
*	@return the covariance mtrx of m
*/
vnl_matrix<double>
cell_mgr::vnl_covariance_mtrx(vnl_matrix<double> m)
{
	vnl_matrix<double> mean(m.rows(),m.columns(),0);
	for(unsigned int i=0;i<m.columns();i++)
	{
		mean.set_column(i,m.get_column(i).mean());
	}
	m = m - mean;
	vnl_matrix<double> cov(m.columns(),m.columns(),0);

	for(unsigned int i=0;i<m.columns();i++)
	{
		cov.set_row(i,m.transpose() * m.get_column(i));
	}
	cov /= m.rows() - 1;
	return cov;
}

/**	@brief	get the coord of pix that has dir and is next to rx and ry
*	@param rx x coord of curr pix
*	@param ry y coord of curr pix
*	@param d direction of next pix
*/
void
cell_mgr::get_chain_dir(int& rx, int& ry, int d)
{
	//check if dir is > 7 if so sub 8 to make it <= 7
	if(d > 7)
		d -= 8;
	switch ( d )
	{
	case 0:
		rx++;
		break;
	case 1:
		rx++;
		ry++;
		break;
	case 2:
		ry++;
		break;
	case 3:
		rx--;
		ry++;
		break;
	case 4:
		rx--;
		break;
	case 5:
		rx--;
		ry--;
		break;
	case 6:
		ry--;
		break;
	case 7:
		rx++;
		ry--;
		break;
	default:
		break;
	} 
}

/**	@brief	clear pts and cell vectors
*/
void
cell_mgr::clear(void)
{
	for(unsigned int i=0;i<mvPix.size();i++)
	{
		delete mvPix[i];
		mvPix[i] = NULL;
	}
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		delete mvCell[i];
		mvCell[i] = NULL;
	}
}


/**	@brief ouputs cell features to file
*	@param sN string name of file
*/
void
cell_mgr::output_cells(std::string sN)
{
	std::vector<int> v;
	std::ofstream outglia(sN.c_str());

	cell *pC;
	outglia << std::setw(13) << "Label" 
		<< std::setw(13) << "Volume" 
		<< std::setw(13) << "Intensity" 
		<< std::setw(13) << "Texture" 
		<< std::setw(13) << "Convexity" 
		<< std::setw(13) << "Shape Fact" 
		<< std::setw(13) << "Bend Eng" 
		<< std::setw(13) << "Bound Grad" 
		<< std::setw(13) << "Vol Grad" 
		<< std::setw(13) << "Ecc" 
		<< std::setw(13) << "rad var"
		<< std::setw(13) << "ints ratio" 
		<< std::setw(13) << "per nbr" 
		<< std::setw(13) << "score" 
		<< std::setw(13) << "Class" 
		<< std::setw(13) << "X" 
		<< std::setw(13) << "Y" 
		<< std::setw(13) << "Z" 
		<< std::setw(13) << "Nbr" 
		<< std::endl;

	for(unsigned int i=0;i<mvCell.size();i++)
	{
		//std::cout << i << std::endl;
		pC = mvCell[i];
		if( pC == NULL )
			continue;
		//if( pC->mClass == -1 || pC->mClass == -1)
		//	continue;

		//std::cout << "volume\n";
		update_volume(pC);
		update_avg_int(pC);
		update_texture(pC);
		update_convexity(pC);
		update_shape_fact(pC);
		update_bend_eng(pC);
		update_center(pC);
		update_vol_grad(pC);
		update_bound_grad(pC);
		update_eccentricity(pC,bla);
		update_ints_ratio(pC);
		update_per_nbr(pC,bla);
		update_nbrs(pC,bla);
		//score_cell(pC);

		outglia << std::setprecision(4);

		outglia << std::setw(13) << pC->mLabel 
			<< std::setw(13) << pC->mVolume 
			<< std::setw(13) << pC->mAvgInt 
			<< std::setw(13) << pC->mTexture 
			<< std::setw(13) << pC->mConvexity 
			<< std::setw(13) << pC->mShapeFact 
			<< std::setw(13) << pC->mBendEng 
			<< std::setw(13) << pC->mBoundGrad 
			<< std::setw(13) << pC->mVolGrad  
			<< std::setw(13) << pC->mEccentricity
			<< std::setw(13) << pC->mRadVar
			<< std::setw(13) << pC->mBoundIntsRatio
			<< std::setw(13) << pC->mPerNbr
			<< std::setw(13) << pC->mScore
			<< std::setw(13) << pC->mClass
			<< std::setw(13) << pC->mCenterX 
			<< std::setw(13) << pC->mCenterY 
			<< std::setw(13) << pC->mCenterZ ;

		for(unsigned int j=0; j< pC->mvNbrs.size(); j++)
		{
			outglia << std::setw(13) << pC->mvNbrs[j];
			outglia << std::setw(13) << (double) pC->mvNbrsPix[j] / (double) pC->mvBoundPix2D.size();
		}
		outglia << std::endl;
	}
}

/**	@brief	update avg vol gradient variance
*	@param	c cell
*/
void
cell_mgr::update_vol_grad_var(cell *c)
{
	c->mVolGradVar = 0.0;
	if(c->mvVolPix3D.empty())
		return;
	mItType it(mImgpGradient,mImgpGradient->GetLargestPossibleRegion());	//image

	update_depth(c);

	//go through each plane
	for(int d=c->mStartPlane;d<=c->mStartPlane+c->mDepth;d++)
	{
		double AvgGrad = 0.0;
		double CellCount = 0; //count cells in plane
		for(unsigned int i=0;i<c->mvVolPix3D.size();i++)	//mean
		{
			if( c->mvVolPix2D[i]->z_ == d)	//check if pix are in same plane
			{
				miIndex[0] = c->mvVolPix3D[i]->x_;
				miIndex[1] = c->mvVolPix3D[i]->y_;
				miIndex[2] = c->mvVolPix3D[i]->z_;
				it.SetIndex(miIndex);
				AvgGrad += it.Get();
				CellCount++;
			}
		}
		AvgGrad /= CellCount;

		for(unsigned int i=0;i<c->mvVolPix3D.size();i++)	//unbaised variance
		{
			if( c->mvVolPix2D[i]->z_ == d)
			{
				miIndex[0] = c->mvVolPix3D[i]->x_;
				miIndex[1] = c->mvVolPix3D[i]->y_;
				miIndex[2] = c->mvVolPix3D[i]->z_;
				it.SetIndex(miIndex);
				c->mVolGradVar += ( double( it.Get() ) - AvgGrad ) * ( double( it.Get() ) - AvgGrad ) / double( CellCount - 1.0 );
				//it.Set(AvgVar);
			}
		}
	}		
	c->mVolGradVar /= double( c->mvVolPix3D.size() );
}

/**	@brief	update the average boundary intensity of cell
*	@param	c	cell
*/
void
cell_mgr::update_avg_bound_int(cell *c)
{
	c->mAvgBoundInt = 0.0;
	mCitType cit(mImgpIntensity,mImgpIntensity->GetLargestPossibleRegion());
	for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)
	{
		miIndex[0] = c->mvBoundPix2D[i]->x_;
		miIndex[1] = c->mvBoundPix2D[i]->y_;
		miIndex[2] = c->mvBoundPix2D[i]->z_;
		cit.SetIndex(miIndex);
		c->mAvgBoundInt += cit.Get();
	}
	c->mAvgBoundInt /= double(c->mvBoundPix2D.size());
}

/**	@brief	update avg vol gradient
*	@param	c cell
*/
void
cell_mgr::update_vol_grad(cell *c)
{
	c->mVolGrad = 0.0;
	if(c->mvVolPix3D.empty())
		return;
	mItType it(mImgpGradient,mImgpGradient->GetLargestPossibleRegion());	//image
	for(unsigned int i=0;i<c->mvVolPix3D.size();i++)
	{
		miIndex[0] = c->mvVolPix3D[i]->x_;
		miIndex[1] = c->mvVolPix3D[i]->y_;
		miIndex[2] = c->mvVolPix3D[i]->z_;
		it.SetIndex(miIndex);
		c->mVolGrad += it.Get();
	}
	c->mVolGrad /= c->mvVolPix3D.size();
}

/**	@brief	update avg bound gradient
*	@param	c cell
*/
void
cell_mgr::update_bound_grad(cell *c)
{
	c->mBoundGrad = 0.0;
	if(c->mvBoundPix3D.empty())	//empty boun
		return;
	mItType it(mImgpGradient,mImgpGradient->GetLargestPossibleRegion());	//image

	for(unsigned int i=0;i<c->mvBoundPix3D.size();i++)	//save each value
	{
		miIndex[0] = c->mvBoundPix3D[i]->x_;
		miIndex[1] = c->mvBoundPix3D[i]->y_;
		miIndex[2] = c->mvBoundPix3D[i]->z_;
		it.SetIndex(miIndex);
		c->mBoundGrad += it.Get();
	}
	c->mBoundGrad /= double( c->mvBoundPix3D.size() );	//avg
}

/**	@brief	loads cells to be merged to create a training sample
*	input file should have the cell labels followed by a -1 at the end
*/

void
cell_mgr::merge_test_cells(void)
{
	std::vector<int> v;
	std::ifstream in("cells_to_merge.txt");
	while( !in.eof() )
	{
		int temp;
		in >> temp;
		if( temp != -1)		//get each cell to be merged from file
			v.push_back(temp);
		else
		{
			update_cells(v,create_merged_cells(v,mvCell));		//create the cell and merge cells
			mvCell[v[0]-1]->mClass = 1;
			v.clear();
		}
	}
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if(mvCell[i] != NULL)
			cell_mgr::update_center(mvCell[i]);		//make sure to update the center for
	}
}

/**	@brief	calcualtes the ecc of cell
*	alg works by starting in center of each plane and at 5 diff ang extend 
*	radius until reach bound of cell, the legnth is then saved and the smallest and
*	largest are used to calc the major and minor axis
*	@param c	cell
*/
void
cell_mgr::update_eccentricity( cell *c ,const std::vector<int> &v)
{
	c->mEccentricity = 0.0;
	c->mRadVar = 0.0;
	c->mMajorAxis = 0;
	c->mMinorAxis = 0;
	update_depth(c);	//need to have cell dpeth
	update_center(c);	//need to have current center
	mItType itU(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());	//image iter
	mItType itL(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	if( !v.empty() )
	{
		//go through each pix of the cell and relabel the cells as one label
		for(unsigned int j=0;j<c->mvPix.size();j++)
		{
			//set miIndex
			miIndex[0] = c->mvPix[j]->x_;
			miIndex[1] = c->mvPix[j]->y_;
			miIndex[2] = c->mvPix[j]->z_;
			itL.SetIndex(miIndex);
			itL.Set(c->mLabel);
		}
	}

	for(int d=c->mStartPlane;d<=c->mStartPlane+c->mDepth;d++)	//go through each plane and find major and minor axis
	{
		std::vector<double> radii;	//radii of the cell
		for( double ang=0; ang<=M_PI; ang+=M_PI/5.0 )	//iter through 10 angles
		{		
			bool goU = true,	//flags to cont upper and lower rad
				goL = true;

			double rU = 0.0		//radius to expand until edge of cell is found,need one for upper and lower parts
				, rL = 0.0;

			while( goU && goL )
			{		
				if( goU )
				{			
					miIndex[0] = (int) ( c->mCenterX + rU * cos(ang) );	//update the position of iter for top part of radius
					miIndex[1] = (int) ( c->mCenterY + rU * sin(ang) );
					miIndex[2] = d;
					if( miIndex[0] < mImageSize[0] && miIndex[1] < mImageSize[1] 
					&& miIndex[0] >= 0 && miIndex[1] >= 0)	//check if radius is out of bounds of image
						itU.SetIndex(miIndex);
					else
						goU = false;


					if( itU.Get() == c->mLabel )	//check if still in cell,if not end
						rU++;
					else
						goU = false;
				}

				if( goL )
				{		
					miIndex[0] = (int) ( c->mCenterX + rL * cos(ang) );	//update the position of iter for top part of radius
					miIndex[1] = (int) ( c->mCenterY + rL * sin(ang) );
					miIndex[2] = d;
					if( miIndex[0] < mImageSize[0] && miIndex[0] < mImageSize[1] 
					&& miIndex[0] >= 0 && miIndex[1] >= 0 )	//check if radius is out of bounds of image
						itL.SetIndex(miIndex);
					else
						goL = false;

					if( itL.Get() == c->mLabel )	//check if still in cell,if not end
						rL--;	
					else
						goL = false;
				}
			}
			if( rU - rL == 0 )	//need to add one to make sure no /0
				rU++;
			radii.push_back(rU - rL);	//add radius to vector
		}
		sort(radii.begin(),radii.end());	//sort radii


		if( radii[ radii.size() - 1 ] == 0 || radii[ 0 ] == 0 )	//check if cell too small to have axis
			continue;

		c->mMajorAxis += (int) radii[ radii.size() - 1 ];	//add largest radius
		c->mMinorAxis += (int) radii[ 0 ];	//add smallest radius

		c->mMajorAxis += (int) radii[ radii.size() - 1 ];	//add largest radius
		c->mMinorAxis += (int) radii[ 0 ];	//add smallest radius

		double mean = 0.0;
		for(unsigned int i=0; i<radii.size(); i++)
		{
			mean += radii[i]  ;
		}

		mean /= 4.0;
		double tRad = 0.0; //temp radius
		for(unsigned int i=0; i<radii.size(); i++)
		{
			tRad += ( radii[i] - mean ) * ( radii[i] - mean ) ;
		}
		if( tRad != 0.0 )
			tRad /= 4.0;
		c->mRadVar += tRad;
	}
	if(c->mMajorAxis == 0 || c->mMinorAxis == 0)	//return if no axis
		return;
	c->mMajorAxis /= c->mDepth + 1;		//take avg of all planes, need to add 1 to depth
	c->mMinorAxis /= c->mDepth + 1;
	c->mRadVar /= c->mDepth + 1;
	c->mRadVar = sqrt( c->mRadVar );
	c->mEccentricity = (double)c->mMajorAxis / (double)c->mMinorAxis;

	for(unsigned int i=0;i<v.size();i++)
	{
		//go through each pix of the cell and relabel the cells to orig label
		for(unsigned int j=0;j<mvCell[v[i]-1]->mvPix.size();j++)
		{
			//set miIndex
			miIndex[0] = mvCell[v[i]-1]->mvPix[j]->x_;
			miIndex[1] = mvCell[v[i]-1]->mvPix[j]->y_;
			miIndex[2] = mvCell[v[i]-1]->mvPix[j]->z_;
			itL.SetIndex(miIndex);
			itL.Set(mvCell[v[i]-1]->mLabel);
		}
	}
}


/**	@brief calc the bound to vol intensity ratio
*	@param c cell ptr
*/
void
cell_mgr::update_ints_ratio( cell *c )
{
	c->mBoundIntsRatio = 0.0;
	update_avg_int(c);
	update_avg_bound_int(c);
	if( c->mAvgInt == 0 )
		return;
	c->mBoundIntsRatio = c->mAvgBoundInt / c->mAvgInt ;
}

/**	@brief	update the percentage of pix that have a nbr
*	@param c cell ptr
*/	
void
cell_mgr::update_per_nbr( cell *c , const std::vector<int> &v )
{
	if( !v.empty() )
	{
		mItType itL(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
		//go through each pix of the cell and relabel the cells as one label
		for(unsigned int j=0;j<c->mvPix.size();j++)
		{
			//set miIndex
			miIndex[0] = c->mvPix[j]->x_;
			miIndex[1] = c->mvPix[j]->y_;
			miIndex[2] = c->mvPix[j]->z_;
			itL.SetIndex(miIndex);
			itL.Set(c->mLabel);
		}
	}
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;	//set kernel rad
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;

	NeighborhoodIteratorType nit( radius, mImgpLabel , mImgpLabel->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{0,1}};		//set offsets
	NeighborhoodIteratorType::OffsetType offset2 = {{0,-1}};
	NeighborhoodIteratorType::OffsetType offset3 = {{-1,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{1,0}};

	double nbrpix = 0.0; 	//nbr pix
	//check each pix nbrhood
	for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)
	{
		miIndex[0] = c->mvBoundPix2D[i]->x_;
		miIndex[1] = c->mvBoundPix2D[i]->y_;
		miIndex[2] = c->mvBoundPix2D[i]->z_;
		nit.SetLocation(miIndex);

		if(nit.GetPixel(offset1) != c->mLabel && nit.GetPixel(offset1) != -1)
		{
			++nbrpix;
			continue;
		}
		if(nit.GetPixel(offset2) != c->mLabel && nit.GetPixel(offset2) != -1)
		{
			++nbrpix;
			continue;
		}
		if(nit.GetPixel(offset3) != c->mLabel && nit.GetPixel(offset3) != -1)
		{
			++nbrpix;
			continue;
		}
		if(nit.GetPixel(offset4) != c->mLabel && nit.GetPixel(offset4) != -1)
		{
			++nbrpix;
			continue;
		}
	}
	c->mPerNbr = nbrpix / (double) c->mvBoundPix2D.size();

	for(unsigned int i=0;i<v.size();i++)
	{
		mItType itL(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
		//go through each pix of the cell and relabel the cells to orig label
		for(unsigned int j=0;j<mvCell[v[i]-1]->mvPix.size();j++)
		{
			//set miIndex
			miIndex[0] = mvCell[v[i]-1]->mvPix[j]->x_;
			miIndex[1] = mvCell[v[i]-1]->mvPix[j]->y_;
			miIndex[2] = mvCell[v[i]-1]->mvPix[j]->z_;
			itL.SetIndex(miIndex);
			itL.Set(mvCell[v[i]-1]->mLabel);
		}
	}
}

void 
cell_mgr::score_background(std::string param)
{
	param = "back_" + param + ".txt";
	std::ofstream out( param.c_str() );
	std::ofstream oscore( "_score.txt",std::ofstream::app );
	mItType itL(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());	//label image iter
	mItType itI(mImgpIntensity,mImgpIntensity->GetLargestPossibleRegion());	//int image iter

	if( mMtrxNeuron.rows() < 12 || mMtrxGlia.rows() < 12)
	{
		mForeScore = -100;
		mBackScore = 100;
		std::cout << "Total Score" << std::endl << -(mForeScore - mBackScore) << std::endl;
		out << "Fore Score " << mForeScore << " Back Score " << mBackScore << std::endl;
		out << "Total Score " << mForeScore - mBackScore << std::endl;

		oscore << param << ' ' << mForeScore << ' ' << mBackScore << std::endl;
		return;
	}

	vnl_vector<int> back_hist(256,0);	//background histogram
	vnl_vector<int> fore_hist(256,0);	//foreground histogram
	vnl_vector<double> fore_score(256,0.0);	//foreground score for histogram

	for(itI.GoToBegin(),itL.GoToBegin();!itI.IsAtEnd();++itI,++itL)
	{
		//out << (int) itL.Get() << ' ' << (int) itI.Get() << std::endl;
		if( itL.Get() == -1 )	//background pix
		{
			back_hist[(int) itI.Get()]++;	//update back hist
		}
		else
		{
			fore_hist[(int) itI.Get()]++;	//update fore hist
			fore_score[ (int) itI.Get() ] += mvCell[ (int) itL.Get() - 1 ]->mScore;
		}
	}
	out << "Fore Hist" << std::endl;
	for(unsigned int i=0;i<fore_hist.size();i++)	//init hist
	{
		out << std::setw(10) << std::setfill('0') << fore_hist[i] << ' ';
	}
	out << std::endl;

	out << "Fore Score" << std::endl;
	for(unsigned int i=0;i<fore_score.size();i++)	//init hist
	{
		out << std::setw(10) << std::setfill('0') << fore_score[i] << ' ';
	}
	out << std::endl;

	out << "Back Hist" << std::endl;
	for(unsigned int i=0;i<back_hist.size();i++)	//init hist
	{
		out << std::setw(10) << std::setfill('0') << back_hist[i] << ' ';
	}
	out << std::endl;


	double back_score = 0.0;
	for(unsigned int i=0;i<fore_hist.size();i++)	//init hist
	{
		if( fore_hist[i] != 0 )
			back_score += ( fore_score[i] / (double)fore_hist[i] ) * ( (double)fore_hist[i] / (double) fore_hist.sum() ) * (double) back_hist[i] * ( (double)fore_hist[i] / (double)back_hist[i] ) ;

//			back_score += ( fore_score[i] / (double)fore_hist[i] ) * ( exp( (double)fore_hist[i] / (double) fore_hist.sum() ) / exp( (double) back_hist[i] / (double) back_hist.sum() ) ) * (double) back_hist[i];
	}
	mBackScore = back_score * ( (double) back_hist.sum() / (double) fore_hist.sum() );
	//mBackScore = back_score;
	std::cout << "Back" << std::endl;
	std::cout << back_score << std::endl;

	std::cout << "Total Score" << std::endl << -(mForeScore - mBackScore) << std::endl;
	out << "Fore Score " << mForeScore << " Back Score " << mBackScore << std::endl;
	out << "Total Score " << mForeScore - mBackScore << std::endl;

	oscore << param << ' ' << mForeScore << ' ' << mBackScore << ' ' << mForeScore - mBackScore << std::endl;
}


/**	@brief scores foreground
*/
void
cell_mgr::score_foreground(void)
{
	if( mMtrxNeuron.rows() < 12 || mMtrxGlia.rows() < 12)
	{
		mForeScore = -5;
		return;
	}
	vnl_matrix<double> x(mMtrxNeuron.cols(),1);	//temp vector to hold resul

	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] == NULL )
			continue;

		update_volume( mvCell[i] );	//load features
		x(0,0) = mvCell[i]->mVolume;

		update_convexity( mvCell[i] );
		x(1,0) = mvCell[i]->mConvexity;	

		update_shape_fact( mvCell[i] );
		x(2,0) = mvCell[i]->mShapeFact;

		update_eccentricity( mvCell[i] ,bla );
		x(3,0) = mvCell[i]->mEccentricity;

		update_bend_eng( mvCell[i] );
		x(4,0) = mvCell[i]->mBendEng;

		update_vol_grad( mvCell[i] );
		x(5,0) = mvCell[i]->mVolGrad;

		update_bound_grad( mvCell[i] );
		x(6,0) = mvCell[i]->mBoundGrad;

		update_texture( mvCell[i] );
		x(7,0) = mvCell[i]->mTexture;		

		update_avg_int( mvCell[i] );
		x(8,0) = mvCell[i]->mAvgInt ;

		x(9,0) = mvCell[i]->mRadVar;

		update_ints_ratio( mvCell[i] );
		x(10,0) = mvCell[i]->mBoundIntsRatio;

		update_per_nbr( mvCell[i] , bla );
		x(11,0) = mvCell[i]->mPerNbr;

		mvCell[i]->mScore = 0.0;

		double h = 4.0/log(double(mMtrxNeuron.rows()));
		vnl_matrix<double> t(mMtrxNeuron.cols(),1);
		if( mvCell[i]->mClass == 1 )
		{
			for(unsigned  int j=0; j<mMtrxNeuron.rows(); j++ )	//and calc score relative to other cells
			{
				t = x - mMtrxNeuron.get_n_rows(j,1).transpose();
				mvCell[i]->mScore += exp( - ( ( t.transpose()*mMtrxInvCovNeuron*t ).get(0,0) / (2.0*h*h)  ) );
			}
			mvCell[i]->mScore /= double(mMtrxNeuron.rows()) * sqrt(2.0 * M_PI * h * h * mDetNeuron);
		}
		else if( mvCell[i]->mClass == 2 )
		{
			h = 4.0/log(double(mMtrxGlia.rows()));
			for(unsigned  int j=0; j<mMtrxGlia.rows(); j++ )	//and calc score relative to other cells
			{
				t = x - mMtrxGlia.get_n_rows(j,1).transpose();
				mvCell[i]->mScore += exp( - ( ( t.transpose()*mMtrxInvCovGlia*t ).get(0,0) / (2.0*h*h)  ) );
			}
			mvCell[i]->mScore /= double(mMtrxGlia.rows()) * sqrt(2.0 * M_PI * h * h * mDetGlia);
		}
	}
	mForeScore = 0.0;
	for(unsigned  int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] == NULL )
			continue;
		std::cout << mForeScore << std::endl;
		mForeScore += mvCell[i]->mScore * mvCell[i]->mvPix.size();
	}

	std::cout << "Score" << std::endl;
	std::cout << mForeScore <<  std::endl;
}

/**	@brief	get the good train cells from image based on convexity bending energy and shape factor, cell must meet each thresh to be a train cell
*	@param	conv	convexity thresh
*	@param	bend 	bending eng thresh
*	@param	shape	shape factor thresh
*	@param	pnbr	percent pix with a nbr
*/

void
cell_mgr::calc_train_data( const double &conv , const double &bend , const double &shape , const double &pnbr , const double &min_vol,const double &max_vol)
{
	mMtrxTrain.set_size( mvCell.size() , 5 );
	//now generate the mtrx
	for(unsigned  int i=0 ; i < mvCell.size(); i++ )
	{
		update_convexity( mvCell[i] );
		mMtrxTrain( i , 0 ) = mvCell[i]->mConvexity ;

		update_shape_fact( mvCell[i] );
		mMtrxTrain( i , 1 ) = mvCell[i]->mShapeFact ;

		update_bend_eng( mvCell[i] );
		mMtrxTrain( i , 2 ) = mvCell[i]->mBendEng ;

		mMtrxTrain( i , 3 ) = mvCell[i]->changebound ;

		//mMtrxTrain( i , 3 ) = mvCell[i]->bendengstd ;

		//mMtrxTrain( i , 3 ) = mvCell[i]->bendengoutliers ;

		update_per_nbr( mvCell[i] , bla );
		mMtrxTrain( i , 4 ) = mvCell[i]->mPerNbr;
	}
	mMtrxCovTrain = vnl_covariance_mtrx(mMtrxTrain);	

	//calc determinant
	mDetTrain = vnl_determinant(mMtrxCovTrain);

	//calc inv
	mMtrxInvCovTrain = vnl_matrix_inverse<double> (mMtrxCovTrain);

}

/**	@brief	classify neurons and glia
*/
void
cell_mgr::sep_glia_neur(const double &param_text , const double &param_int)
{
	double mean_int = 0.0;
	double mean_text = 0.0;
	double count = 0.0;
	//std::cout << "Getting mean\n";
	for(unsigned  int i=0; i<mvCell.size(); i++ )	//calc mean
	{
		if( mvCell[i] == NULL)
			continue;
		update_avg_int( mvCell[i] );
		update_texture( mvCell[i] );
		mean_int += mvCell[i]->mAvgInt ;
		mean_text += mvCell[i]->mTexture;
		++count;
	}
	mean_int /= count;
	mean_text /= count;
	//std::cout << mean_int << ' ' << mean_text << std::endl;

	double std_int = 0.0;
	double std_text = 0.0;
	//std::cout << "Getting std\n";
	for(unsigned int i=0; i<mvCell.size(); i++ )	//calc std
	{
		if( mvCell[i] == NULL)
			continue;
		std_int += ( mvCell[i]->mAvgInt - mean_int ) * ( mvCell[i]->mAvgInt - mean_int );
		std_text += ( mvCell[i]->mTexture - mean_text ) * ( mvCell[i]->mTexture - mean_text );
	}
	std_int /= count;
	std_text /= count;

	std_int = sqrt( std_int );
	std_text = sqrt( std_text );

	//std::cout << std_int << ' ' << std_text << std::endl;

	//std::cout << "Classifying\n";
	for(unsigned  int i=0; i<mvCell.size(); i++ )	//calssify each cell
	{
		if( mvCell[i] == NULL )
			continue;
		if( mvCell[i]->mTexture > ( mean_text + param_text * std_text ) && mvCell[i]->mAvgInt < ( mean_int + param_int * std_int ) )
		{
			mvCell[i]->mClass = 1;
		}
		else
		{
			mvCell[i]->mClass = 2;
		}
	}
}

void
cell_mgr::clear_train_data(void)
{
	mvTrain.clear();
}

/**	@brief	get more train data
*/
void
cell_mgr::calc_more_train_data(const double &conv , const double &bend , const double &shape , const double &pnbr , const double &min_vol ,const double &max_vol)
{
	//merge tree for each rag
	std::vector<std::vector<int>* > ttree;
	//temp pointer
	std::vector<int> *pt = NULL;
	//front of queue
	std::vector<int>* qfront = NULL;

	cell *c;
	for(unsigned int i=0;i<mvCell.size();i++)
	{		
		//if cell is already merged move on to next cell
		if(mvCell[i] == NULL)
			continue;
		if( check_shape(mvCell[i],conv ,bend ,shape ,pnbr ,min_vol) )	//cell already a train cell
			continue;
		std::vector<std::vector<int>* > ptree;	//paths for current rag
		std::queue<std::vector<int> *> q;

		//queue
		pt = new std::vector<int>();

		//create init path
		pt->push_back(mvCell[i]->mLabel);

		//add path to queue
		q.push(pt);
		ttree.push_back(pt);

		while(!q.empty() && (*q.front()).size() < 8 )
		{
			//get the cell from top of queue
			qfront = q.front();
			//delete top cell
			q.pop();

			//get the cell at end of path
			c = mvCell[qfront->at( qfront->size() - 1) - 1];

			if( c == NULL )	//if call already merged move on
				continue;

			update_nbrs(c,bla);

			for(unsigned int j=0;j<c->mvNbrs.size();j++)
			{
//				std::ofstream ofile("_extramerge.txt",std::ofstream::app);
//				ofile << "Trying: " <<  c->mvNbrs[j] << std::endl;

				//check to make sure that the path is not going in a loop and check if the new cell has not been merged before
				if( std::find( qfront->begin(),qfront->end(),c->mvNbrs[j]) != qfront->end() || mvCell[c->mvNbrs[j]-1] == NULL)
					continue;

				cell *tcell = create_merged_cells( *qfront ,mvCell,mvCell[c->mvNbrs[j]-1]);
				if( tcell == NULL )
				{ 
					continue;
				}
				update_avg_int( tcell );
				update_avg_int( mvCell[c->mvNbrs[j]-1] );
				update_texture( tcell );
				update_texture( mvCell[c->mvNbrs[j]-1] );
				if( tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() > 8000 || sqrt( (tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt)*(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) + (tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture)*(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) + ( 20.0*exp(-5.0 *get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1])) - 1 ) ) > 45.0 )	//check if int and text are fairly similar for new segmented to be added to be merged
				{
	//				std::ofstream ofile("_extramerge.txt",std::ofstream::app);
	//				ofile << "Cell Text and Int not Good Enough: " <<  mvCell[c->mvNbrs[j]-1]->mLabel << ' ' << abs(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) << ' ' << abs(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) << ' ' << tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() << ' ' << get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) << ' ' << sqrt( (tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt)*(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) + (tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture)*(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) + ( exp(1-get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1])) - 1 )*2000.0 ) << std::endl;
//					ofile.close();
					if( (*qfront).size() != 1 )	//if its a merged cell delete it
						delete tcell;
					continue;
				}
//				ofile << "Cell Text and Int Good:  " <<  mvCell[c->mvNbrs[j]-1]->mLabel << ' ' << abs(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) << ' ' << abs(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) << ' ' << tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() << ' ' << get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) << ' ' << sqrt( (tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt)*(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) + (tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture)*(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) + ( 0.5 - get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1]) )*400.0 ) << std::endl;


				//add the nbr to the path
				pt = new std::vector<int>(*qfront);
				pt->push_back(c->mvNbrs[j]);
				//check if the path is already in tree
				if(!compare_path_to_tree(ptree,pt))
				{
//					ofile << "Unique path" << std::endl;
					tcell = create_merged_cells( *pt , mvCell );
					if( check_shape_dist(tcell,conv ,bend ,shape ,pnbr ,min_vol) )
					{
						update_cells(*pt,tcell);
						delete pt;
						mvTrain.push_back(tcell);
//						std::ofstream ofile("_extramerge.txt",std::ofstream::app);
//						ofile << "Found new train cell!!" << std::endl << std::endl;
//						ofile.close();
						if( mvTrain.size() > 20 )
							goto jump2;
					}
					else
					{
//						ofile << "features not met" << std::endl;
						delete tcell;
						q.push(pt);
						ttree.push_back(pt);
					}
				}
//				ofile << "Moving on" << std::endl;
//				ofile.close();
			}
		}
	}
jump2:
	for(unsigned int i=0;i<ttree.size();i++)
		delete ttree[i];
}

/**	@brief	checks shape of cell
*	@param	conv convex hull
*	@param	bend bending energy
*	@param	shape shape factor
*	@param	pnbr percent neighboor
*	@param	min_vol minimum volume
*	@return true if cell is good, false cells fails at least one of features
*/

bool
cell_mgr::check_shape( cell* c , const double &conv , const double &bend , const double &shape , const double &pnbr , const double &vol )
{
	if( c == NULL )
		return false;

	//check each feature
	update_per_nbr( c , bla );
	if( c->mPerNbr >= pnbr )
		return false;

	update_volume(  c );
	if( c->mVolume <= vol )
		return false;

	update_convexity( c );
	if( c->mConvexity <= conv )
		return false;

	update_shape_fact( c );
	if(  c->mShapeFact >= shape )
		return false;

	update_bend_eng( c );
	if( c->mBendEng >= bend )
		return false;

	return true;
}

/**	@brief	classifes the cells using k-means clustering
*/
void
cell_mgr::classify_cells(void)
{
	vnl_matrix<double> mean(2,4,0.0);	//0-neuron,1-glia
	int iter = 0;	//iterations of k clustering
	double count1 = 0.0 , count2 = 0.0;		//keep track of how many cells classified into each group

	for(unsigned int i=0; i<mvCell.size(); i++ )	//find initial means, max and min
	{
		if( mvCell[i] == NULL )
			continue;
		update_texture( mvCell[i] );
		update_avg_int( mvCell[i] );
		update_bound_grad( mvCell[i] );
		update_vol_grad( mvCell[i] );

		mean(0,0) += mvCell[i]->mAvgInt;
		mean(0,1) += mvCell[i]->mTexture;
		mean(0,2) += mvCell[i]->mBoundGrad;
		mean(0,3) += mvCell[i]->mVolGrad;

		mean(1,0) += mvCell[i]->mAvgInt;
		mean(1,1) += mvCell[i]->mTexture;
		mean(1,2) += mvCell[i]->mBoundGrad;
		mean(1,3) += mvCell[i]->mVolGrad;

		count1++;
	}	

	mean /= count1;

	mean(0,0) -= 10;
	mean(0,1) += 10;
	mean(0,2) -= 10;
	mean(0,3) -= 10;

	mean(1,0) += 10;
	mean(1,1) -= 10;
	mean(1,2) += 10;
	mean(1,3) += 10;


	vnl_matrix<double> oldmean(2,4,0.0);
	//	std::cout << "Clustering\n";

	count1 = count2 = 1.0;	

	while( iter < 40 && (mean - oldmean).absolute_value_sum() != 0.0  )
	{
		std::cout << iter << std::endl;

		std::cout << mean << std::endl;
		oldmean = mean;
		mean.fill(0.0);
		count1 = count2 = 0.0;

		for(unsigned  int i=0; i<mvCell.size(); i++ )
		{
			//			std::cout << i << std::endl;
			if( mvCell[i] == NULL )
				continue;

			vnl_vector<double> t(4,0.0);
			t(0) = mvCell[i]->mAvgInt;
			t(1) = mvCell[i]->mTexture;
			t(2) = mvCell[i]->mBoundGrad;
			t(3) = mvCell[i]->mVolGrad;

			//std::cout << ( t - oldmean.get_row(0) ).magnitude() << ' ' << ( t - oldmean.get_row(1) ).magnitude() << std::endl;

			//			std::cout << "Comparing\n";
			if( ( t - oldmean.get_row(0) ).magnitude() > ( t - oldmean.get_row(1) ).magnitude() )	//mean2 is closer so classify cell as 2
			{
				//				std::cout << "Glia cell\n";
				mvCell[i]->mClass = 2;
				mean(1,0) += mvCell[i]->mAvgInt;
				mean(1,1) += mvCell[i]->mTexture;
				mean(1,2) += mvCell[i]->mBoundGrad;
				mean(1,3) += mvCell[i]->mVolGrad;
				count2++;
			}
			else
			{
				//				std::cout << "Neuron cell\n";
				mvCell[i]->mClass = 1;
				mean(0,0) += mvCell[i]->mAvgInt;
				mean(0,1) += mvCell[i]->mTexture;
				mean(0,2) += mvCell[i]->mBoundGrad;
				mean(0,3) += mvCell[i]->mVolGrad;
				count1++;
			}
		}
		iter++;
		mean(0,0) /= count1;
		mean(0,1) /= count1;
		mean(0,2) /= count1;
		mean(0,3) /= count1;

		mean(1,0) /= count2;
		mean(1,1) /= count2;
		mean(1,2) /= count2;
		mean(1,3) /= count2;
	}
	std::cout << mean << std::endl;
}

/**	@brief	merges the cells, builds rag scores cells 
*/
void
cell_mgr::build_rag(const int &mode,std::string &file,const int &rDepth)
{
	//temp pointer
	std::vector<int> *pt = NULL;
	//front of queue
	std::vector<int>* qfront = NULL;
	std::ofstream clearfile("_mergecells.txt");
	clearfile.close();
	std::vector<std::vector<int>* > ttree;

	cell *c;
	for(unsigned int i=0;i<mvCell.size();i++)
	{		
		std::vector<std::vector<int>* > ptree;	//paths for current rag
		cell * pMxScC = NULL;	//ptr to cell with max score
		std::queue<std::vector<int> *> q;

		std::cout << "On Cell: " << i+1 << " of " << mvCell.size() << std::endl;

		//if cell is already merged move on to next cell
		if(mvCell[i] == NULL)
		{
			std::ofstream ofile("_mergecells.txt",std::ofstream::app);
			ofile << "skipping: " << i+1 << std::endl << std::endl;
			ofile.close();
			continue;
		}
		//queue
		pt = new std::vector<int>();

		//create init path
		pt->push_back(mvCell[i]->mLabel);

		//add path to queue
		q.push(pt);
		ptree.push_back(pt);
		ttree.push_back(pt);

		while(!q.empty() && (*q.front()).size() < rDepth)
		{
			//get the cell from top of queue
			qfront = q.front();
			//delete top cell
			q.pop();

			//get the cell at end of path
			c = mvCell[qfront->at( qfront->size() - 1) - 1];

			if( c == NULL )	//if call already merged move on
				continue;

			update_nbrs(c,bla);

			if(c->mvNbrs.empty())	//if a loner cell score it and move on
			{
				std::ofstream ofile("_mergecells.txt",std::ofstream::app);
				ofile << "loner: " << i+1 << std::endl << std::endl;
				ofile.close();
				score_cell( mvCell[i] ,bla );
				continue;
			}
			//check each nbr
			for(unsigned int j=0;j<c->mvNbrs.size();j++)
			{
				//check to make sure that the path is not going in a loop and check if the new cell has not been merged before
				if( std::find( qfront->begin(),qfront->end(),c->mvNbrs[j]) != qfront->end() || mvCell[c->mvNbrs[j]-1] == NULL)
				{
					continue;
				}	

				//add the nbr to the path
				pt = new std::vector<int>(*qfront);
				pt->push_back(c->mvNbrs[j]);
				//check if the path is already in tree
				if(!compare_path_to_tree(ptree,pt))
				{
					cell *tcell = create_merged_cells( *pt ,mvCell,mvCell[c->mvNbrs[j]-1]);
					if( tcell != NULL )
					{
						if( tcell->mvPix.size() < 8500 )
						{
							//add new path
							ptree.push_back(pt);
						}
					}
					delete tcell;
				}
				q.push(pt);
				ttree.push_back(pt);
			}
		}
		//std::cout << "Built tree\n";
		if(ptree.size() < 2)	//if there is nothing to merge move on
		{
			std::ofstream ofile("_mergecells.txt",std::ofstream::app);
			ofile << "nothing to merge: " << i+1 << ' ' << mvCell[i]->mScore << std::endl << std::endl;
			ofile.close();
			score_cell( mvCell[i] ,bla );
			continue;
		}

		double max_score = 0.0 ;	//max score of merged cells, init to current cell score

		int max_miIndex = -1;
		std::ofstream out("_mergecells.txt",std::ofstream::app);
		score_cell( mvCell[ i ] , bla );
		c = mvCell[ i ];
		out << "trying: " << (*ptree[0]).at(0) << " F: " << c->mVolume << ' ' << c->mConvexity << ' ' <<  c->mShapeFact << ' ' << c->mEccentricity << ' ' << c->mBendEng << ' ' << c->mVolGrad << ' ' << c->mBoundGrad << ' ' << c->mTexture << ' ' << c->mAvgInt << ' ' << c->mRadVar << ' ' << c->mBoundIntsRatio << ' ' << c->mPerNbr << ' ' << max_score << std::endl;

		c = NULL;
		//create new cells with each new path
		for(unsigned int j=1;j<ptree.size();j++)
		{

			out << "trying: ";
			for(unsigned int k=0;k<(*ptree[j]).size();k++)
			{
				out << (*ptree[j]).at(k) << ' ';
			}

			c = create_merged_cells( *ptree[j] , mvCell );

			//score each cell
			score_cell(c,*ptree[j]);


			out << "F: " << c->mVolume << ' ' << c->mConvexity << ' ' <<  c->mShapeFact << ' ' << c->mEccentricity << ' ' << c->mBendEng << ' ' << c->mVolGrad << ' ' << c->mBoundGrad << ' ' << c->mTexture << ' ' << c->mAvgInt << ' ' << c->mRadVar << ' ' << c->mBoundIntsRatio << ' ' << c->mPerNbr;
			out << " S: " << c->mScore;
			//check if new cell has a higher score than n
			if(c->mScore > max_score)
			{
				for(unsigned int k=0;k<(*ptree[j]).size();k++)
				{
					score_cell( mvCell[ (*ptree[j]).at(k) - 1 ],bla );
					if( mvCell[ (*ptree[j]).at(k) - 1 ]->mScore > c->mScore )
					{
						out << std::endl << "Merging Object Still Better: " << (*ptree[j]).at(k) << ' ' << mvCell[ (*ptree[j]).at(k) - 1 ]->mScore << std::endl;
						delete c;
						c = NULL;
						break;
					}
				}
				if( c == NULL )
					continue;
				max_score = c->mScore;	//if yes set new max score
				max_miIndex = j;	//set max index
				delete pMxScC;	//delete prev max score cell
				pMxScC = c;	//set ptr to new max cell
				c = NULL;
			}
			else
			{
				delete c;
				c = NULL;
			}
			//out << std::endl;
			out << ' ' << max_score << ' ' << max_miIndex << std::endl;

		}

		//if( max_miIndex != -1 )
		//{
		//	for(int k=0;k<(*ptree[max_miIndex]).size();k++)
		//	{
		//		score_cell( mvCell[ (*ptree[max_miIndex]).at(k) - 1 ],bla );
		//		cell *t = mvCell[ (*ptree[max_miIndex]).at(k) - 1 ];
		//		if( mvCell[ (*ptree[max_miIndex]).at(k) - 1 ]->mScore > max_score )
		//		{
		//			out << "Merging Object Still Better: " << (*ptree[max_miIndex]).at(k) << ' ' << "F: " << t->mVolume << ' ' << t->mConvexity << ' ' <<  t->mShapeFact << ' ' << t->mEccentricity << ' ' << t->mBendEng << ' ' << t->mVolGrad << ' ' << t->mBoundGrad << ' ' << t->mTexture << ' ' << t->mAvgInt << ' ' << t->mRadVar << ' ' << t->mBoundIntsRatio << ' ' << t->mPerNbr << ' ' << mvCell[ (*ptree[max_miIndex]).at(k) - 1 ]->mScore << std::endl;
		//			max_miIndex = -1;
		//			delete pMxScC;
		//			pMxScC = NULL;
		//			break;
		//		}
		//	}
		//}

		//if there is a cell with higher score than current update cells
		if( max_miIndex != -1 )
		{

			out << "merging: ";
			for(unsigned int k=0;k<(*ptree[max_miIndex]).size();k++)
			{
				out << (*ptree[max_miIndex]).at(k) << ' ';
			}
			update_cells(*ptree[max_miIndex],pMxScC);	//del rest of merged cells
			pMxScC = NULL;
			mvCell[i]->mIsMerged = true;
			out << std::endl;
			out << mvCell[i]->mLabel << ' ' << mvCell[i]->mScore << std::endl;
			--i;	//still need to build rag for new merged cell
		}
		out << std::endl;
		out.close();
	}
	for(unsigned int j=0; j<ttree.size(); j++)
		delete ttree[j];
}

/**	@brief	extracts the training matrices from cells, calc cov and inv cov and det
*	@param	calcMtrx calculate covariance matrices
*/
void
cell_mgr::calc_ng_data(const bool &calcMtrx)
{
	int mNcount = 0 , mGcount = 0, count = 0;
	for(unsigned int i=0 ; i < mvCell.size(); i++ )	//count actual amount of neuron cells, some are NULL due to merging
	{
		if( mvCell[i] != NULL)
		{
			if( mvCell[i]->mClass == 1 ) //count neurons
				mNcount++;
			else if( mvCell[i]->mClass == 2 )
				mGcount++;
			count++;
		}
	}
	mMtrxNeuron.clear();
	mMtrxGlia.clear();
	mNoclass.clear();
	mMtrxNeuron.set_size( mNcount,14 );	//matrices for scoring
	mMtrxGlia.set_size( mGcount,14 );
	mNoclass.set_size(count,14);

	mNcount = mGcount = count = 0;
	for(unsigned int i=0 ; i < mvCell.size(); i++ )
	{
		if( mvCell[i] == NULL)
			continue;
		update_volume(mvCell[i]);	//load features
		mNoclass(count,0) = mvCell[i]->mVolume;

		update_convexity(mvCell[i]);
		mNoclass(count,1) = mvCell[i]->mConvexity;	

		update_shape_fact(mvCell[i]);
		mNoclass(count,2) = mvCell[i]->mShapeFact;

		update_eccentricity(mvCell[i],bla);
		mNoclass(count,3) = mvCell[i]->mEccentricity;

		update_bend_eng(mvCell[i]);
		mNoclass(count,4) = mvCell[i]->mBendEng;

		update_vol_grad(mvCell[i]);
		mNoclass(count,5) = mvCell[i]->mVolGrad;

		update_bound_grad(mvCell[i]);
		mNoclass(count,6) = mvCell[i]->mBoundGrad;

		update_texture( mvCell[i] );
		mNoclass(count,7) = mvCell[i]->mTexture;		

		update_avg_int( mvCell[i] );
		mNoclass(count,8) = mvCell[i]->mAvgInt ;

		mNoclass(count,9) = mvCell[i]->mRadVar;

		update_ints_ratio(  mvCell[i] );
		mNoclass(count,10) = mvCell[i]->mBoundIntsRatio;

		update_per_nbr(  mvCell[i] , bla );
		mNoclass(count,11) = mvCell[i]->mPerNbr;

		update_depth(  mvCell[i] );
		mNoclass(count,12) = mvCell[i]->mDepth;

		mNoclass(count,13) = (double) mvCell[i]->mvBoundPix2D.size();

		count++;

		mvCell[i]->mScore = 0.0;	//reset score

		if( mvCell[i]->mClass == 1 )
		{
			mMtrxNeuron(mNcount,0) = mvCell[i]->mVolume;

			mMtrxNeuron(mNcount,1) = mvCell[i]->mConvexity;	
			
			mMtrxNeuron(mNcount,2) = mvCell[i]->mShapeFact;

			mMtrxNeuron(mNcount,3) = mvCell[i]->mEccentricity;
			
			mMtrxNeuron(mNcount,4) = mvCell[i]->mBendEng;
			
			mMtrxNeuron(mNcount,5) = mvCell[i]->mVolGrad;

			mMtrxNeuron(mNcount,6) = mvCell[i]->mBoundGrad;

			mMtrxNeuron(mNcount,7) = mvCell[i]->mTexture;		

			mMtrxNeuron(mNcount,8) = mvCell[i]->mAvgInt ;

			mMtrxNeuron(mNcount,9) = mvCell[i]->mRadVar;

			mMtrxNeuron(mNcount,10) = mvCell[i]->mBoundIntsRatio;

			mMtrxNeuron(mNcount,11) = mvCell[i]->mPerNbr;

			mMtrxNeuron(mNcount,12) = mvCell[i]->mDepth;

			mMtrxNeuron(mNcount,13) = (double) mvCell[i]->mvBoundPix2D.size();

			mNcount++;

			mvCell[i]->mScore = 0.0;	//reset score
			//mvCell[i]->mIsTrain = true;
		}
		else if( mvCell[i]->mClass == 2 )
		{
			mMtrxGlia(mGcount,0) = mvCell[i]->mVolume;

			mMtrxGlia(mGcount,1) = mvCell[i]->mConvexity;	

			mMtrxGlia(mGcount,2) = mvCell[i]->mShapeFact;

			mMtrxGlia(mGcount,3) = mvCell[i]->mEccentricity;

			mMtrxGlia(mGcount,4) = mvCell[i]->mBendEng;

			mMtrxGlia(mGcount,5) = mvCell[i]->mVolGrad;

			mMtrxGlia(mGcount,6) = mvCell[i]->mBoundGrad;

			mMtrxGlia(mGcount,7) = mvCell[i]->mTexture;		

			mMtrxGlia(mGcount,8) = mvCell[i]->mAvgInt ;

			mMtrxGlia(mGcount,9) = mvCell[i]->mRadVar;

			mMtrxGlia(mGcount,10) = mvCell[i]->mBoundIntsRatio;

			mMtrxGlia(mGcount,11) = mvCell[i]->mPerNbr;
			
			mMtrxGlia(mGcount,12) = mvCell[i]->mDepth;

			mMtrxGlia(mGcount,13) = (double) mvCell[i]->mvBoundPix2D.size();
			
			mGcount++;

			mvCell[i]->mScore = 0.0;	//reset score
			//mvCell[i]->mIsTrain = true;
		}

	}
	if( !calcMtrx )
		return;
	mMtrxCovNeuron.clear();		//calc matrices
	mMtrxCovGlia.clear();
	mMtrxCovAll.clear();
		
	//calc all cov, inv and det of the matrices
	mMtrxCovNeuron = vnl_covariance_mtrx(mMtrxNeuron);
	mMtrxCovGlia = vnl_covariance_mtrx(mMtrxGlia);
	mMtrxCovAll = vnl_covariance_mtrx(mNoclass);

	//calc determinant
	mDetNeuron = vnl_determinant(mMtrxCovNeuron);
	mDetGlia = vnl_determinant(mMtrxCovGlia);
	mDetAll = vnl_determinant(mMtrxCovAll);

	mMtrxInvCovNeuron.clear();
	mMtrxInvCovGlia.clear();
	mMtrxInvCovAll.clear();

	mMtrxInvCovNeuron = vnl_matrix_inverse<double> ( mMtrxCovNeuron );
	mMtrxInvCovGlia = vnl_matrix_inverse<double> ( mMtrxCovGlia );
	mMtrxInvCovAll = vnl_matrix_inverse<double> ( mMtrxCovAll );

	if( mDetNeuron < 1.0 )
	{
		mMtrxInvCovNeuron = mMtrxInvCovAll;
		mDetNeuron = mDetAll;
	}
	if( mDetGlia < 1.0 )
	{
		mMtrxInvCovGlia = mMtrxInvCovAll;
		mDetGlia = mDetAll;
	}

}

/**	@brief	score cell based on features
*	@param	c	cell to score
*	@param	v	vector of marged cells cell c is made of, empty if cell not merged
*/

void
cell_mgr::score_cell(cell *c,std::vector<int> v)
{
	vnl_matrix<double> x(mMtrxNeuron.cols(),1);	//temp vector to hold resul

	update_volume( c );	//load features
	x(0,0) = c->mVolume;

	update_convexity( c );
	x(1,0) = c->mConvexity;	

	update_shape_fact( c );
	x(2,0) = c->mShapeFact;

	update_eccentricity( c ,v);
	x(3,0) = c->mEccentricity;

	update_bend_eng( c );
	x(4,0) = c->mBendEng;

	update_vol_grad( c );
	x(5,0) = c->mVolGrad;

	update_bound_grad( c );
	x(6,0) = c->mBoundGrad;

	update_texture( c );
	x(7,0) = c->mTexture;		

	update_avg_int( c );
	x(8,0) = c->mAvgInt ;

	x(9,0) = c->mRadVar;

	update_ints_ratio( c );
	x(10,0) = c->mBoundIntsRatio;

	update_per_nbr( c , v );
	x(11,0) = c->mPerNbr;

	update_depth( c );
	x(12,0) = c->mDepth;

	x(13,0) = (double) c->mvBoundPix2D.size();


	double h = 4.0/log(double(mMtrxNeuron.rows()));
	vnl_matrix<double> t(mMtrxNeuron.cols(),1);
	for(unsigned int i=0; i<mMtrxNeuron.rows(); i++ )	//and calc score relative to other cells
	{
		t = x - mMtrxNeuron.get_n_rows(i,1).transpose();
		c->mScore += exp( - ( ( t.transpose()*mMtrxInvCovNeuron*t ).get(0,0) / (2.0*h*h)  ) );
	}
	c->mScore /= double(mMtrxNeuron.rows()) * sqrt(2.0 * M_PI * h * h * mDetNeuron);
	c->mClass = 1;

	h = 4.0/log(double(mMtrxGlia.rows()));
	double temp_score = 0.0;
	for(unsigned int i=0; i<mMtrxGlia.rows(); i++ )	//and calc score relative to other cells
	{
		t = x - mMtrxGlia.get_n_rows(i,1).transpose();
		temp_score += exp( - ( ( t.transpose()*mMtrxInvCovGlia*t ).get(0,0) / (2.0*h*h)  ) );
	}
	temp_score /= double(mMtrxGlia.rows()) * sqrt(2.0 * M_PI * h * h * mDetGlia);

	if( temp_score > c->mScore)
	{
		c->mScore = temp_score;
		c->mClass = 2;
	}
}

/**	@brief	relabel cells and image and erase NULL cells
*/
void
cell_mgr::relabel_orig_cells(void)
{
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	//relabel the cells and the image
	for(unsigned int i=0;i<mvOrigCell.size();i++)
	{
		if( mvCell[i] != mvOrigCell[i] )
		{
			if( mvCell[i] != NULL )
				if( mvCell[i]->mIsMerged == true )
					delete mvCell[i];
			mvCell[i] = mvOrigCell[i];
		}
		for(unsigned int j=0;j<mvOrigCell[i]->mvPix.size();j++)
		{
			miIndex[0] = mvOrigCell[i]->mvPix[j]->x_;
			miIndex[1] = mvOrigCell[i]->mvPix[j]->y_;
			miIndex[2] = mvOrigCell[i]->mvPix[j]->z_;
			it.SetIndex(miIndex);
			it.Set(mvOrigCell[i]->mLabel);
		}
		update_center(mvOrigCell[i]);
	}
}

/**	@brief	gives the smaller cell percent of pix nbring with bigger cell
*	@param	c1	cell ptr to first cell
*	@param	c2	cell ptr to second cell
*	@return	percent of pix between cells for smaller cell
*/
double
cell_mgr::get_perc_shared( cell *c1 , cell *c2 )
{
	if( c1 == NULL || c2 == NULL || c1 == c2 )
		return 0.0;
	cell *cs , *cb;	//small cell , big cell
	if( c1->mvPix.size() < c2->mvPix.size() )	//see which cell is smaller
	{
		cs = c1;
		cb = c2;
	}
	else
	{
		cs = c2;
		cb = c1;
	}
	std::vector<int>::iterator it = std::find( cs->mvNbrs.begin() , cs->mvNbrs.end() , cb->mLabel );	//find the cell in nbrs
	return (double) cs->mvNbrsPix[ it - cs->mvNbrs.begin() ] / (double) cs->mvBoundPix2D.size();	//return the perc
	return 0.0;
}

/**	@brief	calculates avg mah dist of each class. Then based on param treats cells within
*	dist as training cells
*	@param	param	paramter for cutoff (avg_mah_dist*param)
*/
void
cell_mgr::calc_mean_vect(const double &param)
{
	mMeanN.set_size(mMtrxNeuron.cols(),1);
	mMeanG.set_size(mMtrxGlia.cols(),1);

	for(unsigned int i=0; i<mMtrxNeuron.cols(); i++)		//calc mean vectors for features
	{
		mMeanN(i,0) = ( mMtrxNeuron.get_column(i) ).mean();
	}
	for(unsigned int i=0; i<mMtrxGlia.cols(); i++)
	{
		mMeanG(i,0) = ( mMtrxGlia.get_column(i) ).mean();
	}
	double rN = 0.0, rG = 0.0;	//avg mahobis distance, neuron
	for(unsigned int i=0; i<mMtrxNeuron.rows(); i++)
	{
		rN += sqrt( ( ( mMtrxNeuron.get_n_rows(i,1) - mMeanN ) * mMtrxInvCovNeuron * ( mMtrxNeuron.get_n_rows(i,1) - mMeanN ).transpose() ).get(0,0) );
	}
	rN /= (double) mMtrxNeuron.rows();
	rN *= param;

	for(unsigned int i=0; i<mMtrxGlia.rows(); i++)	//avg mahobis distance, glia
	{
		rG += sqrt( ( ( mMtrxGlia.get_n_rows(i,1) - mMeanG ) * mMtrxInvCovGlia * ( mMtrxGlia.get_n_rows(i,1) - mMeanG ).transpose() ).get(0,0) );
	}
	rG /= (double) mMtrxGlia.rows();
	rG *= param;
	unsigned int count = 0;
	while( count < 30)
	{
		for(unsigned int i=0; i<mMtrxNeuron.rows(); i++)
		{
			if( sqrt( ( ( mMtrxNeuron.get_n_rows(i,1) - mMeanN ) * mMtrxInvCovNeuron * ( mMtrxNeuron.get_n_rows(i,1) - mMeanN ).transpose() ).get(0,0) ) < rN)
				count++;
		}
		if( count >= mMtrxNeuron.rows() )
			break;
		if( count < 30 )
		{
			rN *= 1.1;
			count = 0;
		}
	}
	count = 0;
	while( (count < 30) )
	{
		for(unsigned int i=0; i<mMtrxGlia.rows(); i++)	//avg mahobis distance, glia
		{
			if( sqrt( ( ( mMtrxGlia.get_n_rows(i,1) - mMeanG ) * mMtrxInvCovGlia * ( mMtrxGlia.get_n_rows(i,1) - mMeanG ).transpose() ).get(0,0) ) < rG)
				count++;
		}
		if( count >= mMtrxGlia.rows() )
			break;
		if( count < 30 )
		{
			rG *= 1.1;
			count = 0;
		}
	}

	std::ofstream out("mahdist.txt");
	out << "Avg N: " << rN << std::endl;
	out << "Avg G: " << rG << std::endl;

	vnl_matrix<double> x(mMtrxNeuron.cols(),1);	//temp vector
	for(unsigned int i=0 ; i < mvCell.size(); i++ )
	{
		if( mvCell[i] == NULL)
			continue;
		update_volume( mvCell[i] );	//load features
		x(0,0) = mvCell[i]->mVolume;

		update_convexity( mvCell[i] );
		x(1,0) = mvCell[i]->mConvexity;	

		update_shape_fact( mvCell[i] );
		x(2,0) = mvCell[i]->mShapeFact;

		update_eccentricity( mvCell[i] ,bla);
		x(3,0) = mvCell[i]->mEccentricity;

		update_bend_eng( mvCell[i] );
		x(4,0) = mvCell[i]->mBendEng;

		update_vol_grad( mvCell[i] );
		x(5,0) = mvCell[i]->mVolGrad;

		update_bound_grad( mvCell[i] );
		x(6,0) = mvCell[i]->mBoundGrad;

		update_texture( mvCell[i] );
		x(7,0) = mvCell[i]->mTexture;		

		update_avg_int( mvCell[i] );
		x(8,0) = mvCell[i]->mAvgInt ;

		x(9,0) = mvCell[i]->mRadVar;

		update_ints_ratio( mvCell[i] );
		x(10,0) = mvCell[i]->mBoundIntsRatio;

		update_per_nbr( mvCell[i] , bla );
		x(11,0) = mvCell[i]->mPerNbr;

		update_depth( mvCell[i] );
		x(12,0) = mvCell[i]->mDepth;
		
		x(13,0) = (double) mvCell[i]->mvBoundPix2D.size();

		out << i+1 << ' ' << sqrt( ( ( x - mMeanN ).transpose() * mMtrxInvCovNeuron * ( x - mMeanN ) ).get(0,0) ) << std::endl;

		if( mvCell[i]->mClass == 1 )	//check if mah dist bigger than avg, if yes get rid of cell
		{
			if( sqrt( ( ( x - mMeanN ).transpose() * mMtrxInvCovNeuron * ( x - mMeanN ) ).get(0,0) ) > rN )
			{
				mvCell[i]->mClass = -1;
			}
		}
		else if( mvCell[i]->mClass == 2 )
		{
			if( sqrt( ( ( x - mMeanG ).transpose() * mMtrxInvCovGlia * ( x - mMeanG ) ).get(0,0) ) > rG )
			{
				mvCell[i]->mClass = -1;
			}
		}
	}
}

/**	@brief	find train cells that meet criteria
*	@param	conv convex hull
*	@param	bend bending energy
*	@param	shape shape factor
*	@param	pnbr percent neighboor
*	@param	min_vol minimum volume
*/
void
cell_mgr::check_all_cells(const double &conv , const double &bend , const double &shape , const double &pnbr , const	double &min_vol )
{
	for(unsigned int i=0 ; i < mvCell.size(); i++ )	//go through each cell and find train cells
	{ 
		if( mvCell[i] == NULL )
			continue;
		if( !check_shape(mvCell[i],conv, bend,shape ,pnbr ,min_vol) )
		{
			//cell is a train cell, meets all feature req
			if( mvCell[i]->mIsMerged )
				delete mvCell[i];
			mvCell[i] = NULL;
		}
	}
}

void
cell_mgr::build_rag_mah(void)
{
	//temp pointer
	std::vector<int> *pt = NULL;
	//front of queue
	std::vector<int>* qfront = NULL;
	std::ofstream clearfile("_mergecells.txt");
	clearfile.close();
	std::vector<std::vector<int>* > ttree;
	cell *c;
	for(unsigned int i=0;i<mvCell.size();i++)
	{		
		std::vector<std::vector<int>* > ptree;	//paths for current rag
		cell * pMxScC = NULL;	//ptr to cell with max score
		std::queue<std::vector<int> *> q;

		std::cout << "On Cell: " << i+1 << " of " << mvCell.size() << std::endl;

		//if cell is already merged move on to next cell
		if(mvCell[i] == NULL)
		{
			std::ofstream ofile("_mergecells.txt",std::ofstream::app);
			ofile << "skipping: " << i+1 << std::endl << std::endl;
			ofile.close();
			continue;
		}
		//queue
		pt = new std::vector<int>();

		//create init path
		pt->push_back(mvCell[i]->mLabel);

		//add path to queue
		q.push(pt);
		ptree.push_back(pt);
		ttree.push_back(pt);

		while(!q.empty() && (*q.front()).size() < 5 )
		{
			//get the cell from top of queue
			qfront = q.front();
			//delete top cell
			q.pop();

			//get the cell at end of path
			c = mvCell[qfront->at( qfront->size() - 1) - 1];

			if( c == NULL )	//if call already merged move on
				continue;

			update_nbrs(c,bla);

			if(c->mvNbrs.empty())	//if a loner cell score it and move on
			{
				std::ofstream ofile("_mergecells.txt",std::ofstream::app);
				ofile << "loner: " << i+1 << std::endl << std::endl;
				ofile.close();
				score_cell_mah( mvCell[i] ,bla );
				continue;
			}
			//check each nbr
			for(unsigned int j=0;j<c->mvNbrs.size();j++)
			{
				//check to make sure that the path is not going in a loop and check if the new cell has not been merged before
				if( std::find( qfront->begin(),qfront->end(),c->mvNbrs[j]) != qfront->end() || mvCell[c->mvNbrs[j]-1] == NULL)
				{
					continue;
				}	

				cell *tcell = create_merged_cells( *qfront ,mvCell,mvCell[c->mvNbrs[j]-1]);
				if( tcell == NULL )
				{ 
					continue;
				}
				update_avg_int( tcell );
				update_avg_int( mvCell[c->mvNbrs[j]-1] );
				update_texture( tcell );
				update_texture( mvCell[c->mvNbrs[j]-1] );
				if( tcell->mvPix.size() + mvCell[c->mvNbrs[j]-1]->mvPix.size() > 8000 || sqrt( (tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt)*(tcell->mAvgInt - mvCell[c->mvNbrs[j]-1]->mAvgInt) + (tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture)*(tcell->mTexture - mvCell[c->mvNbrs[j]-1]->mTexture) + ( (1-get_perc_shared(tcell,mvCell[c->mvNbrs[j]-1])) - 1 )*60.0 ) > 150.0 )	//check if int and text are fairly similar for new segmented to be added to be merged
				{
					if( (*qfront).size() != 1 )	//if its a merged cell delete it
						delete tcell;
					continue;
				}

				//add the nbr to the path
				pt = new std::vector<int>(*qfront);
				pt->push_back(c->mvNbrs[j]);
				//check if the path is already in tree
				if(!compare_path_to_tree(ptree,pt))
				{
					//add new path
					ptree.push_back(pt);
				}
				q.push(pt);
				ttree.push_back(pt);
			}
		}
		//std::cout << "Built tree\n";
		if(ptree.size() < 2)	//if there is nothing to merge move on
		{
			std::ofstream ofile("_mergecells.txt",std::ofstream::app);
			ofile << "nothing to merge: " << i+1 << ' ' << mvCell[i]->mScore << std::endl << std::endl;
			ofile.close();
			score_cell_mah( mvCell[i] ,bla );
			continue;
		}
		score_cell_mah( mvCell[ i ] ,bla);	//score first cell
		double max_score = mvCell[ i ]->mScore ;	//max score of merged cells, init to current cell score

		int max_miIndex = -1;
		std::ofstream out("_mergecells.txt",std::ofstream::app);
		c = mvCell[ i ];
		out << "trying: " << (*ptree[0]).at(0) << " F: " << c->mVolume << ' ' << c->mConvexity << ' ' <<  c->mShapeFact << ' ' << c->mEccentricity << ' ' << c->mBendEng << ' ' << c->mVolGrad << ' ' << c->mBoundGrad << ' ' << c->mTexture << ' ' << c->mAvgInt << ' ' << c->mRadVar << ' ' << c->mBoundIntsRatio << ' ' << c->mPerNbr << ' ' << max_score << std::endl;

		c = NULL;
		//create new cells with each new path
		for(unsigned int j=1;j<ptree.size();j++)
		{

			out << "trying: ";
			for(unsigned int k=0;k<(*ptree[j]).size();k++)
			{
				out << (*ptree[j]).at(k) << ' ';
			}

			c = create_merged_cells( *ptree[j] , mvCell );

			//score each cell
			score_cell_mah(c,*ptree[j]);


			out << "F: " << c->mVolume << ' ' << c->mConvexity << ' ' <<  c->mShapeFact << ' ' << c->mEccentricity << ' ' << c->mBendEng << ' ' << c->mVolGrad << ' ' << c->mBoundGrad << ' ' << c->mTexture << ' ' << c->mAvgInt << ' ' << c->mRadVar << ' ' << c->mBoundIntsRatio << ' ' << c->mPerNbr;
			out << ' ' << c->mScore;
			//check if new cell has a higher score than n
			if(c->mScore < max_score)
			{
				max_score = c->mScore;	//if yes set new max score
				max_miIndex = j;	//set max index
				if( pMxScC != NULL)	//delete prev max score cell
					delete pMxScC;
				pMxScC = c;	//set ptr to new max cell
				c = NULL;
			}
			else
			{
				delete c;
				c = NULL;
			}
			out << ' ' << max_score << ' ' << max_miIndex << std::endl;

		}


		//if there is a cell with higher score than current update cells
		if( max_miIndex != -1 )
		{

			out << "merging: ";
			for(unsigned int k=0;k<(*ptree[max_miIndex]).size();k++)
			{
				out << (*ptree[max_miIndex]).at(k) << ' ';
			}
			update_cells(*ptree[max_miIndex],pMxScC);	//del rest of merged cells
			pMxScC = NULL;
			mvCell[i]->mIsMerged = true;
			out << std::endl;
			out << mvCell[i]->mLabel << ' ' << mvCell[i]->mScore << std::endl;
			--i;	//still need to build rag for new merged cell
		}
		out << std::endl;
		out.close();
	}
	for(unsigned int j=0; j<ttree.size(); j++)
		delete ttree[j];
}

/**	@brief	score cell
*	@param	c	cell to score
*	@param	v	cells that cell c is merged from, if not merged will be empty
*/
void
cell_mgr::score_cell_mah(cell *c,std::vector<int> v)
{
	vnl_matrix<double> x(mMtrxNeuron.cols(),1);	//temp vector to hold resul

	update_volume( c );	//load features
	x(0,0) = c->mVolume;

	update_convexity( c );
	x(1,0) = c->mConvexity;	

	update_shape_fact( c );
	x(2,0) = c->mShapeFact;

	update_eccentricity( c ,v);
	x(3,0) = c->mEccentricity;

	update_bend_eng( c );
	x(4,0) = c->mBendEng;

	update_vol_grad( c );
	x(5,0) = c->mVolGrad;

	update_bound_grad( c );
	x(6,0) = c->mBoundGrad;

	update_texture( c );
	x(7,0) = c->mTexture;		

	update_avg_int( c );
	x(8,0) = c->mAvgInt ;

	x(9,0) = c->mRadVar;

	update_ints_ratio( c );
	x(10,0) = c->mBoundIntsRatio;

	update_per_nbr( c , v );
	x(11,0) = c->mPerNbr;

	c->mScore = sqrt( ( ( x - mMeanN ).transpose() * mMtrxInvCovNeuron * ( x - mMeanN ) ).get(0,0) );
	c->mClass = 1;

	if( sqrt( ( ( x - mMeanG ).transpose() * mMtrxInvCovGlia * ( x - mMeanG ) ).get(0,0) ) < c->mScore )
	{
		c->mScore = sqrt( ( ( x - mMeanG ).transpose() * mMtrxInvCovGlia * ( x - mMeanG ) ).get(0,0) );
		c->mClass = 2;
	}
}

bool
cell_mgr::check_shape_dist( cell* c , const double &conv , const double &bend , const double &shape , const double &pnbr , const double &vol )
{
	if( c == NULL )
		return false;

	//check each feature

	update_volume(  c );
	if( c->mVolume <= vol )
		return false;

	update_per_nbr( c , bla);
	update_convexity( c );
	update_shape_fact( c );
	update_bend_eng( c );

	double celldist = sqrt( c->mPerNbr * c->mPerNbr + (1-c->mConvexity)*(1-c->mConvexity) + c->mShapeFact * c->mShapeFact + c->mBendEng*c->mBendEng );

	if( celldist > sqrt( pnbr*pnbr + (1-conv)*(1-conv) + shape*shape + bend*bend ) )
		return false;

	return true;
}

void
cell_mgr::check_score_cells( const double &per )
{
	double gscore = 0.0,nscore =0.0;
	double ncount = 0, gcount=0;
	for(unsigned int i=0 ; i < mvCell.size(); i++ )	//go through each cell and find train cells
	{ 
		if( mvCell[i] == NULL )
			continue;
		if( mvCell[i]->mClass == 1 )
		{
			nscore += mvCell[i]->mScore;
			ncount++;
		}
		else if( mvCell[i]->mClass == 2 )
		{
			gscore += mvCell[i]->mScore;
			gcount++;
		}
	}

	nscore /= ncount;
	gscore /= gcount;
	nscore *= per;
	gscore *= per;
	for(unsigned int i=0 ; i < mvCell.size(); i++ )	//go through each cell and find train cells
	{ 
		if( mvCell[i] == NULL )
			continue;
		if( mvCell[i]->mClass == 1 &&  mvCell[i]->mScore < nscore)
		{
			mvCell[i]->mClass = -1;
		}
		else if( mvCell[i]->mClass == 2 &&  mvCell[i]->mScore < gscore)
		{
			mvCell[i]->mClass = -1;
		}
	}
}

/**	@brief scores foreground
*/
void
cell_mgr::score_foreground_mah( std::string param )
{
	double totalScore = 0.0;
	double h = 4.0/log(double(mMtrxNeuron.rows()));
	vnl_matrix<double> x(mMtrxNeuron.cols(),1);	//temp vector to hold resul
	vnl_matrix<double> t(mMtrxNeuron.cols(),1);

	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] == NULL )
			continue;

		update_volume( mvCell[i] );	//load features
		x(0,0) = mvCell[i]->mVolume;

		update_convexity( mvCell[i] );
		x(1,0) = mvCell[i]->mConvexity;	

		update_shape_fact( mvCell[i] );
		x(2,0) = mvCell[i]->mShapeFact;

		update_eccentricity( mvCell[i] ,bla );
		x(3,0) = mvCell[i]->mEccentricity;

		update_bend_eng( mvCell[i] );
		x(4,0) = mvCell[i]->mBendEng;

		update_vol_grad( mvCell[i] );
		x(5,0) = mvCell[i]->mVolGrad;

		update_bound_grad( mvCell[i] );
		x(6,0) = mvCell[i]->mBoundGrad;

		update_texture( mvCell[i] );
		x(7,0) = mvCell[i]->mTexture;		

		update_avg_int( mvCell[i] );
		x(8,0) = mvCell[i]->mAvgInt ;

		x(9,0) = mvCell[i]->mRadVar;

		update_ints_ratio( mvCell[i] );
		x(10,0) = mvCell[i]->mBoundIntsRatio;

		update_per_nbr( mvCell[i] , bla );
		x(11,0) = mvCell[i]->mPerNbr;

		update_depth( mvCell[i] );
		x(12,0) = mvCell[i]->mDepth;

		x(13,0) = (double) mvCell[i]->mvBoundPix2D.size();


		mvCell[i]->mScore = 0.0;

		double tempG = 0.0, tempN = 0.0;

		h = 4.0/log(double(mMtrxNeuron.rows()));
		for(unsigned int j=0; j<mMtrxNeuron.rows(); j++ )	//and calc score relative to other cells
		{
			t = x - mMtrxNeuron.get_n_rows(j,1).transpose();
			tempN += exp( - ( ( t.transpose()*mMtrxInvCovNeuron*t ).get(0,0) / (2.0*h*h)  ) );
		}
		tempN /= double(mMtrxNeuron.rows()) * sqrt(2.0 * M_PI * h * h * mDetNeuron);

		h = 4.0/log(double(mMtrxGlia.rows()));
		for(unsigned int j=0; j<mMtrxGlia.rows(); j++ )	//and calc score relative to other cells
		{
			t = x - mMtrxGlia.get_n_rows(j,1).transpose();
			tempG += exp( - ( ( t.transpose()*mMtrxInvCovGlia*t ).get(0,0) / (2.0*h*h)  ) );
		}
		tempG /= double(mMtrxGlia.rows()) * sqrt(2.0 * M_PI * h * h * mDetGlia);

		if( mvCell[i]->mClass == 1 )
		{
			mvCell[i]->mScore = tempN - tempG;
		}
		else if( mvCell[i]->mClass == 2 )
		{
			mvCell[i]->mScore = tempG - tempN;
		}
			
		if( (mvCell[i]->mScore >= 0 && mvCell[i]->mScore < 1.7E308 ) )
		{
			totalScore += mvCell[i]->mScore * mvCell[i]->mVolume;
		}
		

	}
	std::ofstream oscore( "_score.txt",std::ofstream::app );
	oscore << param << ' ' << totalScore <<  std::endl;
	oscore.close();
}

/**	@brief	update avg vol gradient variance
*	@param	c cell
*/
void
cell_mgr::update_bound_grad_var(cell *c)
{
	c->mBoundGradVar = 0.0;
	if(c->mvBoundPix2D.empty())
		return;
	mItType it(mImgpGradient,mImgpGradient->GetLargestPossibleRegion());	//image

	update_bound_grad(c);

	for(unsigned int i=0;i<c->mvBoundPix2D.size();i++)	//unbaised variance
	{

		miIndex[0] = c->mvBoundPix2D[i]->x_;
		miIndex[1] = c->mvBoundPix2D[i]->y_;
		miIndex[2] = c->mvBoundPix2D[i]->z_;
		it.SetIndex(miIndex);
		c->mBoundGradVar += ( double( it.Get() ) - c->mBoundGrad ) * ( double( it.Get() ) - c->mBoundGrad );
		//it.Set(AvgVar);
	}		
	c->mBoundGradVar /= double( c->mvBoundPix2D.size() - 1 );
}

inline double
cell_mgr::sigmoid( const double &a, const double &b, const double &net )
{
	return a * ( exp( b * net ) - exp( -b * net ) ) / ( exp( b * net ) + exp( -b * net ) );
}

/**	@brief	loads cells to be merged to create a training sample
*	input file should have the cell labels followed by a -1 at the end
*/

void
cell_mgr::del_some_cells(void)
{
	std::vector<int> v;
	std::ifstream in("cells_to_del.txt");
	while( !in.eof() )
	{
		int temp;
		in >> temp;
		std::cout << "Deleting cell: " << temp << std::endl;
		if( mvCell[temp - 1] != NULL )
		{
			if( mvCell[temp - 1]->mIsMerged == true )
				delete mvCell[temp - 1];
			mvCell[temp - 1] = NULL;
		}
	}
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if(mvCell[i] != NULL)
			cell_mgr::update_center(mvCell[i]);		//make sure to update the center for
	}
}

/**	@brief	save labeled image
*/
void
cell_mgr::save_labeled(std::string &file)
{
	std::ifstream in(std::string( "_label" + file ).c_str() );
	if( in.is_open() )
	{
		in.close();
		return;
	}
	mItType itL(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	file = "_label" + file;

	std::ofstream out( file.c_str() );

	for(itL.GoToBegin();!itL.IsAtEnd();++itL)
	{
		out << itL.Get() << ' ';
	}

	for(unsigned int i=0;i<mvCell.size(); i++ )
	{
		out << mvCell[i]->mClass << ' ';
	}
}

/**	@brief project cells onto one plane
*/
void
cell_mgr::project_cells(void)
{
	std::vector<int> labels(mvCell.size(),0);
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] != NULL )
		{
			labels[i] = mvCell[i]->mClass ;
		}
		
	}
	mItType itLabel(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	for ( itLabel.GoToBegin(); !itLabel.IsAtEnd(); ++itLabel)
	{
		itLabel.Set( -1 );
	}
	mInputImageFileType::IndexType ind;
		//loop through all cells
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if(mvCell[i] == NULL)
			continue;

		for(unsigned int j=0;j<mvCell[i]->mvPix.size();j++)
		{
			ind[0] = mvCell[i]->mvPix[j]->x_;
			ind[1] = mvCell[i]->mvPix[j]->y_;
			ind[2] = mvCell[i]->mvPix[j]->z_;
			mImgpLabel->SetPixel(ind,mvCell[i]->mLabel);
		}
		mvCell[i] = NULL;
	}

	mvCell.clear();

	mCitType cit(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	mImageType::Pointer imProj = mImageType::New();

	imProj->SetRegions(mImgpLabel->GetLargestPossibleRegion());
	imProj->CopyInformation( mImgpLabel );
	imProj->Allocate();

	mItType itp(imProj,imProj->GetLargestPossibleRegion());

	for ( itp.GoToBegin(); !itp.IsAtEnd(); ++itp)
	{
		itp.Set( -1 );
	}
	for(int x=0; x<mImageSize[0]; x++ )
	{
		for(int y=0; y<mImageSize[1]; y++ )
		{
			int max = -1;
			for(int z=0; z<mImageSize[2]; z++ )
			{
				miIndex[0] = x;
				miIndex[1] = y;
				miIndex[2] = z;
				cit.SetIndex(miIndex);
				max = ((cit.Get() > max) ? cit.Get() : max);
			}
			miIndex[0] = x;
			miIndex[1] = y;
			miIndex[2] = 0;
			itp.SetIndex(miIndex);
			itp.Set(max);
		}
	}

	mCitType citp(imProj,imProj->GetLargestPossibleRegion());
		//sort labeled image to init cells
	std::sort(mvPix.begin(),mvPix.end(),compare(imProj,mImageSize));

	//get nbr iterator
	typedef itk::ConstNeighborhoodIterator< mImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 1;

	NeighborhoodIteratorType it( radius, imProj , imProj->GetRequestedRegion() );

	NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
	NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
	NeighborhoodIteratorType::OffsetType offset3 = {{0,1,0}};
	NeighborhoodIteratorType::OffsetType offset4 = {{0,-1,0}};
	NeighborhoodIteratorType::OffsetType offset5 = {{0,0,1}};
	NeighborhoodIteratorType::OffsetType offset6 = {{0,0,-1}};

	int istart;
	//find first non background cell
	for(unsigned istart=0;istart<mvPix.size();istart++)
	{
		miIndex[0] = mvPix[istart]->x_;
		miIndex[1] = mvPix[istart]->y_;
		miIndex[2] = mvPix[istart]->z_;
		citp.SetIndex(miIndex);
		if(citp.Get() != -1)
			break;
	}
	//keep track of cell labels
	int count = 1;
	//label of previous pix, to know when we start pix of a new cell
	int prev = -1;
	cell *c=NULL;
	for(unsigned int i=istart;i<mvPix.size();i++)
	{
		if( mvPix[i]->z_ != 0 )
			continue;
		miIndex[0] = mvPix[i]->x_;
		miIndex[1] = mvPix[i]->y_;
		miIndex[2] = mvPix[i]->z_;
		it.SetLocation(miIndex);
		//check if new label found
		if(prev != it.GetCenterPixel() )
		{
			//check to see if c has a cell and push it back on the vector
			if(c != NULL)
				mvCell.push_back(c);
			//create new cell
			c = new cell(it.GetCenterPixel());
			c->mClass = labels[ it.GetCenterPixel() - 1 ];
			//set old lebel as prev
			prev = (int)it.GetCenterPixel();
			count++;
		}

		if( mvPix[i]->x_ == mImageSize[0] - 1  || mvPix[i]->x_ == 0 || mvPix[i]->y_ == 0 || mvPix[i]->y_ == mImageSize[1] - 1 )
		{
			add_pix_bound_2d(mvPix[i],c);
			add_pix_vol_3d(mvPix[i],c);
			continue;
		}
		//check 2d nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
			|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4) )
		{
			add_pix_bound_2d(mvPix[i],c);
		}
		else
			add_pix_vol_2d(mvPix[i],c);

		//check 3d nbr labels, if all same as pix then add pix as volume otherwise as bound
		if( it.GetCenterPixel() != it.GetPixel(offset1) || it.GetCenterPixel() != it.GetPixel(offset2)
			|| it.GetCenterPixel() != it.GetPixel(offset3) || it.GetCenterPixel() != it.GetPixel(offset4)
			|| it.GetCenterPixel() != it.GetPixel(offset5) || it.GetCenterPixel() != it.GetPixel(offset6))
		{
			add_pix_bound_3d(mvPix[i],c);
		}
		else
			add_pix_vol_3d(mvPix[i],c);
	}
	
	//add in last cell
	if(c != NULL)
	{	
		c->mIsMerged = true;
		mvCell.push_back(c);
	}
}

/**	@brief	remove any cell touching border
*/
void
cell_mgr::remove_border_cells(void)
{
	for(unsigned int i=0;i<mvCell.size();i++)
	{
		if( mvCell[i] == NULL )
			continue;
		for(unsigned int j=0;j<mvCell[i]->mvPix.size();j++)
		{
			if( mvCell[i]->mvPix[j]->x_ == 0 || mvCell[i]->mvPix[j]->x_ == mImageSize[0] - 1 || mvCell[i]->mvPix[j]->y_ == 0 || mvCell[i]->mvPix[j]->y_ == mImageSize[1] - 1 )
			{
				if( mvCell[i]->mIsMerged )
					delete mvCell[i];
				mvCell[i] = NULL;
				break;
			}
		}
	}
}

/**	@brief	save training data
*/
void
cell_mgr::save_train_data(void)
{
	std::ofstream outglia("glia.txt");
	std::ofstream outneuron("neuron.txt");
	std::ofstream outall("allcell.txt");
	
	outneuron << mMtrxNeuron;
	outglia << mMtrxGlia;
	outall << mNoclass;
	outall << "Neuron " << mDetNeuron <<std::endl;
	outall << mMtrxInvCovNeuron;
	outall << "Glia " << mDetGlia <<std::endl;
	outall << mMtrxInvCovGlia;
}

/**	@brief	remove cells not touching medain planes
*/
void
cell_mgr::remove_non_median_cells(void)
{
	//calc lower and upper planes, median plane + 10% z depth are valid med planes
	int lower = (int) ( (double)mImageSize[2] / 2.0 - 0.1*(double)mImageSize[2] );
	int upper = (int) ( (double)mImageSize[2] / 2.0 + 0.1*(double)mImageSize[2] );
	std::cout << "lower: " << lower << std::endl;
	std::cout << "upper: " << upper << std::endl;
	for(unsigned int i=0;i<mvCell.size();i++)	//check each cell if touching med planes
	{
		if( mvCell[i] == NULL )
			continue;
		update_depth(mvCell[i]);
		if( mvCell[i]->mStartPlane < lower && mvCell[i]->mStartPlane + mvCell[i]->mDepth < lower )
		{
			if( mvCell[i]->mIsMerged )
					delete mvCell[i];
			mvCell[i] = NULL;
		}
		else if( mvCell[i]->mStartPlane > upper )
		{
			if( mvCell[i]->mIsMerged )
					delete mvCell[i];
			mvCell[i] = NULL;
		}
	}
}
