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

/** @file wtd_lg.cpp
*   @brief main class for watershed algorithm
*
*   @author Maciej Wotjon
*/



#include "wtd_lg.h"
#include "itkImage.h"
#include "pix_t.h"
#include <vector>
#include <numeric>

#include <math.h>

wtd_lg::wtd_lg(int ims[],std::vector<pix_t*> p,mImageType::Pointer ip)
{
	//set image size
	mImageSize[0] = ims[0];
	mImageSize[1] = ims[1];
	mImageSize[2] = ims[2];
	//init distance 
	for(unsigned int i=0;i<p.size();i++)
	{
		mvDist.push_back(0);
	}
	//init pix
	mvPix = p;
	//create labeled image
	mImgpFilt = mImageType::New();
	mImgpFilt = ip;
	mImgpLabel = mImageType::New();
	mImgpLabel->SetRegions( mImgpFilt->GetLargestPossibleRegion() );
	mImgpLabel->CopyInformation( mImgpFilt );
	mImgpLabel->Allocate();
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	for(it.GoToBegin();!it.IsAtEnd();++it)
	{
		it.Set(-1);
	}
}

wtd_lg::~wtd_lg()
{
}

/**	@brief setsup and runs watershed alg
*	@param value_min
*	@param value_max
*	@param step_size
*	@param set_nb6_flag
*/
void
wtd_lg::run_watershed(short value_min,short value_max, short step_size, unsigned char set_nb6_flag)
{
	nb6_flag=set_nb6_flag;

	mCitType cit(mImgpFilt,mImgpFilt->GetLargestPossibleRegion());

	// 3D watershed algorithm
	std::vector<pix_t*> pix;

	for( int z=0;z<mImageSize[2];++z )
	{
		for( int y=0;y<mImageSize[1];++y )
		{
			for( int x=0;x<mImageSize[0];++x )
			{
				miIndex[0] = x;
				miIndex[1] = y;
				miIndex[2] = z;
				cit.SetIndex(miIndex);
				if(cit.Get()>=value_min && cit.Get()<=value_max)
				{
					pix.push_back(mvPix[id_image(x,y,z)]);
				}
			}
		}
	}
	// NOTE ALEX: unused variable
        // short max_label=wtd3d(pix, value_min, value_max, step_size );

	pix.clear();
}

/**	@brief get neighbors of the given pixel depending on the flags set
*	@param p nbs found around pix p
*	@param nbs nbrs of p
*	@param tdflag 2dflag -> true if 4 2d nbrs, flase -> 6 3d nbrs
*/
void 
wtd_lg::get_neighbors( pix_t* p, std::vector<pix_t*>& nbs , bool tdflag)
{
	nbs.clear();
	pix_t* ptr=NULL;
	if(tdflag)
	{
		if(nb6_flag==1)
		{
			nbs.reserve(6);
			if( p->x_>0 )
			{
				ptr=mvPix[id_image(p->x_-1,  p->y_,   p->z_)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
			if( p->x_<mImageSize[0]-1 )
			{
				ptr=mvPix[id_image(p->x_+1,  p->y_,   p->z_)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
			if( p->y_>0 )
			{
				ptr=mvPix[id_image(p->x_,p->y_-1, p->z_)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
			if( p->y_<mImageSize[1]-1 )
			{
				ptr=mvPix[id_image(p->x_,    p->y_+1, p->z_)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
			if( p->z_>0 )
			{
				ptr=mvPix[id_image(p->x_,    p->y_,   p->z_-1)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
			if( p->z_<mImageSize[2]-1 )
			{
				ptr=mvPix[id_image(p->x_,    p->y_,   p->z_+1)];
				if(ptr!=NULL) 
					nbs.push_back(ptr);
			}
		}
		else
		{
			nbs.reserve(26);
			for( int z=p->z_-1; z<=p->z_+1; ++z )
			{
				for( int y=p->y_-1; y<=p->y_+1; ++y )
				{
					for( int x=p->x_-1; x<=p->x_+1; ++x )
					{
						if( x<0 || x>mImageSize[0]-1 || y<0 || y>mImageSize[1]-1 || z<0 || z>mImageSize[2]-1 ) 
							continue;
						if( x==p->x_ && y==p->y_ && z==p->z_ )
							continue;
						ptr=mvPix[id_image(x,y,z)];
						if(ptr!=NULL) 
							nbs.push_back(ptr);
					}
				}
			}
		}
	}
	else	//get 4 nbrs in 2d
	{
		nbs.reserve(4);
		if( p->x_>0 )
		{
			ptr=mvPix[id_image(p->x_-1,  p->y_,   p->z_)];
			if(ptr!=NULL) 
				nbs.push_back(ptr);
		}
		if( p->x_<mImageSize[0]-1 )
		{
			ptr=mvPix[id_image(p->x_+1,  p->y_,   p->z_)];
			if(ptr!=NULL) 
				nbs.push_back(ptr);
		}
		if( p->y_>0 )
		{
			ptr=mvPix[id_image(p->x_,p->y_-1, p->z_)];
			if(ptr!=NULL) 
				nbs.push_back(ptr);
		}
		if( p->y_<mImageSize[1]-1 )
		{
			ptr=mvPix[id_image(p->x_,    p->y_+1, p->z_)];
			if(ptr!=NULL) 
				nbs.push_back(ptr);
		}
	}
}

/**	@brief after running the watershed procedure, obtain the real tesselation by removing the WSHED pixels
*/
void 
wtd_lg::remove_wtd(std::vector<pix_t*>& pix)
{
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	mItType it2(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	std::vector<pix_t*> nbs;
	int num=1;
	while(num>0)
	{
		num=0;
		for(unsigned int j=0;j<pix.size();j++)
		{
			pix_t* ptr=pix[j];
			if(ptr==NULL) 
				continue;
			miIndex[0] = ptr->x_;
			miIndex[1] = ptr->y_;
			miIndex[2] = ptr->z_ ;
			it.SetIndex(miIndex);
			if(it.Get()!=WSHED) 
				continue;

			get_neighbors(ptr, nbs );
			for(unsigned int i=0;i<nbs.size();i++ )
			{
				if( nbs[i]!=NULL)
				{
					miIndex[0] = nbs[i]->x_;
					miIndex[1] = nbs[i]->y_;
					miIndex[2] = nbs[i]->z_;
					it2.SetIndex(miIndex);
					if(it2.Get()>0)
					{
						it.Set(it2.Get());
						num++;
					}
				}
			}
		}
	}
}

/**	@brief main routinue for 3D watershed
*/
short 
wtd_lg::wtd3d(std::vector<pix_t*>& pix, short h_min, short h_max, short step_size)
{
	mCitType cit(mImgpFilt,mImgpFilt->GetLargestPossibleRegion());
	mItType it(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());
	mItType it2(mImgpLabel,mImgpLabel->GetLargestPossibleRegion());

	short current_label = 0; //minimum label ???
	unsigned char current_dist;
	pix_t* fictitious_pix = new pix_t();

	//sorting
	std::sort( pix.begin(), pix.end(),compare(mImgpFilt,mImageSize));

	int size = (int)pix.size();

	//store the neighbors of a pixel
	std::vector<pix_t*> nbs;

	//queue structure
	std::deque<pix_t*> q;

	int start_i=0, end_i=0;		
	miIndex[0] = pix[end_i]->x_;
	miIndex[1] = pix[end_i]->y_;
	miIndex[2] = pix[end_i]->z_;
	cit.SetIndex(miIndex);
	for( short h=h_min; h<=h_max; h+=step_size )
	{
		start_i = end_i;

		for(;end_i<size && cit.Get()<=h;end_i++)
		//while(end_i<size && cit.Get()<=h )
		{	
			//end_i++;
			miIndex[0] = pix[end_i]->x_;
			miIndex[1] = pix[end_i]->y_;
			miIndex[2] = pix[end_i]->z_;
			cit.SetIndex(miIndex);
		}

		for( int i=start_i;i<end_i;i++ )
		{
			miIndex[0] = pix[i]->x_;
			miIndex[1] = pix[i]->y_;
			miIndex[2] = pix[i]->z_;
			it.SetIndex(miIndex);
			it.Set(MASK);
			get_neighbors(pix[i], nbs );
			bool flag=false;
			for(unsigned int j=0;j<nbs.size();++j )
			{
				miIndex[0] = nbs[j]->x_;
				miIndex[1] = nbs[j]->y_;
				miIndex[2] = nbs[j]->z_;
				it.SetIndex(miIndex);
				if( it.Get()>0 || it.Get() == WSHED )
				{
					flag=true;
					break;
				}
			}
			if(flag)
			{
				mvDist[ id_image( pix[i]->x_ , pix[i]->y_ , pix[i]->z_ ) ] = 1;
				q.push_back( pix[i] );
			}
		}
		current_dist = 1;
		q.push_back( fictitious_pix );
		while( true )
		{
			pix_t* p = q.front();
			q.pop_front();
			if( p->z_ == -1 )
			{
				if( q.empty() ) break;
				else
				{
					q.push_back( fictitious_pix );
					current_dist++;
					p = q.front();
					q.pop_front();
				}
			}

			get_neighbors(p, nbs );
			for(unsigned int i=0;i<nbs.size();++i )
			{
				miIndex[0] = nbs[i]->x_;
				miIndex[1] = nbs[i]->y_;
				miIndex[2] = nbs[i]->z_;
				it.SetIndex(miIndex);

				if( mvDist[ id_image( nbs[i]->x_ , nbs[i]->y_ , nbs[i]->z_ ) ]<current_dist && ( it.Get()>0 
					|| it.Get()==WSHED ) )
				{
					miIndex[0] = p->x_;
					miIndex[1] = p->y_;
					miIndex[2] = p->z_;
					it2.SetIndex(miIndex);
					if( it.Get()>0 )
					{
						if( it2.Get()==MASK || it2.Get() == WSHED )
						{

							it2.Set( it.Get() );
						}
						else if( it2.Get() != it.Get() )
						{
							it2.Set(WSHED);
						}
					}
					else if( it2.Get() == MASK )
						it2.Set( WSHED );
				}
				else if( it.Get()==MASK && mvDist[ id_image( nbs[i]->x_ , nbs[i]->y_ , nbs[i]->z_ ) ] == 0 )
				{
					mvDist[ id_image( nbs[i]->x_ , nbs[i]->y_ , nbs[i]->z_ ) ] = current_dist+1;
					q.push_back( nbs[i] );
				}
			}
		}

		for( int i=start_i;i<end_i;i++ )
		{
			mvDist[ id_image( pix[i]->x_ , pix[i]->y_ , pix[i]->z_ ) ] = 0;
			miIndex[0] = pix[i]->x_;
			miIndex[1] = pix[i]->y_;
			miIndex[2] = pix[i]->z_;
			it.SetIndex(miIndex);
			if( it.Get() == MASK )
			{
				current_label++;
				q.push_back( pix[i] );
				it.Set(current_label);
				while( !q.empty() )
				{
					pix_t* p1 = q.front();
					q.pop_front();
					get_neighbors(p1, nbs );
					for(unsigned int j=0;j<nbs.size();++j )
					{
						miIndex[0] = nbs[j]->x_;
						miIndex[1] = nbs[j]->y_;
						miIndex[2] = nbs[j]->z_;
						it2.SetIndex(miIndex);
						if( it2.Get() == MASK )
						{
							q.push_back( nbs[j] );
							it2.Set( current_label );
						}
					}
				}

			}
		}
	}

	//remove wtd
	remove_wtd(pix);

	//release memory
	delete fictitious_pix;
	nbs.clear();
	q.clear();

	return current_label;
}

/**	@brief returns labeled image
*	@return ITK image pointer to labeled image
*/
wtd_lg::mImageType::Pointer 
wtd_lg::getlabel(void)
{
	return mImgpLabel;
}
