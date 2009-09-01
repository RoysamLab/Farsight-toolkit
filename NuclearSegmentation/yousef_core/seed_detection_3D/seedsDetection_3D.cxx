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

//seedsDetection_2D.cxx

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include <stdio.h>
#include<stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>

#include "itkImage.h"
#include "itklaplacianrecursivegaussianimagefilternew.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkDanielssonDistanceMapImageFilter.h> 
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"

//added by Yousef on 8/26/2009
#include "itkExtractImageFilter.h"


////
#include "itkImageFileWriter.h"

using namespace std;


typedef    unsigned short     MyInputPixelType;
typedef itk::Image< MyInputPixelType,  3 >   MyInputImageType;
typedef itk::Image< MyInputPixelType,  2 >   MyInputImageType2D;

int detect_seeds(itk::SmartPointer<MyInputImageType>, int , int , int, const double, float*, int);
//int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, int r, int c, int z,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, int unsigned short* dImg, int* minIMout, int UseDistMap);
int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, int rmin, int rmax, int cmin, int cmax, int zmin, int zmax,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, int unsigned short* dImg, int* minIMout, int UseDistMap);
float get_maximum_3D(float* A, int r1, int r2, int c1, int c2, int z1, int z2, int R, int C);
unsigned short get_maximum_3D(unsigned short* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C);
void Detect_Local_MaximaPoints_3D(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1, unsigned short* bImg);
int distMap(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG);
int distMap_SliceBySlice(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG);
MyInputImageType2D::Pointer extract2DImageSlice(itk::SmartPointer<MyInputImageType> im, int plane, int slice);
MyInputImageType::Pointer extract3DImageRegion(itk::SmartPointer<MyInputImageType> im, int sz_x, int sz_y, int sz_z, int start_x, int start_y, int start_z);
void estimateMinMaxScales(itk::SmartPointer<MyInputImageType> im, unsigned short* distIm, double* minScale, double* maxScale, int r, int c, int z);
int computeMedian(std::vector< std::vector<unsigned short> > scales, int cntr);

int Seeds_Detection_3D( float* IM, float** IM_out, unsigned short** IM_bin, int r, int c, int z, double sigma_min, double sigma_max, double scale_xy, double scale_z, int sampl_ratio, unsigned short* bImg, int UseDistMap, int* minIMout)
{	
	//Create an itk image
	MyInputImageType::Pointer im;
	im = MyInputImageType::New();
	MyInputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    im->SetOrigin( origin );

    MyInputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    MyInputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z
  
    MyInputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = sampl_ratio; //spacing along z

    im->SetRegions( region );
	im->SetSpacing(spacing);
    im->Allocate();
    im->FillBuffer(0);
	im->Update();	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< MyInputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	
	unsigned short* dImg = NULL;	
	if(UseDistMap == 1)
	{
		for(int i=0; i<r*c*z; i++)
		{				
			if(bImg[i]>0)
				iterator1.Set(0.0);//IM[i]);
			else
				iterator1.Set(255.0);
			++iterator1;
		
		}
		std::cout<<"Computing distance transform...";
		dImg = (unsigned short *) malloc(r*c*z*sizeof(unsigned short));
		if(!dImg)
		{
			std::cerr<<"Failed to allocate memory for the distance image"<<std::endl;
			return 0;
		}
		//distMap(im, r, c, z,dImg);
		distMap_SliceBySlice(im, r, c, z,dImg);
		std::cout<<"done"<<std::endl;
		iterator1.GoToBegin();		
	}
	

	for(int i=0; i<r*c*z; i++)
	{				
	    iterator1.Set((unsigned short)IM[i]);	
		++iterator1;
	}

	//try this: estimate the min and max scales
	if(UseDistMap == 1)
	{
		std::cout<<"Estimating parameters..."<<std::endl;
		estimateMinMaxScales(im, dImg, &sigma_min, &sigma_max, r, c, z);
		scale_xy = sigma_min+1;
		scale_z = ceil(scale_xy / sampl_ratio);
		std::cout<<"    Minimum scale = "<<sigma_min<<std::endl;
		std::cout<<"    Maximum scale = "<<sigma_max<<std::endl;
		std::cout<<"    Clustering Resolution = "<<scale_xy<<std::endl;
	}
	//By Yousef (8/28/2009)
	//In some situations the image is very larg and we cannot allocate memory for the LoG filter (20xthe size of the image in bytes)
	//In such cases, we can divide the image into small tiles, process them independently, and them combine the results
	int block_divisor;	
	//see if we have enought memory for the LoG step
	//approximately, we need (20~21)ximage size in bytes
	//try to allocate memory for an unsigned char* of the 23ximage size
	unsigned char *tmpp = (unsigned char*) malloc(23*r*c*z*sizeof(unsigned char));
	if(tmpp)
		block_divisor = 1;
	else
		block_divisor = 2;
	free(tmpp); //delete it
	tmpp = NULL;
	
	
	int min_x, min_y, max_x, max_y;
	for(int i=0; i<r; i+=r/block_divisor)
	{
		for(int j=0; j<c; j+=c/block_divisor)
		{						
			min_x = j; 
			max_x = 40+(int)j+c/block_divisor; //40 is the size of the overlapping between the two tiles along x
			min_y = i;
			max_y = 40+(int)i+r/block_divisor; //40 is the size of the overlapping between the two tiles along y
			if(max_x >= c)
				max_x = c-1;
			if(max_y >= r)
				max_y = r-1;

			//Create an itk image to hold the sub image (tile) being processing
			MyInputImageType::Pointer im_Small = extract3DImageRegion(im, max_x-min_x+1, max_y-min_y+1, z, min_x, min_y, 0);
						
			//By Yousef (8/27/2009): multi-scale LoG is done in one function now
			multiScaleLoG(im_Small, r, c, z, min_y, max_y, min_x, max_x, 0, z-1, sigma_min, sigma_max, IM, sampl_ratio, dImg, minIMout, UseDistMap);
			//
		}
	}
	
		
	free(dImg);
   
	IM_out[0] = new float[r*c*z];
	if(!IM_out[0])
	{
		std::cerr<<"could not allocate memory for the LoG response image"<<std::endl;
		return 0;
	}
	for(int i=0; i<r*c*z; i++)
	{
		IM_out[0][i] = IM[i];
	}
	//Detect the seed points (which are also the local maxima points)	
	std::cout<<"Detecting Seeds"<<std::endl;
	IM_bin[0] = new unsigned short[r*c*z];
	Detect_Local_MaximaPoints_3D(IM_out[0], r, c, z, scale_xy, scale_z, IM_bin[0], bImg);	
	
	return 1;
}


int detect_seeds(itk::SmartPointer<MyInputImageType> im, int r, int c, int z,const double sigma, float* IMG, int sampl_ratio)
{
  
  //  Types should be selected on the desired input and output pixel types.
  typedef    float     OutputPixelType;


  //  The input and output image types are instantiated using the pixel types.
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;


  //  The filter type is now instantiated using both the input image and the
  //  output image types.
  typedef itk::LaplacianRecursiveGaussianImageFilterNew<MyInputImageType, OutputImageType >  FilterType;
  FilterType::Pointer laplacian = FilterType::New();
  

  //  The option for normalizing across scale space can also be selected in this filter.
  laplacian->SetNormalizeAcrossScale( true );

  //  The input image can be obtained from the output of another
  //  filter. Here the image comming from the calling function is used as the source
  laplacian->SetInput( im);
  

  
  //  It is now time to select the $\sigma$ of the Gaussian used to smooth the
  //  data.  Note that $\sigma$ must be passed to both filters and that sigma
  //  is considered to be in millimeters. That is, at the moment of applying
  //  the smoothing process, the filter will take into account the spacing
  //  values defined in the image.
  //
  laplacian->SetSigma(sigma);
 
  //  Finally the pipeline is executed by invoking the \code{Update()} method.
  //
 try
    {
    laplacian->Update();
    }
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 
 
  //   Copy the resulting image into the input array
  long int i = 0;
  typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType;
  IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
  while ( i<r*c*z)
  {
    IMG[i] = /*sigma*sigma*/iterate.Get();
    ++i;
	++iterate;
  }

  return EXIT_SUCCESS;
}

int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, int rmin, int rmax, int cmin, int cmax, int zmin, int zmax,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, int unsigned short* dImg, int* minIMout, int UseDistMap)
{
  
  //  Types should be selected on the desired input and output pixel types.
  typedef    float     OutputPixelType;
  //  The input and output image types are instantiated using the pixel types.
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;


  for(int sigma=sigma_min; sigma<=sigma_max; sigma++)
  {
    std::cout<<"Processing scale "<<sigma<<"...";
	//  The filter type is now instantiated using both the input image and the
	//  output image types.
	typedef itk::LaplacianRecursiveGaussianImageFilterNew<MyInputImageType, OutputImageType >  FilterType;
	FilterType::Pointer laplacian = FilterType::New();
	//  The option for normalizing across scale space can also be selected in this filter.
	laplacian->SetNormalizeAcrossScale( true );

	//  The input image can be obtained from the output of another
	//  filter. Here the image comming from the calling function is used as the source
	laplacian->SetInput( im);
  
	//  It is now time to select the $\sigma$ of the Gaussian used to smooth the
	//  data.  Note that $\sigma$ must be passed to both filters and that sigma
	//  is considered to be in millimeters. That is, at the moment of applying
	//  the smoothing process, the filter will take into account the spacing
	//  values defined in the image.
	//
	laplacian->SetSigma(sigma);
 
	//  Finally the pipeline is executed by invoking the \code{Update()} method.
	//
	try
	{
		laplacian->Update();
    }
	catch( itk::ExceptionObject & err ) 
    { 
		std::cout << "ExceptionObject caught !" << std::endl; 
		std::cout << err << std::endl; 
		return EXIT_FAILURE;
    } 
 
	//   Copy the resulting image into the input array
	//long int i = 0;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType;
	IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
	long int II;
	for(int k1=zmin; k1<=zmax; k1++)
	{
		for(int i1=rmin; i1<=rmax; i1++)	
		{			
			for(int j1=cmin; j1<=cmax; j1++)
			{
				II = (k1*r*c)+(i1*c)+j1;
				if(sigma==sigma_min)
				{
					IMG[II] = /*sigma*sigma*/iterate.Get();			
				}
				else
				{
					float lgrsp = iterate.Get();
					if(UseDistMap == 1)
					{
						if(sigma<=dImg[II])
						{
							IMG[II] = (IMG[II]>=lgrsp)? IMG[II] : lgrsp;				
							if(IMG[II]<minIMout[0])
								minIMout[0] = IMG[II];
						}
					}
					else
					{
						IMG[II] = (IMG[II]>=lgrsp)? IMG[II] : lgrsp;				
						if(IMG[II]<minIMout[0])
							minIMout[0] = IMG[II];
					}
				}				
				++iterate;		
			}
		}
	}	

	std::cout<<"done"<<std::endl;
  }
	
  return EXIT_SUCCESS;
}

float get_maximum_3D(float* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C)
{
	
   	float mx = A[(z1*R*C)+(r1*C)+c1];
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
			for(int k=z1; k<=z2; k++)
			{				
				if(A[(k*R*C)+(i*C)+j]>mx)
					mx = A[(k*R*C)+(i*C)+j];
			}
        }
    }
    return mx;
}

unsigned short get_maximum_3D(unsigned short* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C)
{
	
   	unsigned short mx = A[(z1*R*C)+(r1*C)+c1];
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
			for(int k=z1; k<=z2; k++)
			{				
				if(A[(k*R*C)+(i*C)+j]>mx)
					mx = A[(k*R*C)+(i*C)+j];
			}
        }
    }
    return mx;
}


void Detect_Local_MaximaPoints_3D(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1, unsigned short* bImg)
{  
    int min_r, min_c, max_r, max_c, min_z, max_z;    
       
    //start by getting local maxima points
    //if a point is a local maximam give it a local maximum ID
        
	//int IND = 0;
	int II = 0;
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {				
			for(int k=0; k<z; k++)
			{									
				min_r = (int) max(0.0,i-scale_xy);
				min_c = (int) max(0.0,j-scale_xy);
				min_z = (int) max(0.0,k-scale_z);
				max_r = (int)min((double)r-1,i+scale_xy);
				max_c = (int)min((double)c-1,j+scale_xy);                         
				max_z = (int)min((double)z-1,k+scale_z);                         
				float mx = get_maximum_3D(im_vals, min_r, max_r, min_c, max_c, min_z, max_z,r,c);
				II = (k*r*c)+(i*c)+j;
				if(im_vals[II] == mx)    
				{
					//IND = IND+1;
					//if(bImg[II] > 0)
						out1[II]=255;                 
					//else
					//	out1[II] = -1;
				}
				else
					out1[(k*r*c)+(i*c)+j]=0;

			}			
        }
    }  

}

//added by Yousef on 8/29/2009
//Estimate the min and max scales based on the local maxima points of the distance map
void estimateMinMaxScales(itk::SmartPointer<MyInputImageType> im, unsigned short* distIm, double* minScale, double* maxScale, int r, int c, int z)
{
	int min_r, min_c, max_r, max_c, min_z, max_z;    
	int II = 0;
	minScale[0] = 1000.0;
	maxScale[0] = 0.0;
	int cent_slice = (int) z/2;
	std::vector< std::vector<unsigned short> > scales;
	double mean = 0.0;
	double stdv = 0.0;
	int cnt = 0;
	//ofstream p;
	//p.open("checkme.txt");
	for(int i=1; i<r-1; i+=2)
    {
        for(int j=1; j<c-1; j+=2)
        {				
			//for(int k=1; k<z-1; k++)
			for(int k=cent_slice; k<=cent_slice; k++)
			{									
				min_r = (int) max(0.0,(double)i-2);
				min_c = (int) max(0.0,(double)j-2);
				min_z = (int) max(0.0,(double)k-2);
				max_r = (int)min((double)r-1,(double)i+2);
				max_c = (int)min((double)c-1,(double)j+2);                         
				max_z = (int)min((double)z-1,(double)k+2);                         
				unsigned short mx = get_maximum_3D(distIm, min_r, max_r, min_c, min_c, min_z, max_z,r,c);
				if(mx <= 1)
					continue; //background or edge point
				II = (k*r*c)+(i*c)+j;
				if(distIm[II] == mx)    
				{																																						
					//add the selected scale to the list of scales
					std::vector <unsigned short> lst;
					lst.push_back(mx);
					lst.push_back(i);
					lst.push_back(j);
					lst.push_back(k);
					//p<<j<<" "<<i<<" "<<k<<std::endl;
					scales.push_back(lst);
					//mean +=mx;
					cnt++;										
				}				
			}			
        }
    } 
	//p.close();
	/*mean /= cnt;
	for(int i=0; i<cnt; i++)
		stdv+= ((scales[i][0]-mean)*(scales[i][0]-mean));
	stdv = sqrt(stdv/cnt);		
	if (ceil(mean-stdv)-(mean-stdv)<=.5)
		minScale[0] = ceil((mean-stdv));///sqrt(2.0));
	else	
		minScale[0] = floor((mean-stdv));///sqrt(2.0));	 

	if (ceil(mean+stdv)-(mean+stdv)<=.5)
		maxScale[0] = ceil((mean+stdv));///sqrt(2.0));	 
	else
		maxScale[0] = floor((mean+stdv));///sqrt(2.0));	 */

	//get the median of the scales(distances)
	int medianS = computeMedian(scales, cnt);

	//then compute the Median absolute deviation
	std::vector< std::vector<unsigned short> > madList;
	for(int i=0; i<cnt; i++)
	{
		std::vector<unsigned short> tmp;
		tmp.push_back(abs(scales[i][0]-medianS));
		madList.push_back(tmp);
	}
	int MAD = computeMedian(madList, cnt);
	minScale[0] = medianS-MAD;
	maxScale[0] = medianS+MAD;
	

	//For each local maximum point,try to find the best LoG scale
	//To do that, suppose the distance at a given local maximum point is d, 
	//then compute the its LoG responses at scales from d/2 to d
	//Then, select the scale the gave us the maximum LoG response
	//Create an itk image to hold the sub image (tile) being processing
	int mnScl = 10000;
	int mxScl = 0;
	int cnt2 = 0;
	for(int ind=0; ind<cnt; ind++)
	{
		int mx = scales[ind][0];
		int i = scales[ind][1];
		int j = scales[ind][2];
		int k = scales[ind][3];
		if(mx<minScale[0] || mx>maxScale[0])
			continue;

		int smin = (int) std::ceil(mx/2.0);
		if(smin == mx)
			continue;
		if(smin == 1)
			smin++;
		cnt2++;
		min_r = (int) max(0.0,(double)i-1.5*mx);
		min_c = (int) max(0.0,(double)j-1.5*mx);
		min_z = (int) max(0.0,(double)k-1.5*mx);
		max_r = (int)min((double)r-1,(double)i+1.5*mx);
		max_c = (int)min((double)c-1,(double)j+1.5*mx);                         
		max_z = (int)min((double)z-1,(double)k+1.5*mx);                         
					
		int sub_r = i-min_r;
		int sub_c = j-min_c;
		int sub_z = k-min_z;
		int sz_r = (max_r-min_r+1);
		int sz_c = (max_c-min_c+1);
		int sz_z = (max_z-min_z+1);
		int ind_i = sub_z*sz_r*sz_c+sub_r*sz_c+sub_c;
																	
		MyInputImageType::Pointer im_Small = extract3DImageRegion(im, sz_c, sz_r, sz_z, min_c, min_r, min_z);															
		float* IMG = new float[sz_c*sz_r*sz_z];
		float max_resp = -100000.0;	
		int best_scale = 0.0;
		double sigma;
		for(int kk=smin; kk<=mx; kk++)
		{						
			sigma = kk;
			detect_seeds(im_Small, sz_r, sz_c, sz_z, sigma, IMG, 0);
			if(IMG[ind_i]>=max_resp)
			{
				max_resp = IMG[ind_i];
				best_scale = kk;
			}
		}
		mx = best_scale;
		
		if(mx<mnScl)
			mnScl = mx;
		if(mx>mxScl)
			mxScl = mx;

		delete [] IMG;
	}
	minScale[0] = mnScl;
	maxScale[0] = mxScl;
}

int distMap(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG)
{
  
  //  Types should be selected on the desired input and output pixel types.  
  typedef unsigned short             InputPixelType2;
  typedef float          OutputPixelType2;

  //  The input and output image types are instantiated using the pixel types.
  typedef itk::Image< InputPixelType2,  3 >   InputImageType2;
  typedef itk::Image< OutputPixelType2, 3 >   OutputImageType2;


  //  The filter type is now instantiated using both the input image and the
  //  output image types.
  //typedef itk::ApproximateSignedDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;    
  //typedef itk::DanielssonDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;  
  typedef itk::SignedMaurerDistanceMapImageFilter<InputImageType2, OutputImageType2>  DTFilter;
  DTFilter::Pointer dt_obj= DTFilter::New() ;
  //dt_obj->UseImageSpacingOn();

  typedef itk::CastImageFilter< MyInputImageType, InputImageType2> myCasterType;
  myCasterType::Pointer potCaster = myCasterType::New();
  potCaster->SetInput( im );
  potCaster->Update();

  dt_obj->SetInput(potCaster->GetOutput()) ;
  dt_obj->SetSquaredDistance( false );
  dt_obj->SetUseImageSpacing( true );
  dt_obj->SetInsideIsPositive( false );

  //dt_obj->SetInsideValue(0.0);
  //dt_obj->SetOutsideValue(255.0);
  try{
	 dt_obj->Update() ;
  }
  catch( itk::ExceptionObject & err ){
	std::cerr << "Error calculating distance transform: " << err << endl ;
    return -1;
  }
 
  //   Copy the resulting image into the input array
  long int i = 0;
  typedef itk::ImageRegionIteratorWithIndex< OutputImageType2 > IteratorType;
  IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());
  while ( i<r*c*z)
  {	  
	  double ds = iterate.Get();
	  if(ds<=0)
		  IMG[i] = 0;
	  else
		IMG[i] = (unsigned short) ds;
      ++i;
 	  ++iterate;
  }	
 
  return EXIT_SUCCESS;
}


int distMap_SliceBySlice(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG)
{
  
  //  Types should be selected on the desired input and output pixel types.  
  typedef unsigned short             InputPixelType2;
  typedef float          OutputPixelType2;

  //  The input and output image types are instantiated using the pixel types.  
  typedef itk::Image< OutputPixelType2, 2 >   OutputImageType2;
  long int k = 0;
  for(int i=0; i<z; i++)
  {
	  MyInputImageType2D::Pointer image2D = extract2DImageSlice(im, 2, i);
	  typedef itk::SignedMaurerDistanceMapImageFilter<MyInputImageType2D, OutputImageType2>  DTFilter;
	  DTFilter::Pointer dt_obj= DTFilter::New() ;
	  dt_obj->SetInput(image2D) ;
	  dt_obj->SetSquaredDistance( false );      
	  dt_obj->SetInsideIsPositive( false );
	  try{
		dt_obj->Update() ;
      }
      catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
		return -1;
	  }

	  //By Yousef: try to write out the output at the central slice
	  /*int cent_slice = (int) z/2;
	  if(i==cent_slice)
	  {
		  typedef itk::CastImageFilter< OutputImageType2, MyInputImageType2D> myCasterType;
		  myCasterType::Pointer potCaster = myCasterType::New();
		  potCaster->SetInput( dt_obj->GetOutput() );
		  typedef itk::ImageFileWriter< MyInputImageType2D > WriterType;
		  WriterType::Pointer writer = WriterType::New();
		  writer->SetFileName("dist_test.tif");
		  writer->SetInput( potCaster->GetOutput() );
		  writer->Update();		  
	  }*/


	  //   Copy the resulting image into the input array  
      typedef itk::ImageRegionIteratorWithIndex< OutputImageType2 > IteratorType;
      IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());
	  int j=0;	  

      while ( j<r*c)
      {	  
		  double ds = iterate.Get();
		  if(ds<=0)
			 IMG[k] = 0;
		  else
			 IMG[k] = (unsigned short) ds;
		  ++k;
		  ++j;
 		  ++iterate;
	  }	  
  } 
 
  return EXIT_SUCCESS;
}


MyInputImageType2D::Pointer extract2DImageSlice(itk::SmartPointer<MyInputImageType> im, int plane, int slice) 
{
    typedef itk::ExtractImageFilter< MyInputImageType, MyInputImageType2D > FilterType2D;
	FilterType2D::Pointer filter = FilterType2D::New();
    
    MyInputImageType::RegionType inputRegion = im->GetLargestPossibleRegion();
    
    MyInputImageType::SizeType size = inputRegion.GetSize();
    size[plane] = 0;
    
    MyInputImageType::IndexType start = inputRegion.GetIndex();
    const unsigned int sliceNumber = slice;
    start[plane] = sliceNumber;
    
    MyInputImageType::RegionType desiredRegion;
    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( start );
    
    filter->SetExtractionRegion( desiredRegion );
    
    filter->SetInput( im );

    MyInputImageType2D::Pointer img = filter->GetOutput();
    img->Update();

    return img;
}

MyInputImageType::Pointer extract3DImageRegion(itk::SmartPointer<MyInputImageType> im, int sz_x, int sz_y, int sz_z, int start_x, int start_y, int start_z)
{
    typedef itk::ExtractImageFilter< MyInputImageType, MyInputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();   
    
    MyInputImageType::SizeType size;
	size[0] = sz_x;
	size[1] = sz_y;
	size[2] = sz_z;
    
    MyInputImageType::IndexType start;
    start[0] = start_x;
	start[1] = start_y;
	start[2] = start_z;
    
    MyInputImageType::RegionType desiredRegion;
    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( start );
	
    
    filter->SetExtractionRegion( desiredRegion );
    
    filter->SetInput( im );

    MyInputImageType::Pointer img = filter->GetOutput();
    img->Update();

    return img;
}

int computeMedian(std::vector< std::vector<unsigned short> > scales, int cntr)
{
	if(cntr == 1)
		return scales[0][0];

	unsigned short* srtList = new unsigned short[cntr];
	for(int i=0; i<cntr-1; i++)
		srtList[i] = scales[i][0];

	for(int i=0; i<cntr-1; i++)
	{
		for(int j=i+1; j<cntr; j++)
		{
			if(srtList[j]<srtList[i])
			{
				unsigned short tmp = srtList[i];
				srtList[i] = srtList[j];
				srtList[j] = tmp;
			}
		}
	}

	int mdn = (int) cntr/2;
	mdn = srtList[mdn-1];

	return mdn;
}
