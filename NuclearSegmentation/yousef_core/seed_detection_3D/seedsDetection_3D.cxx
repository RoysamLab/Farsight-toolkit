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

#include "itkImage.h"
#include "itklaplacianrecursivegaussianimagefilternew.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkDanielssonDistanceMapImageFilter.h> 
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"


////
#include "itkImageFileWriter.h"

using namespace std;


typedef    unsigned short     MyInputPixelType;
typedef itk::Image< MyInputPixelType,  3 >   MyInputImageType;

int detect_seeds(itk::SmartPointer<MyInputImageType>, int , int , int, const double, float*, int);
float get_maximum_3D(float* A, int r1, int r2, int c1, int c2, int z1, int z2, int R, int C);
void Detect_Local_MaximaPoints_3D(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1, unsigned short* bImg);
int distMap(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG);

int Seeds_Detection_3D( float* IM, float* IM_out, unsigned short* IM_bin, int r, int c, int z, double sigma_min, double sigma_max, double scale_xy, double scale_z, int sampl_ratio, unsigned short* bImg, int UseDistMap)
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

	///////////////////////
	unsigned short* dImg = NULL;
	unsigned short* minSImg = NULL;
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
		distMap(im, r, c, z,dImg);
		std::cout<<"done"<<std::endl;
		iterator1.GoToBegin();
		//get the smallest acceptable scale at each point
		//Also, change the value in dImg into the largest acceptable scale
		/*minSImg = (unsigned short *) malloc(r*c*z*sizeof(unsigned short));
		for(int i=0; i<r*c*z; i++)
		{			
			int DD = (int) dImg[i];
			int DD2 = 2*DD;
			if(DD<=sigma_min)
				minSImg[i] = (unsigned short) sigma_min;
			else if(DD>=sigma_max)
				minSImg[i] = (unsigned short) sigma_max;
			else
				minSImg[i] = (unsigned short) DD;
			if(DD2<=sigma_min)
				dImg[i] = sigma_min;
			else if(DD2>=sigma_max)
				dImg[i] = sigma_max;
			else
				dImg[i] = DD2;			
		}*/

	}
	//////////////////////////

	for(int i=0; i<r*c*z; i++)
	{				
	    iterator1.Set((unsigned short)IM[i]);	
		++iterator1;
	}

	//Start from sigma_min to sigma sigma_max
	int minIMout = 10000;
	double conv = 0;	
	double sigma = sigma_min;
	//float *IMG_tmp = (float *) malloc(r*c*z*sizeof(float));
	std::cout<<"Processing scale "<<sigma<<"...";
	detect_seeds(im,r,c,z,sigma_min,IM_out,sampl_ratio);
	std::cout<<"done"<<std::endl;
	while(!conv)
	{
		
		sigma = sigma+1;
		
		if(sigma>sigma_max)
		{
			conv=1;
			break;
		}
		std::cout<<"Processing scale "<<sigma<<"...";
//		detect_seeds(im,r,c,z,sigma,IMG_tmp,sampl_ratio);
		detect_seeds(im,r,c,z,sigma,IM,sampl_ratio);
		if(UseDistMap == 1)
		{
			for(int i=0; i<r*c*z; i++)
			{						
				if(sigma<=dImg[i]*2)
					//IM_out[i] = (IM_out[i]>=IMG_tmp[i])? IM_out[i] : IMG_tmp[i];				
					IM_out[i] = (IM_out[i]>=IM[i])? IM_out[i] : IM[i];				
				if(IM_out[i]<minIMout)
					minIMout = IM_out[i];
			}
		}
		else
		{
			for(int i=0; i<r*c*z; i++)
			{			
				//IM_out[i] = (IM_out[i]>=IMG_tmp[i])? IM_out[i] : IMG_tmp[i];				
				IM_out[i] = (IM_out[i]>=IM[i])? IM_out[i] : IM[i];				
				if(IM_out[i]<minIMout)
					minIMout = IM_out[i];
			}
		}
		std::cout<<"done"<<std::endl;
	}

	//Just try this: Multiply the LoG response by a weight based on the distance from the background
	/*for(int i=0; i<r*c*z; i++)
	{
		IM_out[i] *= (dImg[i]/sigma_min);
	}*/
	//free(IMG_tmp);
	free(dImg);
   

	//Detect the seed points (which are also the local maxima points)	
	std::cout<<"Detecting Seeds"<<std::endl;
	Detect_Local_MaximaPoints_3D(IM_out, r, c, z, scale_xy, scale_z, IM_bin, bImg);	
	
	return minIMout;
}

int detect_seeds(itk::SmartPointer<MyInputImageType> im, int r, int c, int z,const double sigma, float* IMG, int sampl_ratio)
{
  
  //  Types should be selected on the desired input and output pixel types.
  //typedef    float     InputPixelType;
  typedef    float     OutputPixelType;


  //  The input and output image types are instantiated using the pixel types.
  //typedef itk::Image< InputPixelType,  3 >   InputImageType;
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