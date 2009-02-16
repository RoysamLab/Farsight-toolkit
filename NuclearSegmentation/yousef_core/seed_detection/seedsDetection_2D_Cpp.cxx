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
//#include "itkImageFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
//#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkSimpleFilterWatcher.h"

using namespace std;

//,.,.

typedef    float     InputPixelType;
typedef itk::Image< InputPixelType,  2 >   InputImageType;

int detect_seeds(itk::SmartPointer<InputImageType>, int , int , const double, float*);

double get_maximum(double** A, int r1, int r2, int c1, int c2);
void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, int* im_bin);
////////////////////////////////////

int detectSeeds2D( float* IM, float* IM_out, int* IM_bin, int r, int c, double sigma_min, double sigma_max, double scale)
{
    //Now, we are expecting the following inputs from IDL:
	//1-The imput image that need to be processed {argv[0]}
	//2-An empty image to hold the output of the multi-scale LoG {argv[1]}
	//3-An empty image to hold the binary image of the local maxima {argv[2]}
	//4-The image dimension r and c {argv[3] & {argv[4]}}
	//5-The filter scales: sigma_min {argv[5]} & sigma_max {argv[6]}
	  
	/*float* IM = (float *) argv[0];
	float* IM_out = (float *) argv[1];
	unsigned char* IM_bin = (unsigned char *) argv[2];
	int r = *(int*) argv[3];
	int c = *(int*) argv[4];
	double sigma_min = *(double*) argv[5];
	double sigma_max = *(double*) argv[6];
	double scale = *(double*) argv[7];*/
	
		
	//Create an itk image
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
    InputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
  
    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(int i=0; i<r*c; i++)
	{		
		iterator1.Set(IM[i]);
		++iterator1;	
	}
    	
	
	//Start from sigma_min to sigma sigma_max	
	double conv = 0;	
	double sigma = sigma_min;
	float *IMG_tmp = (float *) malloc(r*c*sizeof(float));	
	detect_seeds(im,r,c,sigma_min,IM_out);
	while(!conv)
	{		
		sigma = sigma+1;		
		if(sigma>sigma_max)
		{
			conv=1;
			break;
		}
		detect_seeds(im,r,c,sigma,IMG_tmp);

		for(int i=0; i<r*c; i++)
		{
			IM_out[i] = (IM_out[i]>=IMG_tmp[i])? IM_out[i] : IMG_tmp[i];	
		}
	}
	
    //Detect the seed points (which are also the local maxima points)
	Detect_Local_MaximaPoints(IM_out, r, c, scale, IM_bin);

	
	return 1;
}

int detect_seeds(itk::SmartPointer<InputImageType> im, int r, int c, const double sigma, float* IMG)
{  
  typedef    float     InputPixelType;
  typedef    float     OutputPixelType;
  typedef itk::Image< InputPixelType,  2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
 

  //Initialize the laplacian of gaussian filter
  typedef itk::LaplacianRecursiveGaussianImageFilter<InputImageType, OutputImageType >  FilterType;
  FilterType::Pointer laplacian = FilterType::New();
  
  //  Set the option for normalizing across scale space to true.
  laplacian->SetNormalizeAcrossScale( true );
  
  //set the input image
  laplacian->SetInput( im);
  
  //set the current sigma (scale)
  laplacian->SetSigma( sigma );
 
 
  //  Finally the pipeline is executed by invoking the Update() method.
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




  //Now, copy the resulting image into an array and return it  
  long int i = 0;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
  IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
  while ( i<r*c)
  {
    IMG[i] = iterate.Get();
    ++i;
	++iterate;
  }
  return EXIT_SUCCESS;

}

float get_maximum(float** A, int r1, int r2, int c1, int c2)
{
    float mx = A[r1][c1];
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
            if(A[i][j]>mx)
                mx = A[i][j];
        }
    }
    return mx;
}


void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, int* out1)
{  
    float** im;
    int min_r, min_c, max_r, max_c;    
   
    
    
    im = (float **) malloc(r*sizeof(float*));       
    for(int i=0; i<r; i++)
    {
        im[i] = (float *) malloc(c*sizeof(float));
        for(int j=0; j<c; j++)
        {
            im[i][j] = im_vals[(i*c)+j];            
        }
    }    
        

    //start by getting local maxima points
    //if a point is a local maximam give it a local maximum ID
        
	int IND = 0;
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {				
            min_r = (int) max(0.0,i-scale);
            min_c = (int) max(0.0,j-scale);
            max_r = (int)min((double)r-1,i+scale);
            max_c = (int)min((double)c-1,j+scale);                         
            float mx = get_maximum(im, min_r, max_r, min_c, max_c);
            if(im[i][j] == mx)  
			{
				IND = IND+1;
                out1[(i*c)+j]=255;//(unsigned char)IND; //same as sub2ind(size(im),i,j);                                         
			}
        }
    }  
    	        
}

