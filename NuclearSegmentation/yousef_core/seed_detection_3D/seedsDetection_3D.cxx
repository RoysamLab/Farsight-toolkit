/* 
* Copyright 2009 Rensselaer Polytechnic Institute
* This program is free software; you can redistribute it and/or modify 
* it under the terms of the GNU General Public License as published by 
* the Free Software Foundation; either version 2 of the License, or 
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but 
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
* or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
* for more details.
* 
* You should have received a copy of the GNU General Public License along 
* with this program; if not, write to the Free Software Foundation, Inc., 
* 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

//seedsDetection_2D.cxx

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <limits.h>
#include <time.h>

//_OPENMP is automatically defined by any OpenMP compiler that is building with OpenMP support (-fopenmp for GCC or /openmp for MSVC)
#ifdef _OPENMP
#include "omp.h"
#endif

//OpenCL support
#ifdef OPENCL
#include "CL\cl.h"
#include "CL\cl_ext.h"
#ifndef MSTRINGIFY
#define MSTRINGIFY(A) #A
char* stringifiedKernel =
#include "TestKernel.cl"
char* LocalMaximaKernel =
#include "LocalMaximaKernel.cl"
#endif
#endif


#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkDanielssonDistanceMapImageFilter.h> 
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkMultiThreader.h"

//added by Yousef on 8/26/2009
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"

typedef    unsigned short     MyInputPixelType;
typedef itk::Image< MyInputPixelType,  3 >   MyInputImageType;
typedef itk::Image< MyInputPixelType,  2 >   MyInputImageType2D;

int detect_seeds(itk::SmartPointer<MyInputImageType>, int , int , int, const double, float*, int);
//int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, int r, int c, int z,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, int unsigned short* dImg, int* minIMout, int UseDistMap);
int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, size_t r, size_t c, size_t z, int rmin, int rmax, int cmin, int cmax, int zmin, int zmax,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, int unsigned short* dImg, int* minIMout, int UseDistMap);
float get_maximum_3D(float* A, int r1, int r2, int c1, int c2, int z1, int z2, int R, int C);
unsigned short get_maximum_3D(unsigned short* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C);
void Detect_Local_MaximaPoints_3D(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1, unsigned short* bImg);
int distMap(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG);
int distMap_SliceBySlice(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG);
MyInputImageType2D::Pointer extract2DImageSlice(itk::SmartPointer<MyInputImageType> im, int plane, int slice);
MyInputImageType::Pointer extract3DImageRegion(itk::SmartPointer<MyInputImageType> im, int sz_x, int sz_y, int sz_z, int start_x, int start_y, int start_z);
void estimateMinMaxScales(itk::SmartPointer<MyInputImageType> im, unsigned short* distIm, double* minScale, double* maxScale, int r, int c, int z);
int computeMedian(std::vector< std::vector<unsigned short> > scales, int cntr);
void estimateMinMaxScalesV2(itk::SmartPointer<MyInputImageType> im, unsigned short* distIm, double* minScale, double* maxScale, int r, int c, int z);
int computeWeightedMedian(std::vector< std::vector<float> > scales, int cntr);
void queryOpenCLProperties(float* IM, int r, int c, int z);
std::string fileToString(std::string fileName);
void seed_pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data);
void Detect_Local_MaximaPoints_3D_ocl(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1);

//external functions
extern "C" void Detect_Local_MaximaPoints_3D_CUDA(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1);


int Seeds_Detection_3D( float* IM, float** IM_out, unsigned short** IM_bin, size_t r, size_t c, size_t z, double *sigma_min_in, double *sigma_max_in, double *scale_xy_in, double *scale_z_in, int sampl_ratio, unsigned short* bImg, int UseDistMap, int* minIMout, bool paramEstimation)
{	
	// 	std::cout << std::endl << 
	//queryOpenCLProperties(IM, r, c, z); //odd that the OpenCL code will crash the program if it is not part of a function and manually inlined here instead...

	//get this inputs
	double sigma_min = sigma_min_in[0];
	double sigma_max = sigma_max_in[0];
	double scale_xy = scale_xy_in[0];
	double scale_z = scale_z_in[0];

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

	try{
		im->SetRegions( region );
		im->SetSpacing(spacing);
		im->Allocate();
		im->FillBuffer(0);
		im->Update();
	}
	catch( itk::ExceptionObject & excep ){
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< MyInputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());

	unsigned short* dImg = NULL;	
	int max_dist = 1;
	iterator1.GoToBegin();
	if(UseDistMap == 1)
	{
		unsigned long i=0;
		while( !iterator1.IsAtEnd() )
		{
			if(bImg[i]>0)
				iterator1.Set(0.0);//IM[i]);
			else
				iterator1.Set(255.0);
			++iterator1; ++i;

		}
		//By Yousef on 09/08/2009
		//Just for testing purposes: Write out the distance map into an image file
		//typedef itk::ImageFileWriter< MyInputImageType > WriterType;
		//WriterType::Pointer writer = WriterType::New();
		//writer->SetFileName("bin_test.tif");
		//writer->SetInput( im );
		//writer->Update();
		//////////////////////////////////////////////////////////////////////////
		std::cout<<"Computing distance transform...";
		dImg = (unsigned short *) malloc(((unsigned long)r)*((unsigned long)c)*((unsigned long)z)*sizeof(unsigned short));
		if(!dImg)
		{
			std::cerr<<"Failed to allocate memory for the distance image, GOING TO DI"<<std::endl;
 			return 0;
		}
		//max_dist = distMap(im, r, c, z,dImg);
		max_dist = distMap_SliceBySlice(im, r, c, z,dImg);
		std::cout<<"done"<<std::endl;
	}

	float multp = 1.0;
	iterator1.GoToBegin();
	unsigned long i=0;
	while( !iterator1.IsAtEnd() )
	{

		if(UseDistMap == 1)
			multp = 1+((float) dImg[i]/(2*max_dist));
		if(bImg[i] > 0)
			iterator1.Set((unsigned short)IM[i]/multp);			
		else
			iterator1.Set(255);
		++iterator1; ++i;
	}

	std::cout << "About to enter Estimating parameters" << std::endl;
	//By Yousef (8/29/2009)
	//Estimate the segmentation parameters
	if(UseDistMap == 1 && paramEstimation)
	{
		std::cout<<"Estimating parameters..."<<std::endl;
		clock_t start_time_est_params = clock();
		estimateMinMaxScalesV2(im, dImg, &sigma_min, &sigma_max, r, c, z);
		std::cout << "Estimating parameters took " << (clock() - start_time_est_params)/(float)CLOCKS_PER_SEC << " seconds" << std::endl;
		//By Isaac (1/22/2010)
		if( sigma_min < 1 ) sigma_min = 1;

		scale_xy = sigma_min;
		if(scale_xy<3)
			scale_xy = 3; //just avoid very small search boxes
		scale_z = ceil(scale_xy / sampl_ratio);
		std::cout<<"    Minimum scale = "<<sigma_min<<std::endl;
		std::cout<<"    Maximum scale = "<<sigma_max<<std::endl;
		std::cout<<"    Clustering Resolution = "<<scale_xy<<std::endl;
		//write out the parameters
		sigma_min_in[0] = sigma_min;
		sigma_max_in[0] = sigma_max;
		scale_xy_in[0] =  scale_xy;
		scale_z_in[0] = scale_z;
	}

	//By Yousef (8/28/2009)
	//In some situations the image is very larg and we cannot allocate memory for the LoG filter (20xthe size of the image in bytes)
	//In such cases, we can divide the image into small tiles, process them independently, and them combine the results
	//see if we have enought memory for the LoG step
	//approximately, we need (20~21)ximage size in bytes
	//try to allocate memory for an unsigned char* of the 23ximage size
	int block_divisor = 1; //Assumes you will not overflow RAM (according to Yousef above, you need 24 * image size of RAM, otherwise performance suffers) -- this appears to be incorrect for large images -Ho

	int blk = 1;
	int cntr = 0;

	for(int i=0; i<r; i+=r/block_divisor)
		for(int j=0; j<c; j+=c/block_divisor)
			cntr++;

	int min_x, min_y, max_x, max_y;


#ifdef _OPENMP
	omp_set_nested(1);	//This turns on nesting so the omp parallel inside the multiScaleLog can run in parallel. Without this, OpenMP only generates threads for the outermost omp parallel for construct. 
						//It may make sense to turn this off if it causes too much thread contention. THINK CAREFULLY ABOUT IT.
#endif
	clock_t start_time_multiscale_log = clock();

	#pragma omp parallel for private(min_x, min_y, max_x, max_y)
	for(int i=0; i<r; i+=r/block_divisor)
	{
		for(int j=0; j<c; j+=c/block_divisor)
		{
			std::cout<<"LoG block "<<blk++<<" of "<<cntr<<std::endl;
			min_x = j; 
			max_x = 40+(int)j+c/block_divisor; //40 is the size of the overlapping between the two tiles along x
			min_y = i;
			max_y = 40+(int)i+r/block_divisor; //40 is the size of the overlapping between the two tiles along y
			if(max_x >= c)
				max_x = c-1;
			if(max_y >= r)
				max_y = r-1;
			
			if( block_divisor == 1 )
			{
				multiScaleLoG(im, r, c, z, min_y, max_y, min_x, max_x, 0, z-1, sigma_min, sigma_max, IM, sampl_ratio, dImg, minIMout, UseDistMap);
			}
			else
			{

				//Create an itk image to hold the sub image (tile) being processing
				MyInputImageType::Pointer im_Small = extract3DImageRegion(im, max_x-min_x+1, max_y-min_y+1, z, min_x, min_y, 0);

				//By Yousef (8/27/2009): multi-scale LoG is done in one function now
				multiScaleLoG(im_Small, r, c, z, min_y, max_y, min_x, max_x, 0, z-1, sigma_min, sigma_max, IM, sampl_ratio, dImg, minIMout, UseDistMap);
				//
			}
		}
	}

#ifdef _OPENMP
	omp_set_nested(0);
#endif

	std::cout << "Multiscale Log took " << (clock() - start_time_multiscale_log)/(float)CLOCKS_PER_SEC << " seconds" << std::endl;

	free(dImg);

	*IM_out = new float[(((unsigned long)r)*((unsigned long)c)*((unsigned long)z))];
	if(!*IM_out)
	{
		std::cerr<<"could not allocate memory for the LoG response image"<<std::endl;
		return 0;
	}
	for(int i=0; i<((unsigned long)r)*((unsigned long)c)*((unsigned long)z); ++i)
	{
		(*IM_out)[i] = IM[i];
	}
	//Detect the seed points (which are also the local maxima points)	
	std::cout<<"Detecting Seeds"<<std::endl;
	*IM_bin = new unsigned short[(((unsigned long)r)*((unsigned long)c)*((unsigned long)z))];
	if(!*IM_bin)
	{
		std::cerr<<"could not allocate memory for the LoG response calc"<<std::endl;
		return 0;
	}
	//std::cout << "about to call Detect_Local_MaximaPoints_3D" << std::endl;
	clock_t start_time_local_maxima = clock();
#ifdef OPENCL
	Detect_Local_MaximaPoints_3D_ocl(IM_out[0], r, c, z, scale_xy, scale_z, IM_bin[0]);
#elif CUDA
	Detect_Local_MaximaPoints_3D_CUDA(IM_out[0], r, c, z, scale_xy, scale_z, IM_bin[0]);
#else
	Detect_Local_MaximaPoints_3D(IM_out[0], r, c, z, scale_xy, scale_z, IM_bin[0], bImg);
#endif // OPENCL
	std::cout << "Local maxima point detection took " << (clock() - start_time_local_maxima)/(float)CLOCKS_PER_SEC << " seconds" << std::endl;

	std::cout << "done detecting seeds" << std::endl;

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
	typedef itk::LaplacianRecursiveGaussianImageFilter<MyInputImageType, OutputImageType >  FilterType;
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
	unsigned long i = 0;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType;
	IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
	iterate.GoToBegin();
	while ( !iterate.IsAtEnd() )
	{
		IMG[i] = /*sigma*sigma*/iterate.Get()/sqrt(sigma);
		++i;
		++iterate;
	}

	return EXIT_SUCCESS;
}

int multiScaleLoG(itk::SmartPointer<MyInputImageType> im, size_t r, size_t c, size_t z, int rmin, int rmax, int cmin, int cmax, int zmin, int zmax,const double sigma_min, double sigma_max, float* IMG, int sampl_ratio, unsigned short* dImg, int* minIMout, int UseDistMap)
{

    std::cout << "UseDistMap: " << UseDistMap << std::endl;
	//  Types should be selected on the desired input and output pixel types.
	typedef    float     OutputPixelType;
	//  The input and output image types are instantiated using the pixel types.
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	bool failed = false;

    //Initialize the max_log_response array
    for(size_t k1=zmin; k1<=zmax; k1++)
    {
        for(size_t i1=rmin; i1<=rmax; i1++)
        {
            for(size_t j1=cmin; j1<=cmax; j1++)
            {
                size_t image_index = k1 * r * c + i1 * c + j1;
                IMG[image_index] = -std::numeric_limits<float>::max();  //Initialize to the negative of the maximum value
            }
        }
    }
    
	#ifdef _OPENMP
		int num_procs = omp_get_num_procs();
	#else
		int num_procs = 1;
	#endif

	int num_scales = sigma_max - sigma_min + 1;
	int num_openmp_threads = std::min(num_scales, num_procs);
	int num_itk_threads = std::ceil((float)num_procs / num_openmp_threads);

	std::cout << "Thread settings for multiScaleLog()" << std::endl;
	std::cout << "Number of scales: "			<< num_scales << std::endl;
	std::cout << "Number of OpenMP threads: "	<< num_openmp_threads << std::endl; 
	std::cout << "Number of ITK threads: "		<< num_itk_threads << std::endl;

	#pragma omp parallel for firstprivate(im)
	for(int i = sigma_max - sigma_min; i >= 0; i--)
	{
		int sigma = sigma_max - i;
		#pragma omp critical (ProcessingScaleCout)
		{
			std::cout << "Processing scale " << sigma << std::endl;
		}
		//  The filter type is now instantiated using both the input image and the 
		//  output image types.
		typedef itk::LaplacianRecursiveGaussianImageFilter<MyInputImageType, OutputImageType >  FilterType;
		FilterType::Pointer laplacian = FilterType::New();
		//  The option for normalizing across scale space can also be selected in this filter.
		laplacian->SetNormalizeAcrossScale( true );

		//  The input image can be obtained from the output of another
		//  filter. Here the image comming from the calling function is used as the source
		laplacian->SetInput(im);

		//  It is now time to select the $\sigma$ of the Gaussian used to smooth the
		//  data.  Note that $\sigma$ must be passed to both filters and that sigma
		//  is considered to be in millimeters. That is, at the moment of applying
		//  the smoothing process, the filter will take into account the spacing
		//  values defined in the image.
		//
		laplacian->SetSigma(sigma);

		//Set the number of itk threads based on the number of OpenMP threads
		laplacian->SetNumberOfThreads(num_itk_threads);
		
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
			failed = true;
		} 

		//   Copy the resulting image into the input array
		typedef itk::ImageRegionConstIterator< OutputImageType > IteratorType;
		IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
		iterate.GoToBegin();
        
        if (UseDistMap)
        {
            for(size_t k1=zmin; k1<=zmax; k1++)
            {
                for(size_t i1=rmin; i1<=rmax; i1++)
                {
                    for(size_t j1=cmin; j1<=cmax; j1++)
                    {                        
                        size_t image_index = k1 * r * c + i1 * c + j1;
                        float log_response = iterate.Get();
                        
                        if (sigma == sigma_min || sigma <= dImg[image_index] / 100)
                        {
                            #pragma omp critical
                            {
                                IMG[image_index] = std::max(IMG[image_index], log_response);
                                *minIMout = std::min<float>(*minIMout, IMG[image_index]);
                            }
                            ++iterate;
                        }
                    }
                }
            }
        }
        else    //not UseDistMap
        {
            for(size_t k1=zmin; k1<=zmax; k1++)
            {
                for(size_t i1=rmin; i1<=rmax; i1++)
                {
                    for(size_t j1=cmin; j1<=cmax; j1++)
                    {                
                        size_t image_index = k1 * r * c + i1 * c + j1;
                        float log_response = iterate.Get();
                        
                        #pragma omp critical
                        {
                            IMG[image_index] = std::max(IMG[image_index], log_response);
                            *minIMout = std::min<float>(*minIMout, IMG[image_index]);
                        }
                        
                        ++iterate;
                    }
                }
            }
        }

		std::cout<<"Scale " << sigma << " done"<<std::endl;
	}
	if (failed)
		return EXIT_FAILURE;
	else
		return EXIT_SUCCESS;
}

float get_maximum_3D(float* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C)
{
	unsigned long II = ((unsigned long)z1)*((unsigned long)R)*((unsigned long)C)+
			   ((unsigned long)r1)*((unsigned long)C)+((unsigned long)c1);
	float mx = A[II];

	for(int i=r1; i<=r2; i++)
	{
		for(int j=c1; j<=c2; j++)
		{
			for(int k=z1; k<=z2; k++)
			{
				unsigned long III = ((unsigned long)k)*((unsigned long)R)*((unsigned long)C)
						   +((unsigned long)i)*((unsigned long)C)+((unsigned long)j);
				if(A[III]>mx)
					mx = A[III];
			}
		}
	}
	return mx;
}

unsigned short get_maximum_3D(unsigned short* A, int r1, int r2, int c1, int c2, int z1, int z2,int R, int C)
{
	unsigned long II = ((unsigned long)z1)*((unsigned long)R)*((unsigned long)C)+
			   ((unsigned long)r1)*((unsigned long)C)+((unsigned long)c1);
	unsigned short mx = A[II];
	for(int i=r1; i<=r2; i++)
	{
		for(int j=c1; j<=c2; j++)
		{
			for(int k=z1; k<=z2; k++)
			{
				unsigned long III = ((unsigned long)k)*((unsigned long)R)*((unsigned long)C)+
						    ((unsigned long)i)*((unsigned long)C)+((unsigned long)j);
				if(A[III]>mx)
					mx = A[III];
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
	unsigned long II = 0;
	//int itr = 0;
	//std::cout << "In Detect_Local_MaximaPoints_3D, about to plunge in the loop" << std::endl;

	#pragma omp parallel for private(II, min_r, min_c, min_z, max_r, max_c, max_z)
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{				
			for(int k=0; k<z; k++)
			{									
				//calculate bounds
				min_r = (int) std::max(0.0,i-scale_xy);
				min_c = (int) std::max(0.0,j-scale_xy);
				min_z = (int) std::max(0.0,k-scale_z);
				max_r = (int) std::min((double)r-1,i+scale_xy);
				max_c = (int) std::min((double)c-1,j+scale_xy);                         
				max_z = (int) std::min((double)z-1,k+scale_z);                         

				//get the intensity maximum of the bounded im_vals
				float mx = get_maximum_3D(im_vals, min_r, max_r, min_c, max_c, min_z, max_z,r,c);

				//if the current pixel is at the maximum intensity, set it to 255 in out1 (seedImagePtr), else set it to 0
				II = ((unsigned long)k*(unsigned long)r*(unsigned long)c)
					+((unsigned long)i*(unsigned long)c)+(unsigned long)j;
				if(im_vals[II] == mx)    
					out1[II]=255;
				else
					out1[II]=0;
			}			
		}
	}  
	//std::cout << std::endl << "made it out of the loop" << std::endl;
}

void Detect_Local_MaximaPoints_3D_ocl(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1)
{
#ifdef OPENCL

	// START OPENCL BOILERPLATE ----------------------------------------------------------------------------------------------------------------
	cl_platform_id platforms[10];
	cl_uint num_platforms;
	cl_device_id device[10];
	cl_uint num_devices;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel kernel;
	cl_int errorcode;

	// Platform
	clGetPlatformIDs(10, platforms, &num_platforms); //Get OpenCL Platforms
	clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 10, device, &num_devices); //Get the first platform and get a list of devices
	context = clCreateContext(0, 1, device, &seed_pfn_notify, NULL, NULL); //Create context from first device
	queue = clCreateCommandQueue(context, device[0], 0, NULL); //Create a command queue for the first device

	//cout << endl << LocalMaximaKernel << endl << endl; //print out strinfified kernel

	program = clCreateProgramWithSource(context, 1, (const char **) &LocalMaximaKernel, 0, &errorcode); //Read in kernel and create a program

	if (errorcode != CL_SUCCESS)
		cout << "clCreateProgramWithSource Error code: " << errorcode << endl;

	errorcode = clBuildProgram(program, 0, 0, 0, 0, 0); //Build the program

	if (errorcode != CL_SUCCESS) //If there was a build error, print out build_info
	{
		cout << "clBuildProgram Error Code: " << errorcode << endl;
		char build_log[1024*1024];
		errorcode = clGetProgramBuildInfo(program, device[0], CL_PROGRAM_BUILD_LOG, sizeof(build_log), build_log, NULL); //Get the build log
		if (errorcode == CL_SUCCESS)
			cout << "Build Log:" << endl << build_log << endl;
		else
			cout << "clGetProgramBuildInfo Error Code: " << errorcode << endl;
	}

	kernel = clCreateKernel(program, "LocalMaximaKernel", &errorcode); //Create the kernel from the source code

	if (errorcode != CL_SUCCESS)
		cout << "clCreateKernel Error code: " << errorcode << endl;

	//END OPENCL BOILERPLATE ---------------------------------------------------------------------------------------------------------------------

	size_t cnDimension = r * c * z; //array size

	cout << "Allocating " << (sizeof(*im_vals) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for im_vals" << endl;
	cout << "Allocating " << (sizeof(*out1) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for out1" << endl;

	//Allocate device memory
	cl_mem device_mem_im_vals = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * cnDimension, NULL, NULL);
	cl_mem device_mem_out1 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ushort) * cnDimension, NULL, NULL);

	if (device_mem_im_vals == NULL || device_mem_out1 == NULL)
		cout << "Failed to allocate buffer memory on GPU" << endl; 

	//Write memory from host to device
	clEnqueueWriteBuffer(queue, device_mem_im_vals, CL_TRUE, 0, sizeof(cl_float) * cnDimension, im_vals, NULL, NULL, NULL);

	//Set kernel arguments
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &device_mem_im_vals);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &device_mem_out1);
	clSetKernelArg(kernel, 2, sizeof(int), &r);
	clSetKernelArg(kernel, 3, sizeof(int), &c);
	clSetKernelArg(kernel, 4, sizeof(int), &z);
	clSetKernelArg(kernel, 5, sizeof(double), &scale_xy);
	clSetKernelArg(kernel, 6, sizeof(double), &scale_z);	

	//Execute the kernel
	clEnqueueNDRangeKernel(queue, kernel, 1, 0, (const size_t *) &cnDimension, 0, 0, 0, 0);

	//Read the output from the device back into host memory
	clEnqueueReadBuffer(queue, device_mem_out1, CL_TRUE, 0, sizeof(cl_ushort) * cnDimension, out1, NULL, NULL, NULL);

	//Block till all commands are complete
	clFinish(queue);

	cout << endl;

	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clReleaseMemObject(device_mem_im_vals);
	clReleaseMemObject(device_mem_out1);
	clReleaseCommandQueue(queue);
	clReleaseContext(context);

	cout << endl;
#endif
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
	//	double mean = 0.0;
	//	double stdv = 0.0;
	int cnt = 0;
	std::ofstream p;
	//	int max_dist = 0;
	//p.open("checkme.txt");
	for(int i=1; i<r-1; i++)
	{
		for(int j=1; j<c-1; j++)
		{				
			//for(int k=1; k<z-1; k+=2)
			for(int k=cent_slice; k<=cent_slice; k++)
			{									
				min_r = (int) std::max(0.0,(double)i-2);
				min_c = (int) std::max(0.0,(double)j-2);
				min_z = (int) std::max(0.0,(double)k);
				max_r = (int) std::min((double)r-1,(double)i+2);
				max_c = (int) std::min((double)c-1,(double)j+2);                         
				max_z = (int) std::min((double)z-1,(double)k);                         
				unsigned short mx = get_maximum_3D(distIm, min_r, max_r, min_c, max_c, min_z, max_z,r,c);

				if(mx <= 100)
					continue; //background or edge point
				II = (k*r*c)+(i*c)+j;
				if(distIm[II] == mx)    
				{		
					//since we have scaled by 100 earlier, scale back to the original value
					//also, we want to dived by squre root of 2 or approximately 1.4
					mx = mx/140;							
					//add the selected scale to the list of scales
					std::vector <unsigned short> lst;
					lst.push_back(mx);
					lst.push_back(i);
					lst.push_back(j);
					lst.push_back(k);
					p<<j<<" "<<i<<" "<<k<<" "<<mx<<std::endl;
					scales.push_back(lst);
					//mean +=mx;
					cnt++;										
				}				
			}			
		}
	} 
	//p.close();

	//get the median of the scales(distances)
	int medianS = computeMedian(scales, cnt);
	//ofstream p2;
	//p2.open("checkme2.txt");
	//p2<<"med = "<<medianS<<std::endl;

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

	//p2<<"mad = "<<MAD<<std::endl;
	//p2<<"med-mad = "<<minScale[0]<<std::endl;
	//p2<<"med+mad = "<<maxScale[0]<<std::endl;

	std::ofstream p3;
	p3.open("checkme3.txt");	
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
		//if(mx<minScale[0] || mx>maxScale[0])
		//	continue;

		int smin = (int) std::ceil(mx/2.0);
		if(smin == mx)
			continue;
		if(smin == 1)
			smin++;
		cnt2++;
		min_r = (int) std::max(0.0,(double)i-mx);
		min_c = (int) std::max(0.0,(double)j-mx);
		min_z = (int) std::max(0.0,(double)k-mx);
		max_r = (int) std::min((double)r-1,(double)i+mx);
		max_c = (int) std::min((double)c-1,(double)j+mx);                         
		max_z = (int) std::min((double)z-1,(double)k+mx);                         

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
			//Method 1:
			//Get the scale at which the LoG response at our point of interest is maximum
			if(IMG[ind_i]>=max_resp)
			{
				max_resp = IMG[ind_i];								
				best_scale = kk;				
			}
			//Method 2:
			//We need the scale at which the maximum response in the small image region surrouding our point of interest 
			//is maximum over scales
			/*float mx2 = get_maximum_3D(IMG, 0, sz_r-1, 0, sz_c-1, 0, sz_z-1,sz_r,sz_c);
			if(mx2>=max_resp)
			{
			max_resp = mx2;								
			best_scale = kk;				
			}*/
		}
		p3<<j<<" "<<i<<" "<<k<<" "<<mx<<" "<<best_scale<<" "<<max_resp<<std::endl;
		mx = best_scale;

		if(mx<mnScl)
			mnScl = mx;
		if(mx>mxScl)
			mxScl = mx;

		delete [] IMG;
	}
	//I assume at least 4 scales must be used (will be relaxed later)
	if(mxScl<mnScl+3)
		mxScl = mnScl+3;
	minScale[0] = mnScl;
	maxScale[0] = mxScl;	
	//p2<<"min_scale="<<mnScl<<std::endl;
	//p2<<"max_scale="<<mxScl<<std::endl;
	//p2.close();
	p3.close();
}


//added by Yousef on 9/17/2009
//Estimate the min and max scales based on the local maxima points of the distance map
void estimateMinMaxScalesV2(itk::SmartPointer<MyInputImageType> im, unsigned short* distIm, double* minScale, double* maxScale, int r, int c, int z)
{
	int min_r, min_c, max_r, max_c, min_z, max_z;    
	int II = 0;
	minScale[0] = 1000.0;
	maxScale[0] = 0.0;
	int cent_slice = (int) z/2;
	std::vector< std::vector<unsigned short> > scales;
	//double mean = 0.0;
	//double stdv = 0.0;
	int cnt = 0;
	std::ofstream p;
	//int max_dist = 0;
	//p.open("checkme.txt");
	for(int i=1; i<r-1; i++)
	{
		for(int j=1; j<c-1; j++)
		{				
			//for(int k=1; k<z-1; k+=2)
			for(int k=cent_slice; k<=cent_slice; k++)
			{									
				min_r = (int) std::max(0.0,(double)i-2);
				min_c = (int) std::max(0.0,(double)j-2);
				min_z = (int) std::max(0.0,(double)k);
				max_r = (int) std::min((double)r-1,(double)i+2);
				max_c = (int) std::min((double)c-1,(double)j+2);                         
				max_z = (int) std::min((double)z-1,(double)k);                         
				unsigned short mx = get_maximum_3D(distIm, min_r, max_r, min_c, max_c, min_z, max_z,r,c);

				if(mx <= 100)
					continue; //background or edge point
				II = (k*r*c)+(i*c)+j;
				if(distIm[II] == mx)    
				{		
					//since we have scaled by 100 earlier, scale back to the original value
					//also, we want to dived by squre root of 2 or approximately 1.4
					mx = mx/140;							
					//add the selected scale to the list of scales
					std::vector <unsigned short> lst;

					lst.push_back(mx);
					lst.push_back(i);
					lst.push_back(j);
					lst.push_back(k);
					p<<j<<" "<<i<<" "<<k<<" "<<mx<<std::endl;
					scales.push_back(lst);

					//mean +=mx;
					cnt++;
				}				
			}			
		}
	} 
	//p.close();

	//get the median of the scales(distances)
	int medianS = computeMedian(scales, cnt);
	//ofstream p2;
	//p2.open("checkme2.txt");
	//p2<<"med = "<<medianS<<std::endl;				

	//ofstream p3;
	//p3.open("checkme3.txt");	
	//For each local maximum point,try to find the best LoG scale
	//To do that, suppose the distance at a given local maximum point is d, 
	//then compute the its LoG responses at scales from d/2 to d
	//Then, select the scale the gave us the maximum LoG response
	int mnScl = 10000;
	int mxScl = 0;
	int cnt2 = 0;
	std::vector<std::vector<float> > smallScales;
	std::vector<std::vector<float> > largeScales;
	int numSmall = 0;
	int numLarge = 0;

	for(int ind=0; ind<cnt; ind++)
	{
		int mx = scales[ind][0];
		int i = scales[ind][1];
		int j = scales[ind][2];
		int k = scales[ind][3];		

		int smin = (int) std::ceil(mx/2.0);
		if(smin == mx)
			continue;
		if(smin == 1)
			smin++;
		cnt2++;
		min_r = (int) std::max(0.0,(double)i-mx);
		min_c = (int) std::max(0.0,(double)j-mx);
		min_z = (int) std::max(0.0,(double)k-mx);
		max_r = (int) std::min((double)r-1,(double)i+mx);
		max_c = (int) std::min((double)c-1,(double)j+mx);                         
		max_z = (int) std::min((double)z-1,(double)k+mx);                         

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
			//Method 1:
			//Get the scale at which the LoG response at our point of interest is maximum
			if(IMG[ind_i]>=max_resp)
			{
				max_resp = IMG[ind_i];								
				best_scale = kk;				
			}
			//Method 2:
			//We need the scale at which the maximum response in the small image region surrouding our point of interest 
			//is maximum over scales
			/*float mx2 = get_maximum_3D(IMG, 0, sz_r-1, 0, sz_c-1, 0, sz_z-1,sz_r,sz_c);
			if(mx2>=max_resp)
			{
			max_resp = mx2;								
			best_scale = kk;				
			}*/
		}
		std::vector<float> pp;
		pp.push_back(best_scale);
		pp.push_back(max_resp);

		if(mx<=medianS)
		{	
				numSmall++;		
				smallScales.push_back(pp);
		}
		else
		{
			numLarge++;		
			largeScales.push_back(pp);
		}

		//p3<<j<<" "<<i<<" "<<k<<" "<<mx<<" "<<best_scale<<" "<<max_resp<<std::endl;
		/*mx = best_scale;	
		if(mx<mnScl)
		mnScl = mx;
		if(mx>mxScl)
		mxScl = mx;*/

		delete [] IMG;
	}
	//p3.close();
	//set the min and max scales to the LoG-weighted medians of the small and large scale sets
	mnScl =  computeWeightedMedian(smallScales, numSmall);
	mxScl =  computeWeightedMedian(largeScales, numLarge);

	scales.clear();
	smallScales.clear();
	largeScales.clear();

	//I assume at least 4 scales must be used (will be relaxed later)
	if(mxScl<mnScl+3)
		mxScl = mnScl+3;
	minScale[0] = mnScl;
	maxScale[0] = mxScl;	
	//p2<<"min_scale="<<mnScl<<std::endl;
	//p2<<"max_scale="<<mxScl<<std::endl;
	//p2.close();	
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
		std::cerr << "Error calculating distance transform: " << err << std::endl ;
		return -1;
	}

	//   Copy the resulting image into the input array
	unsigned long i = 0;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType2 > IteratorType;
	IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());
	iterate.GoToBegin();

	unsigned long max_dist = 0;
	while (!iterate.IsAtEnd())
	{	  
		double ds = iterate.Get();
		if(ds<=0)
			IMG[i] = 0;
		else
			IMG[i] = (unsigned short) ds;

		if(IMG[i]>max_dist)
			max_dist = IMG[i];

		++i; ++iterate;
	}	
	return max_dist;
}


int distMap_SliceBySlice(itk::SmartPointer<MyInputImageType> im, int r, int c, int z, unsigned short* IMG)
{

	//  Types should be selected on the desired input and output pixel types.  
	typedef unsigned short             InputPixelType2;
	typedef float          OutputPixelType2;

	//  The input and output image types are instantiated using the pixel types.  
	typedef itk::Image< OutputPixelType2, 2 >   OutputImageType2;
	long int k = 0;
	int max_dist = 0;
	for(int i=0; i<z; i++)
	{
		MyInputImageType2D::Pointer image2D = extract2DImageSlice(im, 2, i);
		typedef itk::SignedMaurerDistanceMapImageFilter<MyInputImageType2D, OutputImageType2>  DTFilter;
		DTFilter::Pointer dt_obj= DTFilter::New() ;
		dt_obj->SetInput(image2D) ;
		dt_obj->SetSquaredDistance( false );      
		dt_obj->SetInsideIsPositive( false );
		try {
			dt_obj->Update() ;
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "Error calculating distance transform: " << err << std::endl ;
			return -1;
		}

		//   Copy the resulting image into the input array  
		typedef itk::ImageRegionIteratorWithIndex< OutputImageType2 > IteratorType;
		IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());
		iterate.GoToBegin();
		unsigned long j=0;	  	 
		while ( !iterate.IsAtEnd() )
		{	  
			double ds = iterate.Get();
			ds = ds*100; //enhance the contrast
			if(ds<=0)
			{
				IMG[k] = 0;
				//iterate.Set(0.0);//try to write back 
			}
			else
				IMG[k] = (unsigned short) ds;
			if(IMG[k]>max_dist)
				max_dist = IMG[k];
			++k; ++j; ++iterate;
		}	  

		//By Yousef: try to write out the output at the central slice
		int cent_slice = (int) z/2;
		if(i==cent_slice)
		{
			typedef    unsigned char     MyInputPixelTypeNew;
			typedef itk::Image< MyInputPixelTypeNew,  2 >   MyInputImageType2DNew;
			typedef itk::CastImageFilter< OutputImageType2, MyInputImageType2DNew> myCasterType;
			myCasterType::Pointer potCaster = myCasterType::New();
			potCaster->SetInput( dt_obj->GetOutput() );
			typedef itk::ImageFileWriter< MyInputImageType2DNew > WriterType;
			WriterType::Pointer writer = WriterType::New();
			writer->SetFileName("dist_test.tif");
			writer->SetInput( potCaster->GetOutput() );
			//writer->Update();		  
		}
	} 

	return max_dist;
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
	filter->SetDirectionCollapseToIdentity();

	filter->SetInput( im );

	MyInputImageType2D::Pointer img = filter->GetOutput();
	try
	{
		img->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}



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
	filter->SetDirectionCollapseToIdentity();

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

	//sort the distances
	unsigned short tmp;
	for(int i=0; i<cntr-1; i++)
	{
		for(int j=i+1; j<cntr; j++)
		{
			if(srtList[j]<srtList[i])
			{
				tmp = srtList[i];
				srtList[i] = srtList[j];
				srtList[j] = tmp;
			}
		}
	}

	//if the number of points is odd then the median is the mid point
	//else, it is in between the two mid points
	int res = cntr % 2;
	int mdn;
	if(res!=0)
	{
		mdn = (cntr+1)/2;
		mdn = (int) srtList[mdn-1];
	}
	else
	{
		int mdn1 = cntr/2;
		mdn = (int) (srtList[mdn1]+srtList[mdn1+1])/2;
	}

	return mdn;
}

int computeWeightedMedian(std::vector< std::vector<float> > scales, int cntr)
{
	if(cntr == 1)
		return scales[0][0];

	float* srtList = new float[cntr];
	float* wgtList = new float[cntr];


	float sumWeights = 0;
	for(int i=0; i<cntr; i++)
	{
		srtList[i] = scales[i][0];
		wgtList[i] = scales[i][1];		
		sumWeights+= wgtList[i];
	}

	//normalize the list of weights
	for(int i=0; i<cntr; i++)
		wgtList[i] /= sumWeights;

	//sort the distances
	for(int i=0; i<cntr-1; i++)
	{
		for(int j=i+1; j<cntr; j++)
		{
			if(srtList[j]<srtList[i])
			{
				float tmp = srtList[i];
				srtList[i] = srtList[j];
				srtList[j] = tmp;
				tmp = wgtList[i];
				wgtList[i] = wgtList[j];
				wgtList[j] = tmp;

			}
		}
	}

	//Find the point at which the cummulative sum exceeds .5
	float cumSum = 0;
	int mdn = -1;
	for(int i=0; i<cntr; i++)
	{
		cumSum+=wgtList[i];
		if(cumSum > .5)
		{
			mdn = srtList[i];
			break;
		}

	}	

	delete [] srtList;
	delete [] wgtList;
	return mdn;
}

void queryOpenCLProperties(float* IM, int r, int c, int z)
{
#ifdef OPENCL
	cout << endl;
	cl_platform_id platforms[10];
	cl_uint num_platforms;
	cl_device_id device[10];
	cl_uint num_devices;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel kernel;
	cl_int errorcode;


	char platform_profile[1024], platform_version[1024], platform_name[1024], platform_vendor[1024], platform_extensions[10240];

	// Platform
	clGetPlatformIDs(10, platforms, &num_platforms);
	cout << "Number of OpenCL platforms:\t" << num_platforms << endl;

	for (unsigned int i = 0; i < num_platforms; i++)
	{
		cout << "Platform " << i+1 << endl;
		clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, sizeof(platform_profile), platform_profile, NULL);
		clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, sizeof(platform_version), platform_version, NULL);
		clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(platform_name), platform_name, NULL);
		clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(platform_vendor), platform_vendor, NULL);
		clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, sizeof(platform_extensions), platform_extensions, NULL);

		cout << "OpenCL Profile:\t\t\t" << platform_profile << endl;
		cout << "OpenCL Platform Version:\t" << platform_version << endl;
		cout << "OpenCL Platform Name:\t\t" << platform_name << endl;
		cout << "OpenCL Platform Vendor:\t\t" << platform_vendor << endl;
		cout << "OpenCL Platform Extensions:\t" << platform_extensions << endl;

		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 10, device, &num_devices);

		cout << "Number of devices: " << num_devices << endl << endl;

		for (unsigned int j = 0; j < num_devices; j++)
		{
			cl_uint device_address_bits, device_global_mem_cacheline_size[1024], device_max_clock_frequency, device_max_compute_units, device_max_constant_args, device_max_read_image_args, device_max_samplers, device_max_work_item_dimensions, \
				device_max_write_image_args, device_mem_base_addr_align, device_min_data_type_align_size, device_preferred_vendor_width_char, device_preferred_vendor_width_short, device_preferred_vendor_width_int, device_preferred_vendor_width_long, \
				device_preferred_vendor_width_float, device_preferred_vendor_width_double, device_vendor_id;

			cl_bool device_available, device_compiler_available, device_endian_little, device_ecc_support, device_image_support;

			cl_device_fp_config device_double_fp_config, device_half_fp_config, device_single_fp_config;

			cl_device_exec_capabilities device_exec_capabilities;

			char device_extensions[1024], device_name[1024], device_profile[1024], device_vendor_2[1024], device_version[1024], driver_version[1024];

			cl_ulong device_global_mem_cache_size, device_global_mem_size, device_local_mem_size, device_max_constant_buffer_size, device_max_mem_alloc_size;

			cl_device_mem_cache_type device_global_mem_cache_type;

			size_t device_image2d_max_height, device_image2d_max_width, device_image3d_max_depth, device_image3d_max_height, device_image3d_max_width, device_max_parameter_size, device_max_work_group_size, device_profiling_timer_resolution;

			size_t device_max_work_item_sizes[1024];

			cl_device_local_mem_type device_local_mem_type;

			cl_platform_id device_platform;

			cl_command_queue_properties device_queue_properties;

			cl_device_type device_type;

			clGetDeviceInfo(device[j], CL_DEVICE_ADDRESS_BITS,					1024, &device_address_bits,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_AVAILABLE,						1024, &device_available,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_COMPILER_AVAILABLE,			1024, &device_compiler_available,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_DOUBLE_FP_CONFIG,				1024, &device_double_fp_config,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_ENDIAN_LITTLE,					1024, &device_endian_little,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_ERROR_CORRECTION_SUPPORT,		1024, &device_ecc_support,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_EXECUTION_CAPABILITIES,		1024, &device_exec_capabilities,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_EXTENSIONS,					1024, &device_extensions,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,			1024, &device_global_mem_cache_size,		NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,			1024, &device_global_mem_cache_type,		NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,		1024, &device_global_mem_cacheline_size,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_GLOBAL_MEM_SIZE,				1024, &device_global_mem_size,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_HALF_FP_CONFIG,				1024, &device_half_fp_config,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE_SUPPORT,					1024, &device_image_support,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE2D_MAX_HEIGHT,			1024, &device_image2d_max_height,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE2D_MAX_WIDTH,				1024, &device_image2d_max_width,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE3D_MAX_DEPTH,				1024, &device_image3d_max_depth,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE3D_MAX_HEIGHT,			1024, &device_image3d_max_height,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_IMAGE3D_MAX_WIDTH,				1024, &device_image3d_max_width,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_LOCAL_MEM_SIZE,				1024, &device_local_mem_size,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_CLOCK_FREQUENCY,			1024, &device_max_clock_frequency,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_COMPUTE_UNITS,				1024, &device_max_compute_units,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_CONSTANT_ARGS,				1024, &device_max_constant_args,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,		1024, &device_max_constant_buffer_size,		NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_MEM_ALLOC_SIZE,			1024, &device_max_mem_alloc_size,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_PARAMETER_SIZE,			1024, &device_max_parameter_size,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_READ_IMAGE_ARGS,			1024, &device_max_read_image_args,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_SAMPLERS,					1024, &device_max_samplers,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_GROUP_SIZE,			1024, &device_max_work_group_size,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,		1024, &device_max_work_item_dimensions,		NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_ITEM_SIZES,			1024, &device_max_work_item_sizes,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MAX_WRITE_IMAGE_ARGS,			1024, &device_max_write_image_args,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MEM_BASE_ADDR_ALIGN,			1024, &device_mem_base_addr_align,			NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE,		1024, &device_min_data_type_align_size,		NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_NAME,							1024, &device_name,							NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PLATFORM,						1024, &device_platform,						NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,	1024, &device_preferred_vendor_width_char,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,	1024, &device_preferred_vendor_width_short,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,	1024, &device_preferred_vendor_width_int,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,	1024, &device_preferred_vendor_width_long,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,	1024, &device_preferred_vendor_width_float,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, 1024, &device_preferred_vendor_width_double,NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PROFILE,						1024, &device_profile,						NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_PROFILING_TIMER_RESOLUTION,	1024, &device_profiling_timer_resolution,	NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_QUEUE_PROPERTIES,				1024, &device_queue_properties,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_SINGLE_FP_CONFIG,				1024, &device_single_fp_config,				NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_TYPE,							1024, &device_type,							NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_VENDOR,						1024, &device_vendor_2,						NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_VENDOR_ID,						1024, &device_vendor_id,					NULL);
			clGetDeviceInfo(device[j], CL_DEVICE_VERSION,						1024, &device_version,						NULL);
			clGetDeviceInfo(device[j], CL_DRIVER_VERSION,						1024, &driver_version,						NULL);

			cout << "Device Name:\t\t\t" << device_name << endl;
			cout << "Device Version:\t\t\t" << device_version << endl;
			cout << "Driver Version:\t\t\t" << driver_version << endl;
			cout << "Device Address Bits:\t\t" << device_address_bits << endl;
			cout << "Device ECC Support:\t\t" << device_ecc_support << endl;
			cout << "Device Max Clock Frequency:\t" << device_max_clock_frequency << endl;
			cout << "Device Max Compute Units:\t" << device_max_compute_units << endl;
			cout << "Device Global Memory:\t\t" << device_global_mem_size / (1024*1024) << " MB" << endl;
			cout << "Device Local Memory:\t\t" << device_local_mem_size / 1024 << " KB" << endl;
			cout << "Device Max Mem Alloc:\t\t" << device_max_mem_alloc_size  / (1024 * 1024) << " MB" << endl;
			cout << "Device Max Constant Mem:\t" << device_max_constant_buffer_size / 1024 << " KB" << endl;
			cout << "Device Extensions:\t\t" << device_extensions << endl;
			cout << "Timer Resolution:\t\t" << device_profiling_timer_resolution << " ns" << endl;

			cout << endl;
		}

		cout << stringifiedKernel << endl << endl;

		context = clCreateContext(0, 1, device, &seed_pfn_notify, NULL, NULL);
		queue = clCreateCommandQueue(context, device[0], 0, NULL);
		program = clCreateProgramWithSource(context, 1, (const char **) &stringifiedKernel, 0, &errorcode);

		cout << "clCreateProgramWithSource Error code: " << errorcode << endl;

		errorcode = clBuildProgram(program, 0, 0, 0, 0, 0);

		char build_log[1000];

		clGetProgramBuildInfo(program, device[0], CL_PROGRAM_BUILD_LOG, sizeof(build_log), build_log, NULL);

		cout << "Build Log:" << endl << build_log << endl;

		kernel = clCreateKernel(program, "fibonacci", &errorcode);

		cout << "clCreateKernel Error code: " << errorcode << endl;

		size_t cnDimension = 128;

		unsigned long long* testArray = new unsigned long long[cnDimension];

		for (int i = 0; i < cnDimension; i++)
			testArray[i] = 1;

		cout << "Size of array element = " << sizeof(*testArray) << endl;
		cout << "Size of cl_double = " << sizeof(cl_double) << endl;
		cout << "Allocating " << (sizeof(*testArray) * cnDimension)/(double)(1024) << " KB of memory for matrix" << endl;

		cl_mem deviceMem = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ulong) * cnDimension, NULL, NULL);

		if (deviceMem == NULL)
			cout << "Failed to allocate buffer memory" << endl; 
		else
			cout << "Memory allocation sucessful" << endl;		

		clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &deviceMem);

		clEnqueueNDRangeKernel(queue, kernel, 1, 0, (const size_t *) &cnDimension, 0, 0, 0, 0);

		clEnqueueReadBuffer(queue, deviceMem, CL_TRUE, 0, sizeof(cl_ulong) * cnDimension, testArray, NULL, NULL, NULL);

		clFinish(queue);

		cout << "Results" << endl;
		for (int i = 0; i < cnDimension; i++)
			cout << testArray[i] << endl;

		cout << endl;

		delete[] testArray;
		clReleaseKernel(kernel);
		clReleaseProgram(program);
		clReleaseMemObject(deviceMem);
		clReleaseCommandQueue(queue);
		clReleaseContext(context);

	}

	cout << endl;
#endif
}

void seed_pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
	std::cerr << errinfo << std::endl;
}


