#ifndef __Zernike_h
#define __Zernike_h

#include <iostream>
#include <math.h>

#include "itkImage.h"
#include "itkIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkImageRegionConstIterator.h"
//#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"


///////////////////////
/*THE CLASS ZERNIKE RETURNS THE ZERNIKE MOMENTS FOR GRAY SCALE IMAGES
Now,zernike moments are orthogonal on the unit circle(r==1) */
////////////////////////////


namespace ftk{
class zernike
{
	public:

		#define Dimen 2
		#define VIndexDimension 2

		typedef unsigned char PixelType;
		typedef itk::Image<PixelType , Dimen > ImageType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
		typedef itk::Index< VIndexDimension > IndexType;
		typedef ImageType::SizeType InputSizeType;
		
		int orderofmoments;
		std::string imageFileName;
		ImageType::Pointer inputImage;
		// if the input is an image file name from the working directory
		zernike(std::string imageFileName,int orderofmoments);
		//if the input is  an ITK image pointer then call this constructor
		zernike(ImageType::Pointer inputImage,int orderofmoments);
		zernike();
		~zernike();		
		std::vector< std::vector<double> > GetZernike();
		
    private: 	    	
	   
		int ReadInputImage(std::string imageFileName);
		// this method is used to set the input image parameters using  Image Region Iterators as we need image dimensions (M and N)  for the calculations of Zernike moments and need to access the image pixel data
		void SetInputImage();
		int M , N ;
		unsigned char ** f;
	
		// this methog gives the values of X coordinates in the  image
		double* getX();
		// this methog gives the values of Y coordinates in the mapped image
		double* getY();
		//this method returns the (0,0) order geometrical moment which is equal to sum(sum(f(x,y));
		double G00();
		//this method returns the (1,0) order geometrical moment which is needed to calculate centre of floroscence for the image ( i.e. X component of  centre of mass in image is is G10()/G00()); 

		double COF_X();
		//this method returns the (0,1) order geometrical moment which is needed to calculate the  centre of fluoroscence for the image ( i.e. Y component of  centre of mass in image is is G01()/G00()); 
		double COF_Y();
		double beta();
		//Rho -> go to the normalized POLAR  coordinates (RHO AND THETA)
		double**  RHO();
		double** THETA();
		//to get factorial of integers
		double  factorial (int a);
		//to get the gernike moments of order p and repetition q
		double  B(int p, int q , int k );
		double* CalculateZernike(int p , int q);

};
} //end namespace ftk


#endif
