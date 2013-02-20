#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include <itkIndex.h>
#include "itkImageFileReader.h"


//class Haralick returns  the Haralick texture feature vector


class Haralick
{
	public:
		
		#define Dimension 2
		#define VIndexDimension 2

		typedef unsigned char PixelType;
		typedef itk::Image<PixelType , Dimension > ImageType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
		typedef itk::Index< VIndexDimension > IndexType;
		typedef ImageType::SizeType InputSizeType;
		std::string imageFileName;
		ImageType::Pointer inputImage;
		Haralick(std::string imageFileName);
		Haralick(ImageType::Pointer inputImage);
		Haralick();
		~Haralick();
		double * GetHaralik();
		

	private:		

		
    	unsigned char ** f;
	    int M , N ,Ng;	
		int ReadInputImage(std::string imageFileName);
		void SetInputImage();
		// this function calculates the maximum intensity value of image
		int  max_value();
		// this method calculates minimum value of inensity value in image
		int  min_value();
		// the method sets the size of the cooccurence matrix 
		void set_Ng();
		double f1_ASM,f2_Contrast,f3_Correlation,f4_Variance,f5_IDM,f6_Sum_Avg,f7_Sum_Var,f8_Sum_Entropy,f9_Entropy,f10_diff_var,f11_difference_entropy,f12,f13;
		double f14_maxcorr (int** P,int  N); 
		double *pgm_vector (int nl, int nh);
		double **pgm_matrix (int nrl,int  nrh, int ncl, int nch);
		void results (double *Tp,char *c, double *a);
		void simplesrt ( int n,  double arr[]);
		void mkbalanced (  double **a,  int n);
		void reduction ( double **a,  int n);
		int hessenberg ( double **a,double wr[],double wi[],  int n);
		void  SWAP(double a,double b);
		
		
		int ** calculateco_occurrence( int dir[]);
		//this method calculates the cooccurance matrix for a given angle .As the cooccurrence matrix is symmetric 
		//(used by haralick feature calculation) so the matrices coorsponding to angles theta and theta+180 are combined
		void calculatecooccurrence( int** P ,int angle);
		void calculatefeatures(int** P );
		
	

};