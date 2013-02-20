#ifndef __ftkSpectralUnmixing_h
#define __ftkSpectralUnmixing_h
//includes

#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkMedianImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkVector.h>
#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_hungarian_algorithm.h>
#include <vnl/algo/vnl_qr.h>

#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkUtils.h>
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

#define NUM_ELEMENTS_FOR_UNMIXING 9000000
#define MIN_NORM 30
#define MAX_CHANNS 10


namespace ftk
{

static double square_function (double a);

class SpectralUnmixing
{
	public:
		SpectralUnmixing();
		~SpectralUnmixing();
		
		/*Image and Iterators Typedefs*/
		typedef unsigned char InputPixelType;
		typedef unsigned char OutputPixelType;
		typedef itk::Image<InputPixelType,3> InputImageType;
		typedef itk::Image<OutputPixelType,3> OutputImageType;
		typedef itk::Image<InputPixelType,2> Input2DImageType;
		typedef itk::Image<InputPixelType,2> Output2DImageType;
		typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
		typedef itk::ImageRegionIterator<InputImageType> IteratorType;
		typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
		typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;
		typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
		typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;
		typedef itk::MedianImageFilter<InputImageType,InputImageType> MedianFilterType;
		typedef enum { PIXEL_CLASSIFICATION, LINEAR_UNMIX } UnmixMode;

		void SetInputImage(ftk::Image::Pointer image);
		void SetNumberOfChannels(int m);
		void SetUnmixMode(UnmixMode mode = PIXEL_CLASSIFICATION);			// Default is pixel classification
		void Update(void);
		ftk::Image::Pointer GetOutput(void){return this->UnmixedImage;};

	private:
		// Data:
		typedef enum { UNDER_DETERMINED, OVER_DETERMINED } SystemMode;
		UnmixMode unmixMode;
		SystemMode sysMode;
		ftk::Image::Pointer Image;
		ftk::Image::Pointer UnmixedImage;
		int MChannels;							// number of output channels
		int NChannels;							// number of input channels
		std::vector<std::vector<InputImageType::Pointer> >Unmixed_Images;

		// Functions:
		//vnl_matrix<double> GetFingerPrintMatrix(InputImageType::Pointer im[]);
		vnl_matrix<double> GetFingerPrintMatrix(void);
		void EstimateFingerPrintMatrix(vnl_matrix<double> mixed, vnl_matrix<double> &start, vnl_vector<unsigned char> &indices);
		void UnmixPureChannels(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start);	// this unmixing method uses voxel classification based on maximum projection onto cluster centers
		void UnmixUsingIterations(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start);// this method is for over-determined/full ranked systems MChannels<=NChannels
		void UnmixUsingPseudoInverse(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> &start);// this method is for over-determined/full ranked systems MChannels<=NChannels
	//	void UnmixUsingIterations(void);																			// this method is for under-determined systems where MChannels>NChannels
		void UnmixClustering(InputImageType::Pointer im[],InputImageType::Pointer om[],vnl_matrix<double> start);
		void ConvertOutputToftk(void);
		void getColor(int numChann,std::vector<unsigned char> *channelColors);
};// end of class

}// end of namespace
#endif







