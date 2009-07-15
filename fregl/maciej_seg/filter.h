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

/** @file filter.h
*   @brief class for filtering
*   This is the class that filter the image prior to merging
*
*   @author Maciej Wotjon
*/

#ifndef __FILTER_H_
#define __FILTER_H_

// FarSight include
#include "pix_t.h"
#include "cell.h"

// STL
#include <vector>
#include <cmath>
#include <string>
#include <numeric>
#include <functional>
#include <fstream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkMedianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageAdaptor.h"

#include <stdio.h>
#include <stdlib.h>	//for reading in gangs data

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_det.h>
#include <vnl/algo/vnl_matrix_inverse.h>

/**	@brief helper class for green channel accesor
*/
class GreenChannelPixelAccessor  
{
public:
	typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef               unsigned char    ExternalType;

	static ExternalType Get( const InternalType & input ) 
	{
		return static_cast<ExternalType>( input.GetGreen() );
	}
};

/**	@brief helper class for red channel accesor
*/
class RedChannelPixelAccessor  
{
public:
	typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef               unsigned char    ExternalType;

	static ExternalType Get( const InternalType & input ) 
	{
		return static_cast<ExternalType>( input.GetRed() );
	}
};

/**	@brief helper class for blue channel accesor
*/
class BlueChannelPixelAccessor  
{
public:
	typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef               unsigned char    ExternalType;

	static ExternalType Get( const InternalType & input ) 
	{
		return static_cast<ExternalType>( input.GetBlue() );
	}
};

#if defined(_MSC_VER) && (_MSC_VER < 1300)
static const unsigned int msDimension = 3;
#endif

class filter
{
public:
	//msDimension
#if !(defined(_MSC_VER) && (_MSC_VER < 1300))
	static const unsigned int msDimension = 3;
#endif
	//pixel types for filters
	typedef float  mInputPixelType;
	typedef float  mOutputPixelType;

	//pixel type for adaptor input
	typedef GreenChannelPixelAccessor::InternalType  mInputGreenPixelFileType;
	typedef RedChannelPixelAccessor::InternalType  mInputRedPixelFileType;
	typedef BlueChannelPixelAccessor::InternalType  mInputBluePixelFileType;

	//pixel type for file output
	typedef unsigned char mOutputPixelFileType;

	//image type for filters
	typedef itk::Image< mInputPixelType,  msDimension >   mInputImageType;
	typedef itk::Image< mOutputPixelType, msDimension >   mOutputImageType;

	//image type for file output
	typedef itk::Image< mOutputPixelFileType, msDimension >   mOutputImageFileType;

	//image type for adaptor input
	typedef itk::Image< mInputGreenPixelFileType, msDimension >   mInputImageFileType;

	//reader type
	typedef itk::ImageSeriesReader< mInputImageFileType  >  mReaderType;
        typedef itk::ImageFileReader< mInputImageFileType  >  mReaderFileType;

	//adaptor type for green channel
	typedef itk::ImageAdaptor<  mInputImageFileType, 
		GreenChannelPixelAccessor > mImageGreenAdaptorType;
	typedef itk::ImageAdaptor<  mInputImageFileType, 
		RedChannelPixelAccessor > mImageRedAdaptorType;
	typedef itk::ImageAdaptor<  mInputImageFileType, 
		BlueChannelPixelAccessor > mImageBlueAdaptorType;

	//create rescaler for adaptor
	typedef itk::RescaleIntensityImageFilter< mImageGreenAdaptorType, 
		mOutputImageType 
	>   mRescalerGreenAdaptorType;

	typedef itk::RescaleIntensityImageFilter< mImageBlueAdaptorType, 
		mOutputImageType 
	>   mRescalerBlueAdaptorType;

	typedef itk::RescaleIntensityImageFilter< mImageRedAdaptorType, 
		mOutputImageType 
	>   mRescalerRedAdaptorType;

	//declare iterator types
	typedef itk::ImageRegionConstIteratorWithIndex< mOutputImageType > mConstIndexIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< mOutputImageType > mIndexIteratorType;
	typedef itk::ImageRegionConstIterator< mOutputImageType > mConstIteratorType;
	typedef itk::ImageRegionIterator< mOutputImageType > mIteratorType;

	//medium filter
	typedef itk::MedianImageFilter<
		mOutputImageType, mOutputImageType >  mMedianFilterType;

	//structuring element for open filter
	typedef itk::BinaryBallStructuringElement< 
		mOutputPixelType,
		msDimension  >             mStructuringElementType;

	//morphological opening filter
	typedef itk::GrayscaleMorphologicalOpeningImageFilter<
		mOutputImageType, 
		mOutputImageType, 
		mStructuringElementType >  mOpenFilterType;

	//filter types for erode, dialate, open and subtract filters 
	typedef itk::GrayscaleErodeImageFilter<
		mOutputImageType, 
		mOutputImageType,
		mStructuringElementType >  mErodeFilterType;

	typedef itk::GrayscaleDilateImageFilter<
		mOutputImageType, 
		mOutputImageType, 
		mStructuringElementType >  mDilateFilterType;

	typedef itk::SubtractImageFilter<
		mOutputImageType,
		mOutputImageType,
		mOutputImageType  >       mSubFilterType;

	//region of interest fitler
	typedef itk::RegionOfInterestImageFilter< mOutputImageType, 
		mOutputImageType > mRegionFilterType;

	//binary threshold filter
	typedef itk::BinaryThresholdImageFilter<
		mOutputImageType, mOutputImageType >  mThreshFilterType;

	//std::vector to hold miIndex of background pixels
	std::vector<mOutputImageType::IndexType> mvIndex;

	//create distance map filter
	typedef itk::DanielssonDistanceMapImageFilter<
		mOutputImageType, mOutputImageType >  mDanielFilterType;

	//create gaussian filter
	typedef itk::DiscreteGaussianImageFilter<mOutputImageType
		,mOutputImageType> mGaussFilterType;

	typedef itk::ImageFileWriter< mOutputImageFileType >  mWriterType;
	typedef itk::ImageFileWriter< mOutputImageType >  mWriterImageType;

	//rescale final image to be saved
	typedef itk::RescaleIntensityImageFilter< 
		mOutputImageType, mOutputImageFileType > mRescalerFilterType;

protected:
	//points for segmentation
	std::vector<pix_t*> mvPix;
	//region for ITK images
	mOutputImageType::RegionType mRegion;
	//size of image
	int m_nImage_size[3];
	//temp images
	mInputImageFileType::Pointer mImgpRGBimage;
	mOutputImageType::Pointer mImgpIntenisty;
	mOutputImageType::Pointer mImgpGradient;
	mOutputImageType::Pointer mImgpFilt;

	int m_gridX;	//size of grid for thresholding
	int m_gridY;
	int m_radMed;	//radius of med filter
	int m_radGrad;	//rad of grad filter
	double m_paramT;	//param for threshilding filter
	double m_sigma;	//param for gauss filter

	int mChannel;

public:

	filter(int gX=1, int gY=1, int rM=1, int rG=1, double pT=.4, double s=.5,int channel=1);

	~filter(void);

	mOutputImageType::Pointer load_green_channel(const std::string &crsName,const std::string &crsPath,const std::string &crsType,
		const int &crStart,const int &crEnd,const int &crWidth);
  mOutputImageType::Pointer load_green_channel(const std::string &filename);

	std::vector<std::string> gen_names(const std::string &path, const std::string &file_first,const std::string &type,int start,int end,int width);

	void save_intensity(mOutputImageType::Pointer img);

	mOutputImageType::Pointer median_open_image(mOutputImageType::Pointer img,int radius);

	void get_min_max_mean_std(mOutputImageType::Pointer image, double &max, double &min, double &mean, double &std,bool std_flag);

	mOutputImageType::Pointer morph_gradient(mOutputImageType::Pointer img,int radius);

	void save_gradient(mOutputImageType::Pointer img);

	mOutputImageType::Pointer threshold_distance_map(mOutputImageType::Pointer img,int hor,int ver,int dep,double param);

	mOutputImageType::Pointer combine_invert_distance(mOutputImageType::Pointer imgDaniel,mOutputImageType::Pointer imgMorphGrad,double sigma);

	mInputImageFileType::Pointer save_image_rgb(std::string file);

	void save_image(mOutputImageType::Pointer img,std::string file);

	void run_filter(const std::string &crsName,const std::string &crsPath,const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth);
        void run_filter(const std::string & crsName);

	void label_image(std::vector<cell*> cells);

	void label_bound(std::vector<cell*> cells);

	mOutputImageType::Pointer getintensity(void);

	mOutputImageType::Pointer getfilt(void);

	mOutputImageType::Pointer getgradient(void);

	int* getsize(void);

	void load_images(std::string sFilt, std::string sGrad, std::string sInt, std::string sPath, std::string sRGB,
					std::string sType,const int &start,const int &end,const int &width,int &loadimages);

	void clear(void);

	mOutputImageType::Pointer load_image(std::string sName);

	mOutputImageType::Pointer load_channel(const std::string &crsName,const std::string &crsPath,const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth);
  mOutputImageType::Pointer load_channel(std::string const & filename);

	mOutputImageType::Pointer load_blue_channel(const std::string &crsName,const std::string &crsPath,const std::string &crsType,
		const int &crStart,const int &crEnd,const int &crWidth);
	mOutputImageType::Pointer load_blue_channel(const std::string &filename);

	mOutputImageType::Pointer load_red_channel(const std::string &crsName,const std::string &crsPath,const std::string &crsType,
		const int &crStart,const int &crEnd,const int &crWidth);
	mOutputImageType::Pointer load_red_channel(const std::string &filename);
  

        mInputImageFileType::Pointer project_image(std::string &file,std::vector<cell*> cells);

};

///**	@brief class to get gangs data saved into image
//*/
//class data_rec: public filter
//{
//public:
//	data_rec()
//	{
//
//	}
//	filter::mOutputImageType::Pointer
//	get_image(void)
//	{
//		mOutputImageType::Pointer point = load_green_channel("mcvh1_vh4-24-03_arc_2_3cf_ca3_","../Images/",
//		"tif",7,22,2);	//load image to save
//
//		mOutputImageType::IndexType miIndex;
//
//		itk::ImageRegionIterator< mOutputImageType > it(point,point->GetLargestPossibleRegion());
//		FILE *outfile;
//		char *outname = "mcvh1_vh4-24-03_arc_2_3cf_ca3_seg_final.dat";
//		unsigned short num;
//		outfile = fopen(outname, "rb");	//set read and bin
//
//		for(it.GoToBegin();!it.IsAtEnd();++it)	//reset iamage
//			it.Set(0);
//
//		int count = 0;
//
//
//
//		//for(int y=0;y<512;y++)
//		//{
//		for(int y=511;y>-1;y--)
//		{
//			for(int x=0;x<512;x++)
//			{
//				for(int z=0;z<16;z++)
//				{ 
//					fread(&num, 2, 1, outfile);	//read 2 bytes
//					//in.read(&num,2);
//					miIndex[0] = x;
//					miIndex[1] = y;
//					miIndex[2] = z;
//					it.SetIndex(miIndex);	//save into image
//					if( num == 0 )
//						it.Set(-1);
//					else
//						it.Set(num);
//					count++;
//
//				}
//			}
//		}
//		std::cout << count << std::endl;
//		return point;
//	}
//};
//
#endif
