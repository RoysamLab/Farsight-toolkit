

#include "itkRGBToHSVColorSpacePixelAccessor.h" //This is pulled off the IJ -- NOT AN ACTUAL ITK CLASS
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include <itkIndex.h>
#include "itkImageFileReader.h"
#include <itkScalarImageKmeansImageFilter.h>
#include "itkIsolatedConnectedImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkRelabelComponentImageFilter.h>

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef unsigned char CharPixelType;
typedef unsigned short ShortPixelType;
typedef itk::Image<CharPixelType, 3> CharImageType;
typedef itk::Image<ShortPixelType,3> ShortImageType;
typedef CharImageType ImageType;
//typedef itk::Image<float , 2> floatImageType;
typedef float FloatPixelType;
typedef signed short SShortPixelType;
typedef itk::Vector<float, 3> VectorPixelType;
typedef itk::Image<RGBPixelType, 3> RGBImageType;
typedef itk::Image<FloatPixelType, 3> FloatImageType;
typedef itk::ImageFileReader<RGBImageType> ImageFileReaderType;
typedef itk::ImageFileWriter<ShortImageType> ShortImageWriterType;
typedef itk::ImageFileWriter<CharImageType> CharImageWriterType;


class Daniela
{
	public:
		Daniela();
		Daniela(std::string imageFileName,int numberOfInitialClasses);
		~Daniela();
		
		
		int Get_Daniela();	//call to calculate the bone_label_image
		ImageType::Pointer bone_label_image;		//labeled segmented bone image
		RGBImageType::Pointer rgbImage ;
		void Read_RGBInputImage();

	
	private:
		int M, N;
		float Set_S_Image(ImageType::Pointer inputimage);
		
		ImageType::Pointer GzetImage(ImageType::Pointer inputimage );

		std::string imageFileName;
		int numberOfInitialClasses;
        

};