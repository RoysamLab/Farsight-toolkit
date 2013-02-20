#ifndef _FOCUS_H_
#define _FOCUS_H_

//ITK Includes
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVarianceImageFunction.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkSubtractImageFilter.h"


//STL and Other Includes
#include <vector>

typedef float FloatPixelType;
typedef unsigned char UCharPixelType;
typedef unsigned short UShortPixelType;
typedef itk::RGBPixel< UCharPixelType > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 > RGBImageType;
typedef itk::Image< RGBPixelType, 2 > RGBImageType2D;
typedef itk::Image< UCharPixelType, 3> UCharImageType;
typedef itk::Image< UCharPixelType, 2> UCharImageType2D;
typedef itk::Image< FloatPixelType, 3> FloatImageType;
typedef itk::Image< UShortPixelType, 2> UShortImageType2D;

class Focus
{
public:
	//Constructor
	Focus(UCharImageType::Pointer input);
	Focus(RGBImageType::Pointer input);
	~Focus(){};//Destructor

	void SetRadius(int r){ radius = r; };

	void MakeVarianceImage();
	RGBImageType2D::Pointer MakeProjectionColor();
	UCharImageType2D::Pointer MakeProjection();
	std::pair<UCharImageType2D::Pointer,UShortImageType2D::Pointer> MakeProjection2();

	std::vector<float> FindVariance(int x, int y);
	std::vector<float> FindVariance(float x, float y, float ppi); 

protected:
	UCharImageType::Pointer imgIn;
	RGBImageType::Pointer imgInColor;
	FloatImageType::Pointer varImg;

	int radius;
};

#endif
