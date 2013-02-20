#ifndef __ROLLING_BALL_FILTER_H
#define __ROLLING_BALL_FILTER_H

#include "itkSubtractImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkStreamingImageFilter.h"

class RollingBallFilter
{
public:
	//Some typedefs to make things easier to read later
	const static int ImageDimensions = 3;
	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, ImageDimensions >	ImageType;

private:	
	ImageType::Pointer input_img;
	ImageType::Pointer output_img;
	float radius;

public:	
	RollingBallFilter(ImageType::Pointer, float radius);	//Constructor
	RollingBallFilter(char* fileName, float radius);		//Constructor

	void				RunFilter();																																					//Run the filter
	ImageType::Pointer	GetOutput();	//Get the output, the original image is not affected
	void WriteOutput(char* fileName);	//Write the output to a file

private:
	//Purposely put in private as no one should be calling default constructor/destructor
	RollingBallFilter();	
	~RollingBallFilter();
};

#endif