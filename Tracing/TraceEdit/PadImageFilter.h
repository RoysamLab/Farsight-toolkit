#ifndef PADIMAGEFILTER_H
#define PADIMAGEFILTER_H

#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageRegionIterator.h"

typedef itk::Image< unsigned char, 3 >  ImageType;
typedef itk::ConstantPadImageFilter< ImageType, ImageType > ConstantPadImageFilterType;

class PadImageFilter
{
public:
	PadImageFilter(void);
	~PadImageFilter(void);

	void padImage(ImageType::Pointer input, ImageType::Pointer &output, int extendSize[]);
};
#endif