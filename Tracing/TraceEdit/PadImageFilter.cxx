#include "PadImageFilter.h"

PadImageFilter::PadImageFilter(void)
{
}

PadImageFilter::~PadImageFilter(void)
{
}

void PadImageFilter::padImage(ImageType::Pointer input, ImageType::Pointer &output, int extendSize[])
{
	/*!
	 * a circle kernel is used to evaluate the voxel intensities along the centerline
	 * @author Audrey Cheong
	 * @param input itk image
	 * @param output return padded itk image with black (zero-valued) voxels
	 * @param extendSize set 6 parameters [x_lower_bound x_upper_bound y_lower_bound y_upper_bound z_lower_bound z_upper_bound]
	 */

	ImageType::SizeType lowerExtendRegion;
	lowerExtendRegion[0] = extendSize[0];
	lowerExtendRegion[1] = extendSize[2];
	lowerExtendRegion[2] = extendSize[4];

	ImageType::SizeType upperExtendRegion;
	upperExtendRegion[0] = extendSize[1];
	upperExtendRegion[1] = extendSize[3];
	upperExtendRegion[2] = extendSize[5];

	ImageType::PixelType constantPixel = 0;

	ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
	padFilter->SetInput(input);
	padFilter->SetPadLowerBound(lowerExtendRegion);
	padFilter->SetPadUpperBound(upperExtendRegion);
	padFilter->SetConstant(constantPixel);
	padFilter->Update();

	output = padFilter->GetOutput();
}