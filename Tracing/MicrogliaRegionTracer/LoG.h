#ifndef LoG_H
#define LoG_H

#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageFileWriter.h"

#include "Seed.h"

class LoG
{
public:
	typedef itk::Image<unsigned char, 3> ImageType;
	typedef itk::Image<float, 3> LoGImageType;
public:
	LoG();
	~LoG();

	LoGImageType::Pointer RunLoG(ImageType::Pointer image, float scale);
	void WriteLoGImage(std::string filename, LoGImageType::Pointer image);

	std::vector<LoGImageType::Pointer> RunMultiScaleLoG(Seed* seed, ImageType::Pointer seed_image);
};

#endif
