#ifndef LoG_H
#define LoG_H

#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageFileWriter.h"

#include "Cell.h"

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

	std::vector<LoGImageType::Pointer> RunMultiScaleLoG(Cell* cell, ImageType::Pointer cell_image);
};

#endif
