#ifndef LoG_H
#define LoG_H

#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageFileWriter.h"

#include "fregl/fregl_roi.h"
#include "Cell.h"

class LoG
{
private:
	typedef fregl_roi::ImageType ImageType;
	typedef Cell::LoGImageType LoGImageType;

public:
	LoG();
	~LoG();

	LoGImageType::Pointer RunLoG(ImageType::Pointer image, float scale);
	void WriteLoGImage(std::string filename, LoGImageType::Pointer image);

	std::vector<LoGImageType::Pointer> RunMultiScaleLoG(Cell* cell);
};

#endif
