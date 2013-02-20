#ifndef Hessian_H
#define Hessian_H

#include "itkHessianRecursiveGaussianImageFilter.h"

#include "Cell.h"

class Hessian
{
private:
	typedef Cell::ImageType ImageType;
	typedef Cell::HessianImageType HessianImageType;

public:
	Hessian();
	~Hessian();

	HessianImageType::Pointer RunHessian(ImageType::Pointer image, float scale);
};

#endif
