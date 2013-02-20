#include "Hessian.h"

Hessian::Hessian()
{
}

Hessian::HessianImageType::Pointer Hessian::RunHessian(ImageType::Pointer image, float scale)
{
	typedef itk::HessianRecursiveGaussianImageFilter< ImageType , HessianImageType > HessianFilterType;
	HessianFilterType::Pointer HessianFilter = HessianFilterType::New();
	HessianFilter->SetInput( image );
	HessianFilter->SetNormalizeAcrossScale(true);
	HessianFilter->SetSigma( scale );

	try
	{
		HessianFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "HessianFilter Exception: " << err << std::endl;
	}

	return HessianFilter->GetOutput();
}