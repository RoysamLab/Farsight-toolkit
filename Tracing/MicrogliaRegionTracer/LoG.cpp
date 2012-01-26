#include "LoG.h"

LoG::LoG()
{
}

LoG::LoGImageType::Pointer LoG::RunLoG(ImageType::Pointer image, float scale)
{
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType , LoGImageType> LoGFilterType;
	LoGFilterType::Pointer LoGFilter = LoGFilterType::New();
	LoGFilter->SetInput( image );
	LoGFilter->SetNormalizeAcrossScale(true);
	LoGFilter->SetSigma( scale );

	try
	{
		LoGFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "LoGFilter Exception: " << err << std::endl;
	}


	//Inversion filter, since the peak response in LoG is negative, we want to convert this to positive
	typedef itk::ShiftScaleImageFilter<LoGImageType, LoGImageType> InvertFilterType;
	InvertFilterType::Pointer invertFilter = InvertFilterType::New();
	invertFilter->SetScale(-1.0f);
	invertFilter->SetInput(LoGFilter->GetOutput());

	try
	{
		invertFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "invertFilter Exception: " << err << std::endl;
	}

	return invertFilter->GetOutput();
}

void LoG::WriteLoGImage(std::string filename, LoGImageType::Pointer image)
{
	typedef itk::ImageFileWriter< LoGImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(filename);
	try
	{
		std::cout << "Writing LoG image" << std::endl;
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

std::vector<LoG::LoGImageType::Pointer> LoG::RunMultiScaleLoG(Seed* seed, ImageType::Pointer seed_image)
{
	std::vector<LoGImageType::Pointer> multiscale_LoG_vector;

	for (float scale = 0.5; scale <= 2.0; scale+=0.5)
	{
		LoG *log_obj = new LoG();
		LoGImageType::Pointer LoGimage = log_obj->RunLoG(seed_image, scale);

		multiscale_LoG_vector.push_back(LoGimage);
	}

	return multiscale_LoG_vector;
}