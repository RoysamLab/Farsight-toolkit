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


	//Scale to (-1, 1) range and invert
	typedef itk::ShiftScaleImageFilter<LoGImageType, LoGImageType> InvertFilterType;
	InvertFilterType::Pointer invertFilter = InvertFilterType::New();
	invertFilter->SetScale(-1.0 / std::numeric_limits<ImageType::PixelType>::max());
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

std::vector<LoG::LoGImageType::Pointer> LoG::RunMultiScaleLoG(Cell* cell, ImageType::Pointer cell_image)
{
	std::vector<LoGImageType::Pointer> multiscale_LoG_vector;

	for (float scale = 1; scale <= 3; scale+=0.25)
	{
		LoG *log_obj = new LoG();
		LoGImageType::Pointer LoGimage = log_obj->RunLoG(cell_image, scale);

		multiscale_LoG_vector.push_back(LoGimage);
	}

	return multiscale_LoG_vector;
}