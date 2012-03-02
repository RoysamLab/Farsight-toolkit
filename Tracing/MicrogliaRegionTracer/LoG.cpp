#include "LoG.h"

LoG::LoG()
{
}

LoG::LoGImageType::Pointer LoG::RunLoG(ImageType::Pointer image, float scale)
{
	//typedef itk::MedianImageFilter<ImageType, ImageType> MedianImageFilterType;
	//MedianImageFilterType::InputSizeType radius;
	//radius.Fill(3);

	//MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
	//medianFilter->SetInput(image);
	//medianFilter->SetRadius(radius);
	
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
		throw; //Rethrowing exception, since only RunMultiScaleLoG knows the location of the cell
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

LoG::LoGImageType::Pointer LoG::RunMultiScaleLoG(Cell* cell)
{
	std::vector<LoGImageType::Pointer> LoG_vector;

	//Calculate all the LoG scales
	for (float scale = 1.5; scale <= 2.5; scale+=0.1)
	{
		LoGImageType::Pointer LoGimage;

		try
		{
			LoGimage = RunLoG(cell->image, scale);
		}
		catch (itk::ExceptionObject &err)
		{
			ImageType::PointType origin = cell->image->GetOrigin();
			ImageType::SizeType size = cell->image->GetLargestPossibleRegion().GetSize();
			
			std::cerr << "RunMultiScaleLoG exception: " << std::endl;
			std::cerr << "For cell: " << cell->getX() << ", " << cell->getY() << ", " << cell->getZ() << " at scale: " << scale << " Origin: " << origin << " Size: " << size << std::endl;
		}

		LoG_vector.push_back(LoGimage);
	}

	LoGImageType::SizeType size = cell->image->GetLargestPossibleRegion().GetSize();

	//Make a new image to store the critical points	
	LoGImageType::Pointer multiscale_LoG_image = LoGImageType::New();
	LoGImageType::IndexType start;
	start.Fill(0);

	LoGImageType::RegionType region(start, size);
	multiscale_LoG_image->SetRegions(region);
	multiscale_LoG_image->Allocate();
	multiscale_LoG_image->FillBuffer(0);

	itk::ImageRegionIterator<LoGImageType> multiscale_LoG_image_iter(multiscale_LoG_image, multiscale_LoG_image->GetLargestPossibleRegion());

	//Get the maximum response across all scales for each pixel
	std::vector<LoGImageType::Pointer>::iterator LoG_vector_iter;
	for (LoG_vector_iter = LoG_vector.begin(); LoG_vector_iter != LoG_vector.end(); ++LoG_vector_iter)
	{
		LoGImageType::Pointer LoGimage = *LoG_vector_iter;
		
		itk::ImageRegionIterator<LoGImageType> LoGimage_iter(LoGimage, LoGimage->GetLargestPossibleRegion());
		
		multiscale_LoG_image_iter.GoToBegin();
		LoGimage_iter.GoToBegin();

		while (!multiscale_LoG_image_iter.IsAtEnd())
		{
			if (LoGimage_iter.Get() > multiscale_LoG_image_iter.Get())
				multiscale_LoG_image_iter.Set(LoGimage_iter.Get());

			++multiscale_LoG_image_iter; ++LoGimage_iter;
		}
		
	}

	std::ostringstream logimageFileNameStream;	

	logimageFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_LoG.mhd";

	//WriteLoGImage(logimageFileNameStream.str(), multiscale_LoG_image);
	return multiscale_LoG_image;
}