#include "LoG.h"

#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include "ftkTimeStampOverflowSafeUpdate.h"

LoG::LoGImageType::Pointer LoG::RunLoG(const ImageType::Pointer & image, float scale)
{
	
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType , LoGImageType> LoGFilterType;
    LoGFilterType::Pointer LoGFilter = LoGFilterType::New();    
	LoGFilter->SetInput( image );
	LoGFilter->SetNormalizeAcrossScale(true);
	LoGFilter->SetSigma( scale );

	try
	{
		//std::cerr << "Modified Time: " << LoGFilter->GetMTime() << std::endl;
		ftk::TimeStampOverflowSafeUpdate( LoGFilter.GetPointer() );
		//LoGFilter->Update();
		
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "LoGFilter Exception: " << err << std::endl;
		std::cerr << "LoGFilter: " << LoGFilter << std::endl;
		std::cerr << "image: " << image << std::endl;
		std::cerr << "Requested Region: " << image->GetRequestedRegion() << std::endl;
		std::cerr << "Largest Possible Region: " << image->GetLargestPossibleRegion() << std::endl;
		std::cerr << "Buffered Region: " << image->GetBufferedRegion() << std::endl;
		throw; //Rethrowing exception, since only RunMultiScaleLoG knows the location of the cell
	}
	LoGImageType::Pointer LoG_image = LoGFilter->GetOutput();

	//LoG_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...

	//Scale to (-1, 1) range and invert   
	typedef itk::ShiftScaleImageFilter<LoGImageType, LoGImageType> InvertFilterType;
	InvertFilterType::Pointer invertFilter = InvertFilterType::New();
	invertFilter->SetScale(-1.0 / std::numeric_limits<ImageType::PixelType>::max());
	invertFilter->SetInput(LoG_image);

	try
	{
		ftk::TimeStampOverflowSafeUpdate( invertFilter.GetPointer() );
		//invertFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "invertFilter Exception: " << err << std::endl;
	}

	LoGImageType::Pointer inverted_LoG_image = invertFilter->GetOutput();

	inverted_LoG_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...

    //std::cerr << inverted_LoG_image << std::endl;
	return inverted_LoG_image;
}

void LoG::WriteLoGImage(const std::string & filename, const LoGImageType::Pointer & image)
{
	typedef itk::ImageFileWriter< LoGImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(filename);
	try
	{
		std::cout << "Writing LoG image" << std::endl;
		ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );		
		//writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

LoG::LoGImageType::Pointer LoG::RunMultiScaleLoG(const Cell & cell)
{
	std::vector<LoGImageType::Pointer> LoG_vector;

	//Calculate all the LoG scales
	for (float scale = 1.0; scale <= 4.0; scale+=0.1)
	{
		LoGImageType::Pointer LoGimage;

		try
		{
			LoGimage = RunLoG(cell.isotropic_image, scale);
		}
		catch (itk::ExceptionObject & err)
		{
			ImageType::PointType origin = cell.isotropic_image->GetOrigin();
			ImageType::SizeType size = cell.isotropic_image->GetLargestPossibleRegion().GetSize();
			
			std::cerr << "RunMultiScaleLoG exception: " << std::endl;
			std::cerr << "For cell: " << cell.getX() << ", " << cell.getY() << ", " << cell.getZ() << " at scale: " << scale << " Origin: " << origin << " Size: " << size << std::endl;
			std::cerr << err << std::endl;
		}
        
		LoG_vector.push_back(LoGimage);
	}

	LoGImageType::SizeType size = cell.isotropic_image->GetLargestPossibleRegion().GetSize();

	//Make a new image to store the multiscale LoG image	
	LoGImageType::Pointer multiscale_LoG_image = LoGImageType::New();
	LoGImageType::IndexType start;
	start.Fill(0);

	LoGImageType::RegionType region(start, size);
	multiscale_LoG_image->SetRegions(region);
	multiscale_LoG_image->Allocate();
	multiscale_LoG_image->FillBuffer(0);
	multiscale_LoG_image->SetSpacing(cell.isotropic_image->GetSpacing());

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

			++multiscale_LoG_image_iter; 
			++LoGimage_iter;
		}
		
	}

	return multiscale_LoG_image;
}
