#include <iostream>
#include <limits>

#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkMultiThreader.h"
#ifdef _OPENMP
    #include "omp.h"
#endif
#include "ftkTimeStampOverflowSafeUpdate.h"

int main(int argc, char* argv[])
{
	// Since we are using OpenMP to run a filter in each thread, we do not
	// want the individual filters to also spawn threads -- this will result
	// in too many threads.
	itk::MultiThreader::SetGlobalDefaultNumberOfThreads( 1 );

	unsigned long modified_count = std::numeric_limits< unsigned long >::max();	//Maximum value that an unsigned long can take

	typedef itk::Image< unsigned char, 3> ImageType;
	typedef itk::Image< float, 3> FloatImageType;

	typedef itk::ImageFileReader< ImageType > ImageReaderType;
	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}
	std::cout << "Done reading in image" << std::endl;

	ImageType::Pointer image = reader->GetOutput();
	image->DisconnectPipeline();

	std::cout << "Starting to increment the modified timer" << std::endl;
	while (image->GetMTime() < modified_count - 10000)
	{
		image->Modified();

		//Just so we have some visual output of incrementing the modified timer
		if (image->GetMTime() % 100000000 == 0)
			std::cout << "Modified time: " << image->GetMTime() << std::endl;
	}
	std::cout << "Done incrementing the modified time" << std::endl;

	//And now we execute across the rollover
	bool error_happened = false;
	#pragma omp parallel for
	for (int scale = 2; scale < 50; scale += 1)
	{
		typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType , FloatImageType> LoGFilterType;
		LoGFilterType::Pointer LoGFilter = LoGFilterType::New();
		std::cout << "LoGFilter->GetMTime(): " << LoGFilter->GetMTime() << std::endl;
		std::cout << "image->GetMTime(): " << image->GetMTime() << std::endl;
		LoGFilter->SetInput( image );
		LoGFilter->SetNormalizeAcrossScale(true);
		LoGFilter->SetSigma( scale );

		try
		{
			//LoGFilter->Update();
			ftk::TimeStampOverflowSafeUpdate( LoGFilter.GetPointer() );
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "LoGFilter Exception: " << err << std::endl;
			std::cout << "Image Buffered: " << image->GetBufferedRegion() << std::endl;
			std::cout << "LOG Buffered: " << LoGFilter->GetOutput()->GetBufferedRegion() << std::endl;
			error_happened = true;
		}
	}

	if (error_happened)
		{
		return EXIT_FAILURE;
		}
}
