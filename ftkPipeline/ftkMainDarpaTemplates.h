// ############################################################################################################################################################################
#ifndef _ftkMainDarpaTemplates_h_
#define _ftkMainDarpaTemplates_h_
// ############################################################################################################################################################################

// #include "ftkMainDarpaGlobalInclude.h"

template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
// 	printf("Writing %s ... \n",filename);
	std::cout << std::endl << "Writing ... " << filename;
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	itk::Size<3> inputImageSize = im->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2];
	return EXIT_SUCCESS;
}

template <typename TINPUT, typename TOUTPUT>
int writeImageRescaled(typename TINPUT::Pointer im, const char* filename)
{
// 	printf("Writing %s ... \n",filename);
	std::cout << std::endl << "Writing ... " << filename;
	typedef typename itk::RescaleIntensityImageFilter< TINPUT, TOUTPUT > rescalerfilter;
	typename rescalerfilter::Pointer rescaleFilter = rescalerfilter::New();
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(std::numeric_limits<typename TOUTPUT::PixelType>::max());
	rescaleFilter->SetInput(im);
	
	typedef typename itk::ImageFileWriter<TOUTPUT> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	try
	{
		writer->SetInput(rescaleFilter->GetOutput());
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	itk::Size<3> inputImageSize = im->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2];
	return EXIT_SUCCESS;
}


template <typename T>
typename T::Pointer readImage(const char* filename)
{
// 	printf("Reading %s ... \n",filename);
	std::cout << std::endl << "Reading ... " << filename;
	typedef typename itk::ImageFileReader<T> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	itk::Size<3> inputImageSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2];
	return reader->GetOutput();
}


template <typename T>
typename T::Pointer readImageRegion(const char* filename, typename T::RegionType region)
{
	printf("Reading %s ... \n",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
// 	try
// 	{
// 		reader->Update();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// // 		return EXIT_FAILURE;
// 	}
	typedef typename itk::ExtractImageFilter< T, T > ROIFilterType;
	typename ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	ROIfilter->SetExtractionRegion(region);
	ROIfilter->SetInput( reader->GetOutput() );
	ROIfilter->SetDirectionCollapseToIdentity();
	try
	{
		ROIfilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	typename T::Pointer out = ROIfilter->GetOutput();
	out->DisconnectPipeline();
	std::cout<<" done";
	return out;
}

template <typename TIN,typename TOUT>
typename TOUT::Pointer readerINcasterOUT( const char* name )
{
	typename TIN::Pointer inputImage = readImage<TIN>(name);
	typedef typename itk::CastImageFilter< TIN,TOUT > CasterFilterType;
	typename CasterFilterType::Pointer caster = CasterFilterType::New();
	caster->SetInput(inputImage);
	caster->Update();
	return caster->GetOutput();
}

#endif