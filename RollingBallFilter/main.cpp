#include "itkSubtractImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSliceBySliceImageFilter.h"

#include "itkStreamingImageFilter.h"

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " input_image output_image" << std::endl;
		return -1;
	}

	//Some typedefs to make things easier to read later
	const int ImageDimensions = 3;
	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, ImageDimensions >	ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	//Make reader
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);	
	//std::cout << "Reading " << argv[1] << std::endl;
	//try
	//{
	//	reader->Update();
	//}
	//catch (itk::ExceptionObject & err)
	//{
	//	std::cerr << "Reader exception caught !" << std::endl;
	//	std::cerr << err << std::endl;
	//	return -1;
	//}

	//Slice-by-Slice Image Filter
	typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceBySliceFilterType;
	SliceBySliceFilterType::Pointer SliceBySliceFilter = SliceBySliceFilterType::New();

	//Mean Filter
	typedef itk::MeanImageFilter< SliceBySliceFilterType::InternalInputImageType, SliceBySliceFilterType::InternalOutputImageType > MeanFilterType;
	MeanFilterType::Pointer Meanfilter = MeanFilterType::New(); 
	Meanfilter->SetRadius(50);
	SliceBySliceFilter->SetInput(reader->GetOutput());
	SliceBySliceFilter->SetFilter(Meanfilter);
	//std::cout << "Running Mean Filtering" << std::endl;
	//try
	//{
	//	SliceBySliceFilter->Update();
	//}
	//catch (itk::ExceptionObject & err)
	//{
	//	std::cerr << "Meanfilter exception caught !" << std::endl;
	//	std::cerr << err << std::endl;
	//	return -1;
	//}
	
	//Subtract Filter
	typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractFilterType;
	SubtractFilterType::Pointer Subtractfilter = SubtractFilterType::New();
	Subtractfilter->SetInput1(reader->GetOutput());
	Subtractfilter->SetInput2(SliceBySliceFilter->GetOutput()); //Note that the SliceBySliceFilter output is a 3D image again
	
	/*std::cout << "Running Subtraction Filtering" << std::endl;
	try
	{
		Subtractfilter->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Subtractfilter exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return -1;
	}*/


	//Make the writer
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(Subtractfilter->GetOutput());

	std::cout << "Writing file: " << argv[2] << std::endl;
	writer->SetFileName(argv[2]);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Writer exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return -1;
	}
	


}