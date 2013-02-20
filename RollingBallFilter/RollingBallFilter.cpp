#include "RollingBallFilter.h"

RollingBallFilter::RollingBallFilter(char* fileName, float radius)
{
	//Make reader
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fileName);
	
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Reader exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	
	input_img = reader->GetOutput();
	this->radius = radius;
}

RollingBallFilter::RollingBallFilter(ImageType::Pointer input_img, float radius)
{
	this->input_img = input_img;
	this->radius = radius;
}

void RollingBallFilter::RunFilter()
{
	//Structuring element
	typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 2> StructuringElementType;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius);
	structuringElement.CreateStructuringElement();

	//SliceBySliceFilter
	typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceBySliceFilterType;
	SliceBySliceFilterType::Pointer SliceBySliceFilter = SliceBySliceFilterType::New();

	//Erode Filter
	typedef itk::GrayscaleErodeImageFilter< SliceBySliceFilterType::InternalInputImageType, SliceBySliceFilterType::InternalOutputImageType, StructuringElementType > GrayscaleErodeImageFilter; //This filter takes the minimum of all the neighboring values and puts it into the pixel location
	GrayscaleErodeImageFilter::Pointer ErodeFilter = GrayscaleErodeImageFilter::New(); 
	ErodeFilter->SetKernel(structuringElement);
	SliceBySliceFilter->SetInput(input_img);
	SliceBySliceFilter->SetFilter(ErodeFilter);

	//Subtract Filter
	typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractFilterType;
	SubtractFilterType::Pointer Subtractfilter = SubtractFilterType::New();
	Subtractfilter->SetInput1(input_img);
	Subtractfilter->SetInput2(SliceBySliceFilter->GetOutput()); //Note that the SliceBySliceFilter output is a 3D image again
	
	try
	{
		Subtractfilter->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Filter exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	output_img = Subtractfilter->GetOutput();
}

RollingBallFilter::ImageType::Pointer RollingBallFilter::GetOutput()
{
	return output_img;
}

void RollingBallFilter::WriteOutput(char* fileName)
{
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(output_img);

	std::cout << "Writing file: " << fileName << std::endl;
	writer->SetFileName(fileName);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Writer exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}

