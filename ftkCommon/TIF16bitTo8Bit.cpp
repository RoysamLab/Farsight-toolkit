#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "Usage: "
			<< "<input filename> "
			<< "<output filename> "
			<< std::endl;
		return 1;
	}
	
	typedef unsigned short InputPixelType;
	typedef unsigned char OutputPixelType;
	
	const unsigned Dimension = 3;

	typedef itk::Image<InputPixelType, Dimension> InputImageType;
	typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

	typedef itk::ImageFileReader<InputImageType> InputReaderType;
	typedef itk::ImageFileWriter<OutputImageType> OutputWriterType;

	InputReaderType::Pointer reader = InputReaderType::New();
	OutputWriterType::Pointer writer = OutputWriterType::New();

	const char* outputFilename = argv[2];
	const char* inputFilename = argv[1];

	reader->SetFileName(inputFilename);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Error in reader: " << err << std::endl;
		return -1;
	}
	
	typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();

	rescaleFilter->SetInput(reader->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	
	try
	{
		rescaleFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Error in rescaleFilter: " << err << std::endl;
		return -1;
	}

	writer->SetInput(rescaleFilter->GetOutput());
	
	writer->SetFileName(outputFilename);
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Error in writer: " << err << std::endl;
		return -1;
	}
	writer->Update();	
}