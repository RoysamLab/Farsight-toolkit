//#include <boost/cstdint.hpp>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "ftkLaplacianOfGaussian3D.h"


//Function prototypes
double*** AllocateMemoryForImage(size_t x_size, size_t y_size, size_t z_size);

int main(int argc, char *argv[])
{
	const unsigned int Dimension = 3;
	typedef itk::Image<unsigned char, Dimension> InputImageType;
	typedef itk::Image<double, Dimension> LoGImageType;
	typedef itk::Image<unsigned short, Dimension> OutputImageType;
	typedef itk::ImageFileReader<InputImageType> FileReaderType;
	typedef itk::ImageRegionConstIterator<InputImageType> InputImageIteratorType;
	typedef itk::ImageRegionIterator<LoGImageType> LogImageIteratorType;

	FileReaderType::Pointer reader = FileReaderType::New();
	reader->SetFileName("test.tif");
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "Error in reader: " << err << std::endl;
		return -1;
	}

	InputImageType::Pointer image = reader->GetOutput();

	size_t image_x_size = image->GetLargestPossibleRegion().GetSize()[0];
	size_t image_y_size = image->GetLargestPossibleRegion().GetSize()[1];
	size_t image_z_size = image->GetLargestPossibleRegion().GetSize()[2];

	std::cout << "Image size: " << image_x_size << "x" << image_y_size << "x" << image_z_size << std::endl;

	std::cout << "Allocating memory for image...";
	double ***image_matrix = AllocateMemoryForImage(image_x_size, image_y_size, image_z_size);
	std::cout << "Done" << std::endl;

	InputImageIteratorType input_image_iterator(image, image->GetRequestedRegion());

	std::cout << "Iterating over itkImage to convert it to a 3D C array...";
	size_t index;
	for (input_image_iterator.GoToBegin(), index = 0; !input_image_iterator.IsAtEnd(); ++input_image_iterator, ++index)
	{
		size_t k = index % image_x_size;
		size_t l = (index % (image_x_size * image_y_size)) / image_x_size;
		size_t m = index / (image_x_size * image_y_size);

		image_matrix[k][l][m] = input_image_iterator.Value();
	}
	std::cout << "done" << std::endl;
	
	double scale_X = 1;
	double scale_Y = 1;
	double scale_Z = 1;

	/* 3x3 matrix testing */
	image_x_size = 3;
	image_y_size = 3;
	image_z_size = 1;
	image_matrix = AllocateMemoryForImage(image_x_size, image_y_size, image_z_size);

	int matrixIndex = 0;
	for (int m = 0; m < image_z_size; m++)
		for (int l = 0; l < image_y_size; l++)
			for (int k = 0; k < image_x_size; k++)
				image_matrix[k][l][m] = matrixIndex++;
	
	ftkLaplacianOfGaussian3D<double> *LoGFilter = new ftkLaplacianOfGaussian3D<double>(image_matrix, scale_X, scale_Y, scale_Z, image_x_size, image_y_size, image_z_size);

	LoGFilter->RunFilter();
	double ***LoGImage = LoGFilter->GetOutput();

	LoGImageType::RegionType region;
	LoGImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;

	LoGImageType::SizeType size;
	size[0] = image_x_size;
	size[1] = image_y_size;
	size[2] = image_z_size;

	region.SetSize(size);
	region.SetIndex(start);

	LoGImageType::Pointer log_image = LoGImageType::New();
	log_image->SetRegions(region);
	log_image->Allocate();

	LogImageIteratorType log_image_iterator(log_image, log_image->GetRequestedRegion());
	
	std::cout << "Iterating over C array to convert it back to an ITK image...";

	for (log_image_iterator.GoToBegin(), index = 0; !log_image_iterator.IsAtEnd(); ++log_image_iterator, ++index)
	{
		size_t k = index % image_x_size;
		size_t l = (index % (image_x_size * image_y_size)) / image_x_size;
		size_t m = index / (image_x_size * image_y_size);

		log_image_iterator.Value() = LoGImage[k][l][m]; 
	}
	std::cout << "done" << std::endl;

	typedef itk::RescaleIntensityImageFilter<LoGImageType, OutputImageType> RescaleIntensityImageFilterType;
	RescaleIntensityImageFilterType::Pointer rescale_filter = RescaleIntensityImageFilterType::New();
	rescale_filter->SetInput(log_image);
	try
	{
		rescale_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "Error in RescaleIntensityFilter: " << err << std::endl;
		return -1;
	}

	typedef itk::ImageFileWriter<OutputImageType> OutputImageWriterType;
	OutputImageWriterType::Pointer writer = OutputImageWriterType::New();
	writer->SetFileName("LoGtest.tif");
	writer->SetInput(rescale_filter->GetOutput());
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "Error in writer: " << err << std::endl;
		return -1;
	}
	return 0;
}

double*** AllocateMemoryForImage(size_t x_size, size_t y_size, size_t z_size)
{
	double*** image = new double**[x_size];
	for (size_t k = 0; k < x_size; k++)
	{
		image[k] = new double*[y_size];
		for (size_t l = 0; l < y_size; l++)
		{
			image[k][l] = new double[z_size];
			for (size_t m = 0; m < z_size; m++)
				image[k][l][m] = 0;
		}
	}

	return image;
}