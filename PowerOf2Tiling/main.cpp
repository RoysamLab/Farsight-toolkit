#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include <string>
#include "vul/vul_file.h"
#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " " << "image_to_be_tiled" << std::endl;
		return -1;
	}
	

	char *arg1 = argv[1];

	//Make the file names
	vcl_string fileName = vul_file::strip_directory(argv[1]);
	vcl_string fileName_no_ext = vul_file::strip_extension(fileName);	

	//Some typedefs to make things easier to read later
	#define ImageDimensions 3
	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, ImageDimensions >	ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	
	//Make reader
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);	
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Reader exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return -1;
	}
	

	//SizeType typedef
	typedef ImageType::SizeType SizeType;
	
	//Get the size and print it out to the console
	SizeType input_image_size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	std::cout << "Image Size: " << input_image_size[0] << "x" << input_image_size[1] << "x" << input_image_size[2] << std::endl;
	
	//Padding the image to the nearest size of 2
	typedef itk::ConstantPadImageFilter <ImageType, ImageType> ConstantPadImageFilterType;
	ImageType::SizeType upperExtendRegion;
	upperExtendRegion[0] = 1024 - (input_image_size[0] % 1024);
	upperExtendRegion[1] = 1024 - (input_image_size[1] % 1024);
	upperExtendRegion[2] = 0;

	ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
	padFilter->SetInput(reader->GetOutput());
	padFilter->SetPadUpperBound(upperExtendRegion);

	try
	{
		padFilter->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Reader exception caught !" << std::endl;
		std::cerr << err << std::endl;
		return -1;
	}

	SizeType padded_image_size = padFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
	std::cout << "Image Size: " << padded_image_size[0] << "x" << padded_image_size[1] << "x" << padded_image_size[2] << std::endl;

	//ROIfilter
	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();

	//Make the writer
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	//Write the file names to a txt file
	std::ofstream list_of_tile_filenames;
	list_of_tile_filenames.open("list_of_tile_filenames.txt");

	for (int k = 0; k < input_image_size[0]; k+=1024)
	{
		for (int l = 0; l < input_image_size[1]; l+=1024)
		{

			//Set the starting position of where to extract
			ImageType::IndexType start;
			start[0] = k;
			start[1] = l;
			start[2] = 0;

			//Set the size to extract
			ImageType::SizeType size;
			size[0] = 1024;
			size[1] = 1024;
			size[2] = input_image_size[2];

			//Make the region
			ImageType::RegionType desiredRegion;
			desiredRegion.SetSize(size);
			desiredRegion.SetIndex(start);

			ROIfilter->SetRegionOfInterest(desiredRegion);
			ROIfilter->SetInput(padFilter->GetOutput());
			
			//Set the writer input
			writer->SetInput(ROIfilter->GetOutput());

			std::stringstream output_name_stream;
			output_name_stream << fileName << "_Tile_" << k << "_" << l << ".tif";
			
			std::cout << "Writing file: " << output_name_stream.str() << std::endl;
			writer->SetFileName(output_name_stream.str().c_str());
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
			list_of_tile_filenames << output_name_stream.str() << "\n";
		}
	}
	list_of_tile_filenames.close();
	
	return 0;
}
