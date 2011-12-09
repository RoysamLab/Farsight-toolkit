#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"
#include "itkRegionOfInterestImageFilter.h"

typedef itk::Image<unsigned char, 3> ImageType;

int main(int argc, char* argv[])
{
	std::ifstream centroids_file;
	centroids_file.open(argv[1], std::fstream::in);

	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;	

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[2]);
	reader->Update();

	ImageType::Pointer image = reader->GetOutput();	

	WriterType::Pointer writer = WriterType::New();

	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();


	while (centroids_file.good())
	{
		double x_d, y_d, z_d;
		centroids_file >> x_d >> y_d >> z_d;

		int x, y, z;
		x = x_d;
		y = y_d;
		z = z_d;
		
		//std::cout << x << " " << y << " " << z << std::endl;
		ImageType::IndexType start;

		start[0] = y - 150;
		start[1] = x - 150;
		start[2] = z - 50;

		ImageType::SizeType size;
		size[0] = 300;
		size[1] = 300;
		size[2] = 100;

		//Open file to write the centroids out
		std::ostringstream output_centroid_filename_stream;
		output_centroid_filename_stream << vul_file::strip_extension(argv[2]) << "_" << x << "_" << y << "_" << z << ".txt";
		std::ofstream centroid_file_coords(output_centroid_filename_stream.str().c_str());

		centroid_file_coords << x << "_" << y << "_" << z << std::endl;	//Write global coords
		centroid_file_coords << 150 << "_" << 150 << "_" << "50" << std::endl;	//Write lobal coords
		centroid_file_coords.close();
		
		ImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		ROIfilter->SetRegionOfInterest(desiredRegion);
		ROIfilter->SetInput(image);

		std::ostringstream output_image_filename_stream;
		output_image_filename_stream << vul_file::strip_extension(argv[2]) << "_" << x << "_" << y << "_" << z << ".mhd";	

		writer->SetFileName(std::string(output_image_filename_stream.str()));
		writer->SetInput(ROIfilter->GetOutput());
		writer->Update();		

	}

	centroids_file.close();
}