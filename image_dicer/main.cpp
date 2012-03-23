#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"
#include "itkRegionOfInterestImageFilter.h"
#include <algorithm>

typedef itk::Image<unsigned char, 3> ImageType;

int main(int argc, char* argv[])
{
	std::ifstream centroids_file;
	centroids_file.open(argv[1], std::fstream::in);
	
	std::ofstream file_basenames("cell_list.txt");
	
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;	

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[2]);
	reader->Update();

	ImageType::Pointer image = reader->GetOutput();	


	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();


	while (centroids_file.good())
	{
		WriterType::Pointer writer = WriterType::New();
		double x_d, y_d, z_d;
		centroids_file >> x_d >> y_d >> z_d;

		int x, y, z;
		x = x_d;
		y = y_d;
		z = z_d;

		ImageType::SizeType montage_size = image->GetLargestPossibleRegion().GetSize();
		
		//std::cout << x << " " << y << " " << z << std::endl;
		ImageType::IndexType start;

		start[0] = std::max<long long>(0, x - 150);
		start[1] = std::max<long long>(0, y - 150);
		start[2] = std::max<long long>(0, z - 50);
		
		//std::cout << "start: " << start[0] << ", " << start[1] << ", " << start[2] << std::endl;
		
		ImageType::SizeType size;
		size[0] = std::min<long long>(300, montage_size[0] - start[0]);
		size[1] = std::min<long long>(300, montage_size[1] - start[1]);
		size[2] = std::min<long long>(100, montage_size[2] - start[2]);
	
		int local_centroid_coord_x = 150;
		int local_centroid_coord_y = 150; 
		int local_centroid_coord_z = 50;
		
		if (x-150 < 0)
			local_centroid_coord_x = 150 + (x - 150);
		if (y-150 < 0)
			local_centroid_coord_y = 150 + (y - 150);
		if (z-150 < 0)
			local_centroid_coord_z = 150 + (z - 150);		


		//std::cout << "size: " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;

		//Open file to write the centroids out
		std::ostringstream output_centroid_filename_stream;
		output_centroid_filename_stream << vul_file::strip_directory(vul_file::strip_extension(argv[2])) << "_" << x << "_" << y << "_" << z << ".txt";
	
		//Open file to write the local centroids out
		std::ostringstream output_local_centroid_filename_stream;
		output_local_centroid_filename_stream << vul_file::strip_directory(vul_file::strip_extension(argv[2])) << "_" << x << "_" << y << "_" << z << "_local_coord.txt";
		
	
		std::cout << "Writing " << output_centroid_filename_stream.str() << std::endl;
		std::ofstream centroid_file_coords(output_centroid_filename_stream.str().c_str());

		centroid_file_coords << x << " " << y << " " << z << std::endl;	//Write global coords
		centroid_file_coords << local_centroid_coord_x << " " << local_centroid_coord_y << " " << local_centroid_coord_z << std::endl;	//Write local coords
		centroid_file_coords.close();

		std::ofstream local_centroid_file_coords(output_local_centroid_filename_stream.str().c_str());
		local_centroid_file_coords << local_centroid_coord_x << " " << local_centroid_coord_y << " " << local_centroid_coord_z;
		local_centroid_file_coords.close();

		
		ImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		ROIfilter->SetRegionOfInterest(desiredRegion);
		ROIfilter->SetInput(image);

		std::ostringstream output_image_filename_stream;
		output_image_filename_stream << vul_file::strip_extension(argv[2]) << "_" << x << "_" << y << "_" << z << ".mhd";	
		std::cout << "Writing " << output_image_filename_stream.str() << std::endl;
		
		
		file_basenames << vul_file::strip_directory(vul_file::strip_extension(argv[2])) << "_" << x << "_" << y << "_" << z << std::endl;	//Write the basenames out to cell_list.txt		

		writer->SetFileName(std::string(output_image_filename_stream.str()));
		writer->SetInput(ROIfilter->GetOutput());
		writer->Update();		
	}

	file_basenames.close();
	centroids_file.close();
}


