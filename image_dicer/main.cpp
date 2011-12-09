#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"

typedef itk::Image<unsigned char, 3> ImageType;

int main(int argc, char* argv[])
{
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;	
	
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	//reader->Update();
	
	//ImageType::Pointer image = reader->GetOutput();	
	
	WriterType::Pointer writer = WriterType::New();

	ImageType::RegionType region;
	ImageType::IndexType start;
	
	int x = atoi(argv[2]);
	int y = atoi(argv[3]);
	int z = atoi(argv[4]);
	
	start[0] = y - 150;
	start[1] = x - 150;
	start[2] = z - 50;
		
	ImageType::SizeType size;
	size[0] = 300;
	size[1] = 300;
	size[2] = 100;

	region.SetSize(size);
	region.SetIndex(start);

	std::ostringstream output_filename_stream;

	output_filename_stream << vul_file::strip_extension(argv[2]) << "_" << x << "_" << y << "_" << z << ".TIF";	

	writer->SetFileName(std::string(output_filename_stream.str()));
	writer->SetInput(reader->GetOutput());
	writer->Update();		

}


