#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer()
{

}

void MicrogliaRegionTracer::LoadImage(ImageType::Pointer image)
{
	this->image = image;
}

void MicrogliaRegionTracer::LoadImage(std::string filename)
{
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}

	LoadImage(reader->GetOutput());
}

void MicrogliaRegionTracer::LoadSeedPoints(std::string filename)
{
	std::ifstream seed_point_file;
	seed_point_file.open(filename.c_str());

	while (!seed_point_file.eof())
	{
		size_t seedX, seedY, seedZ;
		seed_point_file >> seedX >> seedY >> seedZ;

		//std::cout << "Reading in seed: (" << seedX << ", " << seedY << ", " << seedZ << ")" << std::endl;
		seeds.push_back(new Seed(seedX, seedY, seedZ));
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, ImageType::Pointer image)
{
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::WriteLoGImage(std::string filename, LoGImageType::Pointer image)
{
	typedef itk::ImageFileWriter< LoGImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		std::cout << "Writing LoG image" << std::endl;
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::RunLoG()
{
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType , LoGImageType> LoGFilterType;
	LoGFilterType::Pointer LoGFilter = LoGFilterType::New();
	LoGFilter->SetInput( image );
	LoGFilter->SetNormalizeAcrossScale(true);
	
	float scales[1] = {1};
	for (int k = 0; k < 1; k++)
	{
		float scale = scales[k];
		LoGFilter->SetSigma( scale );

		try
		{
			LoGFilter->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cerr << "gauss Exception: " << err << std::endl;
		}

		std::stringstream scale_stream;
		scale_stream << scale;
		
		WriteLoGImage("D:/farsight_images/MicrogliaRegionTracer/LoGImage" + scale_stream.str() + ".mhd", LoGFilter->GetOutput());
	}
}