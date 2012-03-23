#include <iostream>
#include <limits>
#include <QString>


#include "itkImageFileReader.h"

int main(int argc, char* argv[])
{
	typedef itk::Image< unsigned char, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	QString filename = "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu.tif";
	//QString filename = argv[1];

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename.toStdString());
	reader->SetNumberOfThreads(16);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}
}
