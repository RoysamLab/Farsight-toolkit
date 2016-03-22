#include <iostream>
#include <limits>
#include <QString>


#include "itkImageFileReader.h"

int main(int argc, char* argv[])
{
	typedef itk::Image< unsigned char, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	QString filename = argv[1];

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename.toStdString());
	reader->SetNumberOfThreads(4);
	try
	{
		reader->Update();
		std::cout << "Read was successful." << std::endl;
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
		return 1;
	}

	return 0;
}
