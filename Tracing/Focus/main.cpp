#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "Focus.h"

int main(int argc, char* argv[]) 
{
	if( argc < 2 )
	{
		std::cerr<<argv[0]<<" InputImage"<<std::endl;
		return EXIT_FAILURE;
	}

	typedef RGBImageType InputImageType;
	typedef RGBImageType2D OutputImageType;

	typedef itk::ImageFileReader< InputImageType >  ReaderType;
	//typedef itk::ImageFileReader< UCharImageType >  ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );
	std::cout<<argv[1]<<std::endl;
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	InputImageType::Pointer in = reader->GetOutput();

	Focus *f = new Focus( in );
	f->SetRadius(10);
	//std::vector<float> vars = f->FindVariance(1.93,4.57,300);
	f->MakeVarianceImage();

	//OutputImageType::Pointer p = f->MakeProjection();
	OutputImageType::Pointer p = f->MakeProjectionColor();

	delete f;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "best.tif" );
	writer->SetInput( p );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	std::cerr << "PRESS ENTER TO EXIT\n";
	getchar();
}
