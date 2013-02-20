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

	typedef UCharImageType InputImageType;
	typedef UCharImageType2D OutputImageType;

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
	f->SetRadius(atoi(argv[2]));
	//std::vector<float> vars = f->FindVariance(1.93,4.57,300);
	f->MakeVarianceImage();

	//OutputImageType::Pointer p = f->MakeProjection();
	typedef std::pair<UCharImageType2D::Pointer, UShortImageType2D::Pointer> PairType;
	PairType p = f->MakeProjection2();

	delete f;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	std::string out = argv[1];
	out = out.substr(0,out.size()-4) + "_BF.tif";
	writer->SetFileName( out.c_str());
	writer->SetInput( p.first );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	std::string outdepth = argv[1];
	outdepth = outdepth.substr(0,outdepth.size()-4) + "_DM.tif";
	typedef itk::ImageFileWriter< UShortImageType2D >  WriterType2;
	WriterType2::Pointer writer2 = WriterType2::New();
	writer2->SetFileName( outdepth.c_str());
	writer2->SetInput( p.second );
	try
	{
		//writer2->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//std::cerr << "PRESS ENTER TO EXIT\n";
	//getchar();
}
