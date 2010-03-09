#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "ColorSegmentation.h"

int main(int argc, char* argv[]) 
{
	if( argc < 2 )
	{
		std::cerr<<argv[0]<<" InputImage"<<std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::ImageFileReader< RGBImageType >  ReaderType;
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

	ColorSegmentation *col_bin = new ColorSegmentation(reader->GetOutput(),1,2,1);
	col_bin->SetTesting();
	col_bin->RunInitialBinarization();
	col_bin->FindArchetypalColors();
	col_bin->ComputeClassWeights();

	delete col_bin;
}
