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

	ColorSegmentation *col_bin = new ColorSegmentation(reader->GetOutput());
	col_bin->SetTesting(true);
	col_bin->SetLightBackground(true);
	col_bin->SetIgnoreBackground(false);

	col_bin->TransformToRLI();
	col_bin->FindArchetypalColors();
	//col_bin->SetArchetypalColors(dh::RLI(160,154,78), dh::RLI(100,90,80), dh::RLI(128,122,182));
	col_bin->ComputeClassWeights();

	//col_bin->ComputeBinary(2,1);

	delete col_bin;
}
