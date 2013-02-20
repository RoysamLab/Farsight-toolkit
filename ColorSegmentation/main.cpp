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
	col_bin->SetTesting(false);
	col_bin->SetLightBackground(true);
	col_bin->SetIgnoreBackground(true);
	col_bin->SetGenerateProjections(false);

	//col_bin->InvertInput();
	//col_bin->SmoothInput();
	//col_bin->ComputeBinary(3,2);
	//col_bin->MaskBackgroundFromInput();

	col_bin->TransformToRLI();
	col_bin->FindArchetypalColors();
	//col_bin->SetArchetypalColors(dh::_RGB(255,0,0), dh::_RGB(0,0,255), dh::_RGB(175,175,175));
	col_bin->ComputeClassWeights();
	//col_bin->VoteBasedOnWeights();

	delete col_bin;

	std::cerr << "PRESS ENTER TO EXIT\n";
	getchar();
}
