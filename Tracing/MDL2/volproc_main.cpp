/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#include "mdlVolumeProcess.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 2)
    {
		std::cerr << "Usage: " << argv[0] << " <inputfilename> <outputfilename> " << std::endl;
		return EXIT_FAILURE;
    }
	std::string inputFileName = argv[1];
	std::string outputFileName = "out.tif";
	if(argc == 3)
		outputFileName = argv[2];

	//*****************************************************************
	// IMAGE READER
	typedef itk::ImageFileReader< mdl::ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inputFileName);
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "READER FAILED: " << err << std::endl ;
		return EXIT_FAILURE;
	}
	mdl::ImageType::Pointer img = reader->GetOutput();
	//******************************************************************
	//******************************************************************
	
	//******************************************************************
	// PROCESSING
	mdl::VolumeProcess *volProc = new mdl::VolumeProcess();
	volProc->SetInput(img);
	volProc->SetDebug(true);
	//volProc->RescaleIntensities(0, 255);
	//volProc->RunCAD(); //Curvature Anisotropic Diffusion
	//volProc->RunOtsuDenoising();
	//volProc->DialateImage(1);
	volProc->MaskUsingGraphCuts();
	//volProc->MaskSmallConnComp(50);

	//******************************************************************
	//******************************************************************
	
	//******************************************************************
	// IMAGE WRITER
	typedef itk::ImageFileWriter< mdl::ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputFileName);
	writer->SetInput( volProc->GetOutput() );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "READER FAILED: " << err << std::endl ;
		return EXIT_FAILURE;
	}
	//******************************************************************
	//******************************************************************

	return EXIT_SUCCESS;
}
