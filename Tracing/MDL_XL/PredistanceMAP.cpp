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
#include "mdlUtils.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 2)
    {
		std::cerr << "Usage: " << argv[0] << " inputfilename" << std::endl;
		return EXIT_FAILURE;
    }

	std::string inputFileName = argv[1];
	int edgeRange = 10;
	int morphStrength =20;
	float VectorMagnitude =0.10;
    float threshold = 20;
	int   boolIntensity = 0;

	if(argc < 8)
	{
		threshold = atoi(argv[2]);  
		VectorMagnitude = atoi(argv[3]);
		edgeRange = atoi(argv[4]);
		morphStrength = atoi(argv[5]);   
	}

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
	std::cerr << "Volume Processing\n";
    //***************************New Test Method ***********************	
    typedef itk::ImageFileWriter< mdl::ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
   	
	mdl::VolumeProcess *volProc = new mdl::VolumeProcess();
	volProc->SetInput(img);
	volProc->SetDebug(true);

   
    mdl::ImageType::Pointer PreImage = volProc->GetOutput();
	
	writer->SetInput(PreImage);
	writer->SetFileName("InputImage.tif");
    writer->Update();
	 volProc->RunIntensityNormalize();
    //volProc->RunFillingZeroOnBouandary(4,4,4);
   
	/*
	volProc->RunBinaryForDistanceMapUsingManualThreshold(threshold);
	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("BinaryResults.tif");
    writer->Update();

	volProc->RunBinayMedianHoleFilling(1,1,1);
    writer->SetInput(volProc->GetOutput());
	writer->SetFileName("RunBinayMedianHoleFillingResults.tif");
    writer->Update();
    

	volProc->RunDanielssonDistanceMap();
	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("DanielssonDistanceMapResults.tif");
    writer->Update();
    
	volProc->RunGaussianSmoothing(7,7,3,0.1);
    */
	std::string outputFileName = inputFileName;
	size_t found = outputFileName.find_last_of(".");
	outputFileName.insert(found,"_CorrectedDTMAP");
	found = outputFileName.find_last_of(".");
	outputFileName.replace(found+1,3,"tif");


	writer->SetInput(volProc->GetOutput());
	writer->SetFileName(outputFileName);
    writer->Update();
   
	delete volProc;
  
	return EXIT_SUCCESS;
}
