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
#include "mdlIntegratedSkeleton.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
    {
		std::cerr << "Usage: " << argv[0] << " inputImage outFilename <vectorMagnitude>" << std::endl;
		std::cerr << " Input Image Should be preprocessed before using this tool\n";
		std::cerr << " The output file should be .vtk - it is a list of skeleton points\n";
		std::cerr << " Parameters:\n";
		std::cerr << "  VectorMagnitude (0.1) Gradiant Magnitude Threshold for critical points\n";
		return EXIT_FAILURE;
    }

	std::string inputImageFilename = argv[1];
	std::string outputFilename = argv[2];

	double vectorMagnitude = .1;
	if(argc == 4)
	{
		vectorMagnitude = atof(argv[3]);
	}

	//*****************************************************************
	// IMAGE READER
	std::cout << "Reading Input...";

	typedef itk::ImageFileReader< mdl::ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inputImageFilename);
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
	reader = 0;

	std::cout << "done\n";
    
	//**********************Old Method **********************************
	//Integrated Skeleton to create skeleton points:
	std::cout << "Skeletonizing...";
	mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( img );
	skel->SetVectorMagnitude(vectorMagnitude);
	skel->SetDebug(false);
	skel->SetUseXiaoLiangMethod(false);
	skel->Update();
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;
	std::cout << "done\n";

	////****************************************************************
	//// TREE WRITER
	std::cout << "Saving...";
	mdl::vtkFileHandler * fhd2 = new mdl::vtkFileHandler();
	fhd2->SetNodes(&skeleton);
	fhd2->Write(outputFilename);
	delete fhd2;
	std::cout << "done\n";

	return EXIT_SUCCESS;
}
