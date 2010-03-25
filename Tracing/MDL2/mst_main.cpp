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
#include "mdlMST.h"
#include "mdlUtils.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
    {
		std::cerr << "Usage: " << argv[0] << " inputImage inputSkeleton < >" << std::endl;
		return EXIT_FAILURE;
    }
	
	//*****************************************************************
	// IMAGE READER
	std::string inputImageFilename = argv[1];

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

	//****************************************************************
	// SKELETON POINTS READER
	std::string inputSkelFilename = argv[2];

	std::vector<mdl::fPoint3D> skeleton;
	mdl::vtkFileHandler * fhdl = new mdl::vtkFileHandler();
	fhdl->SetNodes(&skeleton);
	fhdl->Read(inputSkelFilename);
	delete fhdl;

	//******************************************************************
	// MST
	std::cerr << "Starting MST\n";

	//Minimum spanning tree to create nodes and backbone node pairs (lines):
	mdl::MST *mst = new mdl::MST( img );
	mst->SetDebug(true);
	mst->SetUseVoxelRounding(true);
	mst->SetEdgeRange(25);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	// can choose different weight
	//mst->CreateGraphAndMST(3);
	mst->CreateGraphAndMST(1);
	mst->ErodeAndDialateNodeDegree(0);

	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();

	delete mst;

	std::cerr << "Saving\n";

	mdl::vtkFileHandler * fhd2 = new mdl::vtkFileHandler();
	fhd2->SetNodes(&nodes);
	fhd2->SetLines(&bbpairs);
	fhd2->Write("Backbone.vtk");
	delete fhd2;

	std::cerr << "DONE\n";
  
	std::cerr << "PRESS ENTER TO EXIT\n";
	getchar();

	//******************************************************************
	//******************************************************************

	return EXIT_SUCCESS;
}
