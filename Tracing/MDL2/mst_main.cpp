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
	if(argc < 4)
    {
		std::cerr << "Usage: " << argv[0] << " inputImage inputSkeleton outFilename < edgeRange > < morphStrength >" << std::endl;
		std::cerr << " Input Image is used in determining edge weights for MST\n";
		std::cerr << " Input Skeleton must be a .vtk file containing skeleton points\n";
		std::cerr << " The output file should be .vtk - it is nodes and line-pairs from MST\n";
		std::cerr << " Parameters:\n";
		std::cerr << "  edgeRange (25) Maximum Edge Length (euclidean dist - voxels)\n";
		std::cerr << "  morphStrength (0) Number of iterations in graph erosion/dialation\n";
		return EXIT_FAILURE;
    }

	std::string inputImageFilename = argv[1];
	std::string inputSkelFilename = argv[2];
	std::string outFilename = argv[3];

	int edgeRange = 25;
	int morphStrength = 0;
	if(argc == 5)
	{
		edgeRange = atoi(argv[4]);
	}
	if(argc == 6)
	{
		morphStrength = atoi(argv[5]);
	}


	//*****************************************************************
	// IMAGE READER
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
	mst->SetDebug(false);
	mst->SetUseVoxelRounding(false);
	mst->SetEdgeRange(edgeRange);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	// can choose different weight
	//mst->CreateGraphAndMST(3);
	mst->CreateGraphAndMST(1);
	mst->ErodeAndDialateNodeDegree(morphStrength);

	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();

	delete mst;

	std::cerr << "Saving\n";

	//****************************************************************
	// TREE WRITER
	mdl::vtkFileHandler * fhd2 = new mdl::vtkFileHandler();
	fhd2->SetNodes(&nodes);
	fhd2->SetLines(&bbpairs);
	fhd2->Write(outFilename);
	delete fhd2;

	std::cerr << "DONE\n";
  
	//std::cerr << "PRESS ENTER TO EXIT\n";
	//getchar();

	//******************************************************************
	//******************************************************************

	return EXIT_SUCCESS;
}
