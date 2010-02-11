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
#include "mdlIntegratedSkeleton.h"
#include "mdlMST.h"
#include "mdlBSplineFitting.h"
#include "mdlUtils.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
	if(argc < 2)
    {
		std::cerr << "Usage: " << argv[0] << " <inputfilename>" << std::endl;
		return EXIT_FAILURE;
    }
	std::string inputFileName = argv[1];

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
	//volProc->MaskUsingGraphCuts();
	//volProc->MaskSmallConnComp(50);
	mdl::ImageType::Pointer clean_img = volProc->GetOutput();
	//volProc->RunDistanceTransform();
	//mdl::ImageType::Pointer DT_img = volProc->GetOutput();
	delete volProc;

	/*
	//Integrated Skeleton to create skeleton points:
	mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( clean_img );
	skel->SetVectorMagnitude(.05);
	skel->SetDebug(true);
	skel->SetUseXiaoLiangMethod(true);
	skel->Update();
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;

	//Minimum spanning tree to create nodes and backbone node pairs (lines):
	mdl::MST *mst = new mdl::MST( clean_img );
	mst->SetDebug(true);
	mst->SetEdgeRange(10);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	mst->CreateGraphAndMST();
	mst->ErodeAndDialateNodeDegree(50);
	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	//Note: node 0 in bbpairs is index 0 of nodes!!!
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();
	delete mst;
	*/

	
	std::vector<mdl::fPoint3D> nodes;
	std::vector<mdl::pairE> bbpairs;

	mdl::vtkFileHandler file;
	file.SetNodes(&nodes);
	file.SetLines(&bbpairs);
	
	if(!file.Read("BackboneCandidate.vtk"))
	{
		std::cerr << "READ FAILURE\n";
		return EXIT_FAILURE;
	}
	

	mdl::BSplineFitting *bspline = new mdl::BSplineFitting( clean_img );
	bspline->SetDebug(true);
	bspline->SetLevels(8);
	bspline->SetOrder(3);
	bspline->SetNodes( &nodes );
	bspline->SetBBPairs( &bbpairs );
	bspline->Update();
	delete bspline;
   
	getchar();

	//******************************************************************
	//******************************************************************
	
	//******************************************************************
	// IMAGE WRITER
	/*
	typedef itk::ImageFileWriter< mdl::ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("out.tif");
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
	*/
	//******************************************************************
	//*****************************************************************

	return EXIT_SUCCESS;
}
