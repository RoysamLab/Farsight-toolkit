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
		std::cerr << "Usage: " << argv[0] << " inputfilename < useGraphCuts minConnCompSize vectorMagnitude>" << std::endl;
		return EXIT_FAILURE;
    }

	std::string inputFileName = argv[1];
    
	/*
	if(argc == 5)
	{
		useGraphCuts = atoi(argv[2]);
	}
    */

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
    
	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("InputImage.mhd");
    writer->Update();

    volProc->RunFillingZeroOnBouandary(4,4,4);
	volProc->SetDebug(true);
  
	volProc->RescaleIntensities(0,255);

	volProc->RunBinaryForDistanceMapUsingGraphCuts();
	//volProc->RunBinaryForDistanceMapUsingManualThreshold(20);

	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("BinaryResults.mhd");
    writer->Update();

	volProc->RunBinayMedianHoleFilling(1,1,1);
    writer->SetInput(volProc->GetOutput());
	writer->SetFileName("RunBinayMedianHoleFillingResults.mhd");
    writer->Update();
    

	volProc->RunDanielssonDistanceMap();
	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("DanielssonDistanceMapResults.mhd");
    writer->Update();
    
	volProc->RunGaussianSmoothing(7,7,3,0.1);

	writer->SetInput(volProc->GetOutput());
	writer->SetFileName("RunGaussianSmoothingResults.mhd");
    writer->Update();
	//volProc->RunRecursiveGaussianIIRRilter(7,5,2);
	//volProc->RunCAD();

	mdl::ImageType::Pointer clean_img = volProc->GetOutput();
	delete volProc;
  

		
	//******************************************************************
	// IMAGE WRITER
    writer->SetInput(clean_img);
	writer->SetFileName("MaskUsingDistanceMap.mhd");
    writer->Update();

	
	//******************************************************************
	std::cerr << "Integrated Skeleton\n";
     
    //********************** New Method **********************************
    /*   
	mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( clean_img );
	skel->SetDebug(true);	
	skel->RunXiaoLSkeletonPoints(0,4);
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;
     */

	mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( clean_img );
	skel->SetVectorMagnitude(0.002);
	skel->SetDebug(true);
	skel->SetUseXiaoLiangMethod(false);
	skel->Update();
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;

	std::cerr << "MST 1\n";

	//Minimum spanning tree to create nodes and backbone node pairs (lines):
	mdl::MST *mst = new mdl::MST( clean_img );
	mst->SetDebug(false);
	mst->SetUseVoxelRounding(true);
	mst->SetEdgeRange(8);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	// can choose different weight
	//mst->CreateGraphAndMST(3);
	mst->CreateGraphAndMST(1);
	mst->ErodeAndDialateNodeDegree(50);
	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	//Note: node 0 in bbpairs is index 0 of nodes!!!
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();
	delete mst;
   

	
	std::cerr << "BSPLINE\n";
	
	mdl::BSplineFitting *bspline = new mdl::BSplineFitting( clean_img );
	bspline->SetDebug(false);
	bspline->SetLevels(6);
	bspline->SetOrder(3);
	bspline->SetNodes( &nodes );
	bspline->SetBBPairs( &bbpairs );
	bspline->Update();
	nodes = bspline->GetNodes();
	bbpairs = bspline->GetBBPairs();
	delete bspline;

 
	std::cerr << "MST 2\n";
  
    mdl::MST *mst1 = new mdl::MST( clean_img );
	mst1->SetDebug(false);
	mst1->SetUseVoxelRounding(false);
	mst1->SetEdgeRange(8);
	mst1->SetPower(1);
	mst1->SetSkeletonPoints( &nodes );
	mst1->CreateGraphAndMST(1);
	mst1->ErodeAndDialateNodeDegree(0);
	
	nodes = mst1->GetNodes();
	bbpairs = mst1->BackboneExtract();
	delete mst1;


	std::cerr << "Saving\n";

	mdl::vtkFileHandler * fhdl = new mdl::vtkFileHandler();
	fhdl->SetNodes(&nodes);
	fhdl->SetLines(&bbpairs);
	fhdl->Write("Backbone.vtk");
	delete fhdl;

	mdl::vtkFileHandler * fhdl1 = new mdl::vtkFileHandler();
	fhdl1->SetNodes(&skeleton);
	fhdl1->Write("Skeleton.vtk");
	delete fhdl1;

	std::cerr << "DONE\n";
  
	std::cerr << "PRESS ENTER TO EXIT\n";
	getchar();

	//******************************************************************
	//******************************************************************


	return EXIT_SUCCESS;
}
