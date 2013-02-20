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

	if (argc==7)
	{   boolIntensity = atoi(argv[6]);
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
	writer->SetFileName("InputImage.mhd");
    writer->Update();
	
    volProc->RunFillingZeroOnBouandary(4,4,4);
	volProc->RescaleIntensities(0,255);
	volProc->RunBinaryForDistanceMapUsingManualThreshold(threshold);
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
   
	mdl::ImageType::Pointer DT_img = volProc->GetOutput();
	
	if(boolIntensity ==0)
       PreImage = volProc->GetOutput();

	else
	{  
	   reader->SetFileName(inputFileName);
	   mdl::ImageType::Pointer img1 = reader->GetOutput();
	   mdl::VolumeProcess *volProc1 = new mdl::VolumeProcess();
	   volProc1->SetInput(img1);
	   volProc1->RescaleIntensities(0,255);
	   volProc1->RunManualThreshold(1.9);
	   volProc1->RunCAD(100,0.0425,3);
	   
	   writer->SetInput(volProc1->GetOutput());
	   writer->SetFileName("CAD.mhd");
       writer->Update();


	   volProc1->MaskUsingGraphCuts();
	   volProc1->MaskSmallConnComp(10);
	   PreImage = volProc1->GetOutput();

	   writer->SetInput(PreImage);
	   writer->SetFileName("Segmented.mhd");
       writer->Update();

       delete volProc1;
	}
    
	delete volProc;
  

		
	//******************************************************************
	// IMAGE WRITER
    writer->SetInput(DT_img);
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

   
    mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( PreImage );
	skel->SetVectorMagnitude(VectorMagnitude);
	skel->SetDebug(true);
	skel->SetUseXiaoLiangMethod(false);
	skel->Update();
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;

	std::cerr << "MST 1\n";

	//Minimum spanning tree to create nodes and backbone node pairs (lines):
	mdl::MST *mst = new mdl::MST( DT_img );
	mst->SetDebug(false);
	mst->SetUseVoxelRounding(true);
	mst->SetEdgeRange(edgeRange);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	// can choose different weight
	//mst->CreateGraphAndMST(3);
	mst->CreateGraphAndMST(1);
	mst->ErodeAndDialateNodeDegree(morphStrength);
	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	//Note: node 0 in bbpairs is index 0 of nodes!!!
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();
	delete mst;
   

	std::string outputFileName0 = inputFileName;
	size_t found0 = outputFileName0.find_last_of(".");
	outputFileName0.insert(found0,"Candidate_Backbone");
	found0 = outputFileName0.find_last_of(".");
	outputFileName0.replace(found0+1,3,"vtk");
	
	mdl::vtkFileHandler * fhd0 = new mdl::vtkFileHandler();
	fhd0->SetNodes(&nodes);
	fhd0->SetLines(&bbpairs);
	fhd0->Write(outputFileName0.c_str());
	delete fhd0;


	
	std::cerr << "BSPLINE\n";
	
	mdl::BSplineFitting *bspline = new mdl::BSplineFitting( DT_img );
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
  
    mdl::MST *mst1 = new mdl::MST( DT_img );
	mst1->SetDebug(false);
	mst1->SetUseVoxelRounding(false);
	mst1->SetEdgeRange(edgeRange);
	mst1->SetPower(1);
	mst1->SetSkeletonPoints( &nodes );
	mst1->CreateGraphAndMST(1);
	mst1->ErodeAndDialateNodeDegree(0);
	
	nodes = mst1->GetNodes();
	bbpairs = mst1->BackboneExtract();
	delete mst1;


	std::cerr << "Saving\n";
	
	std::string outputFileName = inputFileName;
	size_t found = outputFileName.find_last_of(".");
	outputFileName.insert(found,"_Backbone");
	found = outputFileName.find_last_of(".");
	outputFileName.replace(found+1,3,"vtk");
	


	mdl::vtkFileHandler * fhdl = new mdl::vtkFileHandler();
	fhdl->SetNodes(&nodes);
	fhdl->SetLines(&bbpairs);
	fhdl->Write(outputFileName.c_str());
	delete fhdl;

	std::string outputFileName1 = inputFileName;
	size_t found1 = outputFileName1.find_last_of(".");
	outputFileName1.insert(found1,"_Skeleton");
	found1 = outputFileName1.find_last_of(".");
	outputFileName1.replace(found1+1,3,"vtk");

	mdl::vtkFileHandler * fhdl2 = new mdl::vtkFileHandler();
	fhdl2->SetNodes(&skeleton);
	fhdl2->SetLines(&bbpairs);
	fhdl2->Write(outputFileName1.c_str());
	delete fhdl2;

	std::cerr << "DONE\n";
  
	//******************************************************************
	//******************************************************************


	return EXIT_SUCCESS;
}
