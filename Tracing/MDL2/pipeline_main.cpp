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
	volProc->RescaleIntensities(0,255);
	//volProc->RunManualThreshold(1.9);
	volProc->RunCAD();
	volProc->MaskUsingGraphCuts();
	volProc->MaskSmallConnComp(10);
	mdl::ImageType::Pointer clean_img = volProc->GetOutput();
	//volProc->RunDistanceTransform();
	//mdl::ImageType::Pointer DT_img = volProc->GetOutput();
	volProc->RunAnisotropicDiffusion(1,false);
	mdl::ImageType::Pointer enhance_img = volProc->GetOutput();
	delete volProc;

	
	//Integrated Skeleton to create skeleton points:
	mdl::IntegratedSkeleton *skel = new mdl::IntegratedSkeleton( clean_img );
	skel->SetVectorMagnitude(.05);
	skel->SetDebug(true);
	skel->SetUseXiaoLiangMethod(false);
	skel->Update();
	std::vector<mdl::fPoint3D> skeleton = skel->GetOutput();
	delete skel;

	//Minimum spanning tree to create nodes and backbone node pairs (lines):
	mdl::MST *mst = new mdl::MST( clean_img );
	mst->SetDebug(true);
	mst->SetEdgeRange(12);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &skeleton );
	mst->CreateGraphAndMST();
	mst->ErodeAndDialateNodeDegree(50);
	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	//Note: node 0 in bbpairs is index 0 of nodes!!!
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();
	//mst->SpineExtract();
	delete mst;
	

	/*
	std::vector<mdl::fPoint3D> nodes;
	std::vector<mdl::pairE> bbpairs;

	mdl::vtkFileHandler file;
	file.SetNodes(&nodes);
	file.SetLines(&bbpairs);
	
	if(!file.Read("BackboneCandidateXiao.vtk"))
	{
		std::cerr << "READ FAILURE\n";
		return EXIT_FAILURE;
	}
    */
	
	mdl::BSplineFitting *bspline = new mdl::BSplineFitting( clean_img );
	bspline->SetDebug(true);
	bspline->SetLevels(7);
	bspline->SetOrder(3);
	bspline->SetNodes( &nodes );
	bspline->SetBBPairs( &bbpairs );
	bspline->Update();
	nodes = bspline->GetNodes();
	bbpairs = bspline->GetBBPairs();
	delete bspline;
    
	/*mdl::vtkFileHandler file;
	file.SetNodes(&nodes);
	file.SetLines(&bbpairs);
	
	if(!file.Read("BBSpline.vtk"))
	{
		std::cerr << "READ FAILURE\n";
		return EXIT_FAILURE;
	}
     */


  
    mdl::MST *mst1 = new mdl::MST( clean_img );
	//mst1->SetInputforSpineExtraction(clean_img,(char *)"RealSpinePrior.txt",(char *)"NonSpinePrior.txt");
	mst1->SetDebug(true);
	mst1->SetEdgeRange(12);
	mst1->SetPower(1);
	mst1->SetVesselMap(enhance_img);
	mst1->SetAlpha(0.7);
	mst1->SetPruneThreshold(4.0);
	//mst1->SetFileofRealSpineFeature((char *)"RealSpinePrior.txt");
	//mst1->SetFileofNonSpineFeature((char *)"NonSpinePrior.txt");
	mst1->SetSkeletonPoints( &skeleton );
	mst1->CreateGraphAndMST();
	mst1->ErodeAndDialateNodeDegree(50);
	
	nodes = mst1->GetNodes();
	//bbpairs = mst1->SearchFirstandSecondLevelBranch();
	mst1->SearchFirstandSecondLevelBranch();
	//bbpairs = mst->BackboneExtract();
	bbpairs = mst1->SpineExtract();
	delete mst1;	
  
	std::cerr << "PRESS ENTER TO EXIT\n";
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
