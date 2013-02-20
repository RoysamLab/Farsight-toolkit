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

#include "MDL2WizardHelper.h"

//-----------------------------------------------------------------------------
MDL2WizardHelper::MDL2WizardHelper()
{
  this->VolProc = 0; 
  this->Skel = 0;
  this->MinSpanTree = 0;
  this->BSpline = 0;
  this->FileHandler = 0;
  this->Reader = ReaderType::New();
  this->ITKtoVTK = ITKtoVTKType::New();
}

//-----------------------------------------------------------------------------
MDL2WizardHelper::~MDL2WizardHelper()
{
  if(this->VolProc)
    {
    delete this->VolProc;
    this->VolProc = 0;
    }
  if(this->Skel)
    {
    delete this->Skel;
    this->Skel = 0;
    }
  if(this->MinSpanTree)
    {
    delete this->MinSpanTree;
    this->MinSpanTree = 0;
    }
  if(this->BSpline)
    {
    delete this->BSpline;
    this->BSpline = 0;
    }
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::run()
{
  this->exec();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::ReadImage(std::string inputFileName)
{
  this->Reader->SetFileName(inputFileName);
  try
    {
    this->Reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "READER FAILED: " << err << std::endl ;
    //emit a signal about failing here
    }
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunMaskUsingGraphCuts()
{
  if(this->VolProc == 0)
    {
    this->VolProc = new mdl::VolumeProcess();
    }
  this->VolProc->SetInput(this->Reader->GetOutput());
  this->VolProc->SetDebug(true);
  this->VolProc->MaskUsingGraphCuts();
  emit this->MaskUsingGraphCutsFinished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunMaskSmallConnComp(int componentsSize)
{
  this->VolProc->MaskSmallConnComp(componentsSize);
  this->CleanImage = this->VolProc->GetOutput();
  this->ITKtoVTK->SetInput(this->CleanImage);
  this->ITKtoVTK->Update();
  emit this->MaskSmallConnCompFinished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunIntegratedskel(double vectorMagnitude)
{
  if(this->Skel != 0)
    {
    delete this->Skel;
    }
  this->Skel = new mdl::IntegratedSkeleton( this->CleanImage );
  this->Skel->SetDebug(true);
  this->Skel->SetUseXiaoLiangMethod(false);
  this->Skel->SetVectorMagnitude(vectorMagnitude);
  this->Skel->Update();
  this->SkeletonPoints = this->Skel->GetOutput(); 
  emit this->IntegratedskelFinished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunCreateGraphAndMST1(int edgeRange)
{
  if(this->MinSpanTree != 0)
    {
    delete this->MinSpanTree;
    }
  this->MinSpanTree = new mdl::MST( this->CleanImage );
  this->MinSpanTree->SetDebug(true);
  this->MinSpanTree->SetUseVoxelRounding(true);
  this->MinSpanTree->SetEdgeRange(edgeRange);
  this->MinSpanTree->SetPower(1);
  this->MinSpanTree->SetSkeletonPoints( &(this->SkeletonPoints) );
  this->MinSpanTree->CreateGraphAndMST();
  emit this->CreateGraphAndMST1Finished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunErodeAndDilateNodeDegree1(int morphStrength)
{
  this->MinSpanTree->ErodeAndDialateNodeDegree(morphStrength);
  this->Nodes = this->MinSpanTree->GetNodes();
  emit this->ErodeAndDilateNodeDegree1Finished(); 

}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunBackboneExtract()
{
  this->BackbonePairs = this->MinSpanTree->BackboneExtract();
  emit this->BackboneExtractFinished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunBSplineFitting(unsigned int order, unsigned int levels)
{
  if(this->BSpline != 0)
    {
    delete this->BSpline;
    }
  this->BSpline = new mdl::BSplineFitting( this->CleanImage );
  this->BSpline->SetDebug(true);
  this->BSpline->SetLevels(levels);
  this->BSpline->SetOrder(order);
  this->BSpline->SetNodes( &(this->Nodes) );
  this->BSpline->SetBBPairs( &(this->BackbonePairs) );
  this->BSpline->Update();
  this->Nodes = this->BSpline->GetNodes();
  this->BackbonePairs = this->BSpline->GetBBPairs();
  emit this->BSplineFittingFinished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunCreateGraphAndMST2(int edgeRange)
{
  if(this->MinSpanTree != 0)
    {
    delete this->MinSpanTree;
    }
  this->MinSpanTree = new mdl::MST( this->CleanImage );
  this->MinSpanTree->SetDebug(true);
  this->MinSpanTree->SetUseVoxelRounding(false);
  this->MinSpanTree->SetEdgeRange(edgeRange);
  this->MinSpanTree->SetPower(1);
  //this->MinSpanTree->SetVesselMap(this->EnhancedImage);
  this->MinSpanTree->SetSkeletonPoints( &(this->SkeletonPoints) );
  this->MinSpanTree->CreateGraphAndMST();
  emit this->CreateGraphAndMST2Finished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::RunErodeAndDilateNodeDegree2(int morphStrength)
{
  this->MinSpanTree->ErodeAndDialateNodeDegree(morphStrength);
  this->Nodes = this->MinSpanTree->GetNodes();
  //this should probably be a separate function...
  this->BackbonePairs = this->MinSpanTree->BackboneExtract();
  emit this->ErodeAndDilateNodeDegree2Finished();
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::WriteBackbone(const char *fileName)
{
  if(this->FileHandler != 0)
    {
    delete this->FileHandler;
    }
  this->FileHandler = new mdl::vtkFileHandler();
	this->FileHandler->SetNodes(& (this->Nodes) );
	this->FileHandler->SetLines(& (this->BackbonePairs) );
	this->FileHandler->Write(fileName);
}

//-----------------------------------------------------------------------------
void MDL2WizardHelper::WriteSkeleton(const char *fileName)
{
  if(this->FileHandler != 0)
    {
    delete this->FileHandler;
    }
  this->FileHandler = new mdl::vtkFileHandler();
	this->FileHandler->SetNodes( &(this->SkeletonPoints) );
	this->FileHandler->Write(fileName);
}



