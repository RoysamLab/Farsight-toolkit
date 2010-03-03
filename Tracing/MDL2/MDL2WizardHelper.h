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

#ifndef MDL2WIZARDHELPER_H_
#define MDL2WIZARDHELPER_H_

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"

#include "MDL2Wizard.h"
#include "mdlVolumeProcess.h"
#include "mdlIntegratedSkeleton.h"
#include "mdlMST.h"
#include "mdlBSplineFitting.h"
#include "mdlUtils.h"

#include "vtkImageData.h"

#include <vector>

#include <QThread>

class MDL2WizardHelper : public QThread
{
  Q_OBJECT

public:
	MDL2WizardHelper();
	~MDL2WizardHelper();
  void run();
  vtkImageData* GetImageData() { return this->ITKtoVTK->GetOutput(); }

signals:
  void MaskUsingGraphCutsFinished();
  void MaskSmallConnCompFinished();
  void IntegratedskelFinished();
  void CreateGraphAndMST1Finished();
  void ErodeAndDilateNodeDegree1Finished();
  void BackboneExtractFinished();
  void BSplineFittingFinished();
  void CreateGraphAndMST2Finished();
  void ErodeAndDilateNodeDegree2Finished();

public slots:
  void ReadImage(std::string inputFileName);
  void RunMaskUsingGraphCuts();
  void RunMaskSmallConnComp(int componentSize);
  void RunIntegratedskel(double vectorMagnitude);
  void RunCreateGraphAndMST1(int edgeRange);
  void RunErodeAndDilateNodeDegree1(int morphStrength);
  void RunBackboneExtract();
  void RunBSplineFitting(unsigned int order, unsigned int levels);
  void RunCreateGraphAndMST2(int edgeRange);
  void RunErodeAndDilateNodeDegree2(int morphStrength);
  void WriteBackbone(const char *fileName);
  void WriteSkeleton(const char *fileName);

private:
  typedef itk::ImageFileReader< mdl::ImageType > ReaderType;
  ReaderType::Pointer Reader;
  typedef itk::ImageToVTKImageFilter< mdl::ImageType > ITKtoVTKType;
  ITKtoVTKType::Pointer ITKtoVTK;
  mdl::VolumeProcess *VolProc;
  mdl::IntegratedSkeleton *Skel;
  mdl::MST *MinSpanTree;
  mdl::BSplineFitting *BSpline;
  std::vector<mdl::fPoint3D> SkeletonPoints;
  mdl::ImageType::Pointer CleanImage;
  mdl::ImageType::Pointer EnhancedImage;
  std::vector<mdl::fPoint3D> Nodes;
  std::vector<mdl::pairE> BackbonePairs;
  mdl::vtkFileHandler *FileHandler;
};

#endif

