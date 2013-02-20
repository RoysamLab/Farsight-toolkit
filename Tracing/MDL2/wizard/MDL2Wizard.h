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

#ifndef MDL2WIZARD_H_
#define MDL2WIZARD_H_

#include "vtkSmartPointer.h"
#include <QtGui>
#include "ui_MDL2Wizard.h"

class MDL2WizardHelper;
class QVTKWidget;
class vtkActor;
class vtkImageData;
class vtkRenderer;
class vtkVolume;
class vtkVolumeProperty;

class MDL2Wizard : public QWizard, public Ui::SkeletonizationWizard2
{

  Q_OBJECT

  enum { Page_Intro, Page_Preprocessing, Page_PhaseOne, Page_PhaseTwo };

public:
	MDL2Wizard();
	~MDL2Wizard();
  void Initialize();
	void SetupSignalsAndSlots();
  void RenderVolume(QFileInfo volumeFile);
  void RenderPolyData(QFileInfo polyDataFile);
  void RenderPoints(QString pointsFileName);
  void RenderImage(vtkImageData *image, const char *imageName);
  vtkSmartPointer<vtkActor> CreateActorFromPolyDataFile(const char *filename);
  vtkSmartPointer<vtkActor> CreateActorFromPointsFile(QString pointsFileName);
  vtkSmartPointer<vtkVolume> ConvertRawToVolume(const char *filename);
  vtkSmartPointer<vtkVolumeProperty> NewRGBVolumeProperty(const double range[]);

  //accessors
  void SetInteractiveExecution(bool b) { this->InteractiveExecution = b; }
  bool GetInteractiveExecution() { return this->InteractiveExecution; }

/*
  QString GetConnectedComponentsSize()
    { return this->ConnectedComponentsSize; }
  void SetConnectedComponentsSize(const char *s);
  QString GetVectorMagnitude() { return this->VectorMagnitude; }
  void SetVectorMagnitude(const char *s);
  QString GetEdgeRange() { return this->EdgeRange; }
  void SetEdgeRange(const char *s);
  QString GetMorphStrength() { return this->MorphStrength; }
  void SetMorphStrength(const char *s);
  QString GetBSplineOrder() { return this->BSplineOrder; }
  void SetBSplineOrder(const char *s);
  QString GetBSplineLevels() { return this->BSplineLevels; }
  void SetBSplineLevels(const char *s);
*/

  //properties and their accessors.  Used to control when the next/back/finish
  //buttons are enabled and disabled.  These are public so the WizardPages can
  //access them, but in general they shouldn't be accessed by the user
  QString InputFileName;
  void SetInputFileName(QString s) { this->InputFileName = s; }
  QString GetInputFileName() { return this->InputFileName; }
  QString BackboneFileName;
  void SetBackboneFileName(QString s) { this->BackboneFileName = s; }
  QString GetBackboneFileName() { return this->BackboneFileName; }
  QString SkeletonFileName;
  void SetSkeletonFileName(QString s) { this->SkeletonFileName = s; }
  QString GetSkeletonFileName() { return this->SkeletonFileName; }
  bool PreprocessingDone;
  bool PhaseOneDone;
  bool PhaseTwoDone;

public slots:
  void SelectInputImage();
  void SelectBackboneFile();
  void SelectSkeletonFile();
  void SetInputImage(const char *filename);
  void SetBackboneFile(const char *filename);
  void SetSkeletonFile(const char *filename);
  void MaskUsingGraphCuts();
  void MaskSmallConnComp();
  void Integratedskel();
  void CreateGraphAndMST1();
  void ErodeAndDilateNodeDegree1();
  void BackboneExtract();
  void BSplineFitting();
  void CreateGraphAndMST2();
  void ErodeAndDilateNodeDegree2();
  void SaveOutput();

signals:
  //signals to tell the helper when to run 
  void InputChanged(std::string inputFileName);
  void StartMaskUsingGraphCuts();
  void StartMaskSmallConnComp(int componentSize);
  void StartIntegratedskel(double vectorMagnitude);
  void StartCreateGraphAndMST1(int edgeRange);
  void StartErodeAndDilateNodeDegree1(int morphStrength);
  void StartBackboneExtract();
  void StartBSplineFitting(unsigned int order, unsigned int levels);
  void StartCreateGraphAndMST2(int edgeRange);
  void StartErodeAndDilateNodeDegree2(int morphStrength);
  void ReadyToSaveBackbone(const char *fileName);
  void ReadyToSaveSkeleton(const char *fileName);

protected slots:
  void DisplayMaskUsingGraphCutsResults();
  void DisplayMaskSmallConnCompResults();
  void DisplayIntegratedskelResults();
  void DisplayCreateGraphAndMST1Results();
  void DisplayErodeAndDilateNodeDegree1Results();
  void DisplayBackboneExtractResults();
  void DisplayBSplineFittingResults();
  void DisplayCreateGraphAndMST2Results();
  void DisplayErodeAndDilateNodeDegree2Results();
  void showHelp();
  void UpdateHelpWindow();

protected:
  void SaveParameters();
	void closeEvent(QCloseEvent *event);

private:
  MDL2WizardHelper *Helper;
	QMenu *FileMenu;
	QAction *ExitAction;
  QFileInfo InputFile;
  QFileInfo BackboneFile;
  QFileInfo SkeletonFile;
  QTime Time;
  QTextEdit *OutputWindow;
  QTextEdit *HelpWindow;
  QVTKWidget *RenderWidget;
  vtkSmartPointer<vtkRenderer> Renderer;

  //bool RawInput;
  /*
  QString ConnectedComponentsSize;
  QString VectorMagnitude;
  QString EdgeRange;
  QString MorphStrength;
  QString BSplineOrder;
  QString BSplineLevels;
  QString DataDir;
  */
  bool InteractiveExecution;
};
#endif

