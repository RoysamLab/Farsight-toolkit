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

#ifndef MDLWIZARD_H_
#define MDLWIZARD_H_

#include "vtkSmartPointer.h"
#include <QtGui>
#include "ui_MDLWizard.h"

class QVTKWidget;
class vtkActor;
class vtkRenderer;
class vtkVolume;
class vtkVolumeProperty;

class MDLWizard : public QWizard, public Ui::SkeletonizationWizard
{

  Q_OBJECT

  enum { Page_Intro, Page_Preprocessing, Page_PhaseOne, Page_PhaseTwo };

public:
	MDLWizard();
	~MDLWizard();
  void Initialize();
	void SetupSignalsAndSlots();
  void RegisterFields();
  int Test(int argc, char **argv);
  void RenderVolume(QFileInfo volumeFile);
  void RenderPolyData(QFileInfo polyDataFile);
  void RenderPoints(QString pointsFileName);
  vtkSmartPointer<vtkActor> CreateActorFromPolyDataFile(const char *filename);
  vtkSmartPointer<vtkActor> CreateActorFromPointsFile(QString pointsFileName);
  vtkSmartPointer<vtkVolume> ConvertRawToVolume(const char *filename);
  vtkSmartPointer<vtkVolumeProperty> NewRGBVolumeProperty(const double range[]);

  //accessors
  void SetInteractiveExecution(bool b) { this->InteractiveExecution = b; }
  bool GetInteractiveExecution() { return this->InteractiveExecution; }
  QString GetConnectedComponentsSize()
    { return this->ConnectedComponentsSize; }
  void SetConnectedComponentsSize(const char *s);
  QString GetVectorMagnitude() { return this->VectorMagnitude; }
  void SetVectorMagnitude(const char *s);
  QString GetEdgeRange() { return this->EdgeRange; }
  void SetEdgeRange(const char *s);
  QString GetGraphPruneSize() { return this->GraphPruneSize; }
  void SetGraphPruneSize(const char *s);
  QString GetMorphStrength() { return this->MorphStrength; }
  void SetMorphStrength(const char *s);
  QString GetWeightFactor() { return this->WeightFactor; }
  void SetWeightFactor(const char *s);
  QString GetBSplineOrder() { return this->BSplineOrder; }
  void SetBSplineOrder(const char *s);
  QString GetBSplineLevels() { return this->BSplineLevels; }
  void SetBSplineLevels(const char *s);

  //properties and their accessors.  Used to control when the next/back/finish
  //buttons are enabled and disabled.  These are public so the WizardPages can
  //access them, but in general they shouldn't be accessed by the user
  QString InputFileName;
  void SetInputFileName(QString s) { this->InputFileName = s; }
  QString GetInputFileName() { return this->InputFileName; }
  QString BackboneFileName;
  void SetBackboneFileName(QString s) { this->BackboneFileName = s; }
  QString GetBackboneFileName() { return this->BackboneFileName; }
  QString SpinesFileName;
  void SetSpinesFileName(QString s) { this->SpinesFileName = s; }
  QString GetSpinesFileName() { return this->SpinesFileName; }
  bool PreprocessingDone;
  bool PhaseOneDone;
  bool PhaseTwoDone;

public slots:
  void SelectInputImage();
  void SelectBackboneFile();
  void SelectSpinesFile();
  void SetInputImage(const char *filename);
  void SetBackboneFile(const char *filename);
  void SetSpinesFile(const char *filename);
  void AppendOutputToDisplay(QObject *mapped);
  void RunVolumeProcess();
  void RunConnCompntwFldfill();
  void RunAnisoDiffuse();
  void RunGradientVecField();
  void RunIntegratedskel();
  void RunBackboneExtract1();
  void RunMDABasedSpineExtraction1();
  void RunBSplineFitting();
  void RunRefiningSkeleton1();
  void RunBackboneExtract2();
  void RunRefiningSkeleton2();
  void RunMDABasedSpineExtraction2();

signals:
  void inputChanged();
  void backboneFileChanged();
  void spinesFileChanged();
  void PreprocessingFinished();
  void PhaseOneFinished();
  void PhaseTwoFinished();

protected slots:
  void UpdateEdgeRangeFrom1();
  void UpdateEdgeRangeFrom2();
  void UpdateMorphStrengthFrom1();
  void UpdateMorphStrengthFrom2();
  void UpdateWeightFactorFrom1();
  void UpdateWeightFactorFrom2();
  void VolumeProcessFinished();
  void ConnCompntwFldfillFinished();
  void AnisoDiffuseFinished();
  void GradientVecFieldFinished();
  void IntegratedskelFinished();
  void BackboneExtract1Finished();
  void MDABasedSpineExtraction1Finished();
  void BSplineFittingFinished();
  void RefiningSkeleton1Finished();
  void BackboneExtract2Finished();
  void RefiningSkeleton2Finished();
  void MDABasedSpineExtraction2Finished();
  void DeleteIntermediaryFiles();
  void showHelp();
  void UpdateHelpWindow();

protected:
	void closeEvent(QCloseEvent *event);

private:
  QSignalMapper *Mapper;
	QMenu *FileMenu;
	QAction *ExitAction;
  QFileInfo InputFile;
  QFileInfo BackboneFile;
  QFileInfo SpinesFile;
  QTime Time;
  QProcess *VolumeProcess;
  QProcess *ConnCompntwFldfill;
  QProcess *AnisoDiffuse;
  QProcess *GradientVecField;
  QProcess *Integratedskel;
  QProcess *BackboneExtract1;
  QProcess *MDABasedSpineExtraction1;
  QProcess *BSplineFitting;
  QProcess *RefiningSkeleton1;
  QProcess *BackboneExtract2;
  QProcess *RefiningSkeleton2;
  QProcess *MDABasedSpineExtraction2;
  QTextEdit *OutputWindow;
  QTextEdit *HelpWindow;
  QVTKWidget *RenderWidget;
  vtkSmartPointer<vtkRenderer> Renderer;
  bool RawInput;
  QString ImageSizeX;
  QString ImageSizeY;
  QString ImageSizeZ;
  QString ConnectedComponentsSize;
  QString VectorMagnitude;
  QString EdgeRange;
  QString GraphPruneSize;
  QString MorphStrength;
  QString WeightFactor;
  QString BSplineOrder;
  QString BSplineLevels;
  QString VolumeProcessedFile;
  QString ComponentsConnectedFile;
  QString AnisoDiffusedFile;
  QString VectorFile;
  QString SeedFile;
  QString DataDir;
  QString SkeletonFile;
  QString BackboneCandidateFile;
  QString MDLFeatureFile;
  QString SpineCandidateFile;
  QString SmoothBackboneFile;
  QString ExtraSpineFile;
  QString RefinedSeedFile;
  QString RefinedSkeletonFile;
  QString ExecutablePath;
  bool InteractiveExecution;
};
#endif

