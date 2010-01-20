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

public:
	MDLWizard();
	~MDLWizard();
  void Initialize();
	void SetupSignalsAndSlots();
  void RegisterFields();
  int Test(int argc, char **argv);
  void RenderVolume(QFileInfo volumeFile);
  void RenderPolyData(QFileInfo polyDataFile);
  vtkSmartPointer<vtkVolume> ConvertRawToVolume(const char *filename);
  vtkSmartPointer<vtkActor> CreateActorFromPolyDataFile(const char *filename);
  vtkSmartPointer<vtkVolumeProperty> NewRGBVolumeProperty();

  //properties and their accessors.  Used to control when the next/back/finish
  //buttons are enabled and disabled.
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
  void MDABaseSpineExtraction1Finished();
  void BSplineFittingFinished();
  void RefiningSkeleton1Finished();
  void BackboneExtract2Finished();
  void RefiningSkeleton2Finished();
  void MDABasedSpineExtraction2Finished();
  void DeleteIntermediaryFiles();

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
  QTextEdit *OutputDisplay;
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
  QString RefinedSkeletonFile;
  QString ExecutablePath;
};
#endif

