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

#ifndef MDLGUI_H_
#define MDLGUI_H_

#include <QtGui>
#include "ui_MDLGUI.h"

class MDLGUI : public QMainWindow, public Ui::MDLForm
{
Q_OBJECT;
public:
	MDLGUI();
	~MDLGUI();
  void Initialize();
	void SetupSignalsAndSlots();

public slots:
  void SelectInputImage();
  void SelectBackboneImage();
  void SelectSpinesFile();
  void CheckInputs();
  void AppendOutputToDisplay(QObject *o);
  void RunSkeletonization();
  void RunvolumeProcess();
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
  void FinishedRunningSkeletonization();

protected:
	void closeEvent(QCloseEvent *event);
  void DeleteIntermediaryFiles();

private:
  QSignalMapper *Mapper;
	QMenu *FileMenu;
	QAction *ExitAction;
  QFileInfo InputFile;
  QFileInfo BackboneFile;
  QFileInfo SpinesFile;
  QTime Time;
  QProcess *volumeProcess;
  QProcess *ConnCompntwFldfill;
  QProcess *AnisoDiffuse;
  QProcess *GradientVecField;
  QProcess *Integratedskel;
  QProcess *BackboneExtract1;
  QProcess *MDABasedSpineExtraction1;
  QProcess *BSplineFitting;;
  QProcess *RefiningSkeleton1;
  QProcess *BackboneExtract2;
  QProcess *RefiningSkeleton2;
  QProcess *MDABasedSpineExtraction2;
  bool RawInput;
  QString ImageSizeX;
  QString ImageSizeY;
  QString ImageSizeZ;
  QString ConnectedComponentsSize;
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
};
#endif

