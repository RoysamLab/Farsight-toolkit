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

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <QAction>
#include "MDLGUI.h"

MDLGUI::MDLGUI()
{
  this->volumeProcess = new QProcess(this);
  this->ConnCompntwFldfill = new QProcess(this);
  this->AnisoDiffuse = new QProcess(this);
  this->GradientVecField = new QProcess(this);
  this->Integratedskel = new QProcess(this);
  this->BackboneExtract1 = new QProcess(this);
  this->MDABasedSpineExtraction1 = new QProcess(this);
  this->BSplineFitting = new QProcess(this);
  this->RefiningSkeleton1 = new QProcess(this);
  this->BackboneExtract2 = new QProcess(this);
  this->RefiningSkeleton2 = new QProcess(this);
  this->MDABasedSpineExtraction2 = new QProcess(this);
  this->Mapper = new QSignalMapper(this);
  this->setObjectName("MDLGUI");
  this->Initialize();
}

MDLGUI::~MDLGUI()
{
}

void MDLGUI::Initialize()
{
  this->setupUi(this);
  this->move(40, 59);
  this->setWindowTitle(tr("MDL"));
  this->SetupSignalsAndSlots();
}

void MDLGUI::SetupSignalsAndSlots()
{
  this->connect(this->SelectInputButton, SIGNAL(pressed()), this,
    SLOT(SelectInputImage()));
  this->connect(this->SelectBackboneButton, SIGNAL(pressed()), this,
    SLOT(SelectBackboneImage()));
  this->connect(this->SelectSpinesButton, SIGNAL(pressed()), this,
    SLOT(SelectSpinesFile()));
  this->connect(this->ComponentsSizeInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->VectorMagnitudeInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->EdgeRangeInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->GraphPruneSizeInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->MorphStrengthInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->WeightFactorInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->RunSkeletonizationButton, SIGNAL(pressed()), this,
    SLOT(RunSkeletonization()));
  
  //set up the QSignalMapper so that each process can append to the
  //DisplayOutput field
  this->Mapper->setMapping(this->AnisoDiffuse, this->AnisoDiffuse);
  this->Mapper->setMapping(this->BSplineFitting, this->BSplineFitting);
  this->Mapper->setMapping(this->ConnCompntwFldfill, this->ConnCompntwFldfill);
  this->Mapper->setMapping(this->GradientVecField, this->GradientVecField);
  this->Mapper->setMapping(this->Integratedskel, this->Integratedskel);
  this->Mapper->setMapping(this->volumeProcess, this->volumeProcess);
  this->Mapper->setMapping(this->BackboneExtract1, this->BackboneExtract1);
  this->Mapper->setMapping(this->MDABasedSpineExtraction1,
                           this->MDABasedSpineExtraction1);
  this->Mapper->setMapping(this->RefiningSkeleton1, this->RefiningSkeleton1);
  this->Mapper->setMapping(this->BackboneExtract2, this->BackboneExtract2);
  this->Mapper->setMapping(this->RefiningSkeleton2, this->RefiningSkeleton2);
  this->Mapper->setMapping(this->MDABasedSpineExtraction2,
                           this->MDABasedSpineExtraction2);
  connect(this->AnisoDiffuse, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
  connect(this->BSplineFitting, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
  connect(this->ConnCompntwFldfill, SIGNAL(readyReadStandardOutput()),
    this->Mapper, SLOT(map()));
  connect(this->GradientVecField, SIGNAL(readyReadStandardOutput()),
    this->Mapper, SLOT(map()));
  connect(this->Integratedskel, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
  connect(this->volumeProcess, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
  connect(this->BackboneExtract1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->MDABasedSpineExtraction1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->RefiningSkeleton1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->BackboneExtract2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->RefiningSkeleton2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->MDABasedSpineExtraction2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));

  connect(this->Mapper, SIGNAL(mapped(QObject *)), this,
          SLOT(AppendOutputToDisplay(QObject *)));

  //connect the process so that they run sequentially
  this->connect(this->volumeProcess, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunConnCompntwFldfill()));
  this->connect(this->ConnCompntwFldfill,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunAnisoDiffuse()));
  this->connect(this->AnisoDiffuse, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunGradientVecField()));
  this->connect(this->GradientVecField, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunIntegratedskel()));
  this->connect(this->Integratedskel, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunBackboneExtract1()));
  this->connect(this->BackboneExtract1, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunMDABasedSpineExtraction1()));
  this->connect(this->MDABasedSpineExtraction1, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunBSplineFitting()));
  this->connect(this->BSplineFitting, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunRefiningSkeleton1()));
  this->connect(this->RefiningSkeleton1, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunBackboneExtract2()));
  this->connect(this->BackboneExtract2, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunRefiningSkeleton2()));
  this->connect(this->RefiningSkeleton2, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunMDABasedSpineExtraction2()));
  this->connect(this->MDABasedSpineExtraction2, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(FinishedRunningSkeletonization()));
}

void MDLGUI::SelectInputImage()
{
  QString filePath = QFileDialog::getOpenFileName(this,
    tr("Select Input Image"), "", tr("Image Files (*.*)"));
  this->InputFile = QFileInfo(filePath);
  this->InputImageLabel->setText(this->InputFile.fileName());
  this->CheckInputs();
}

void MDLGUI::SelectBackboneImage()
{
  QString filePath = QFileDialog::getSaveFileName(this,
    tr("Select Output File"), "", tr("VTK files (*.vtk)"));
  this->BackboneFile = QFileInfo(filePath);
  this->DataDir = this->BackboneFile.dir().absolutePath() + "/";
  this->OutputImageLabel->setText(this->BackboneFile.fileName());
  this->CheckInputs();
}

void MDLGUI::SelectSpinesFile()
{
  QString filePath = QFileDialog::getSaveFileName(this,
    tr("Select Spines File"), "", tr("VTK files (*.vtk)"));
  this->SpinesFile = QFileInfo(filePath);
  this->SpinesLabel->setText(this->SpinesFile.fileName());
  this->CheckInputs();
}

void MDLGUI::CheckInputs()
{
  if(this->InputFile.filePath() != "" && this->BackboneFile.filePath() != "" &&
     this->ComponentsSizeInput->text() != "" &&
     this->VectorMagnitudeInput->text() != "" &&
     this->EdgeRangeInput->text() != "" &&
     this->GraphPruneSizeInput->text() != "" &&
     this->MorphStrengthInput->text() != "" &&
     this->WeightFactorInput->text() != "")
    {
    this->RunSkeletonizationButton->setEnabled(true);
    }
  else
    {
    this->RunSkeletonizationButton->setEnabled(false);
    }
}

void MDLGUI::AppendOutputToDisplay(QObject *o)
{
  QProcess *process = static_cast<QProcess *>(o);
  this->OutputDisplay->append(QString(process->readAllStandardOutput())); 
}

void MDLGUI::RunSkeletonization()
{
  //disable the "run" button while a skeletonization is in progress
  this->RunSkeletonizationButton->setEnabled(false);
  this->OutputDisplay->clear();

  //initialize some variables before we get started
  if(this->InputFile.suffix() == "raw")
    {
    this->RawInput = true;
    QRegExp rx("(\\d+)x(\\d+)x(\\d+)");
    if(this->InputFile.fileName().indexOf(rx) != -1)
      {
      QStringList list = rx.capturedTexts();
      this->ImageSizeX = list[1];
      this->ImageSizeY = list[2];
      this->ImageSizeZ = list[3];
      }
    else
      {
      QMessageBox msgBox;
      msgBox.setText(".raw files must have the image dimensions specified in the filename in the following format: img.640x480x24.raw");
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
      this->InputFile = QFileInfo();
      this->InputImageLabel->setText("");
      return;
      }
    }
  else
    {
    this->RawInput = false;
    }
  this->ConnectedComponentsSize = this->ComponentsSizeInput->text();
  this->RunvolumeProcess();
}

void MDLGUI::RunvolumeProcess()
{
  //construct the list of arguments
  QStringList arguments;
  this->VolumeProcessedFile =
    this->BackboneFile.dir().absolutePath() + "/volume_Processed.raw";
  if(this->RawInput)
    {
    arguments << this->InputFile.filePath() << this->ImageSizeX
              << this->ImageSizeY << this->ImageSizeZ
              << this->VolumeProcessedFile;
    }
  else
    {
    arguments << this->InputFile.filePath() << this->VolumeProcessedFile;
    }

  //provide some feedback to the user
  this->OutputDisplay->append("Executing volumeProcess\n");

  //start a timer
  this->Time.start();
  //run the executable
  this->volumeProcess->start("./volumeProcess", arguments);
}

void MDLGUI::RunConnCompntwFldfill()
{
  //provide information about how long the last step took.
  int secondsElapsed = this->Time.elapsed() / 1000.0;

  this->OutputDisplay->append(
    "volumeProcess completed in " + QString::number(secondsElapsed) + " seconds.\n");

  if(!(this->RawInput))
    {
    //also have to grab image size out of the output
    QRegExp rx("input image size: \\[(\\d+), (\\d+), (\\d+)\\]"); 
    if(this->OutputDisplay->document()->toPlainText().indexOf(rx) != -1)
      {
      QStringList list = rx.capturedTexts();
      this->ImageSizeX = list[1];
      this->ImageSizeY = list[2];
      this->ImageSizeZ = list[3];
      }
    else
      {
      QMessageBox msgBox;
      msgBox.setText("BUG: Didn't find the image size in output...");
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
      }
    }

  //construct a list of arguments
  this->ComponentsConnectedFile = 
    this->BackboneFile.dir().absolutePath() + "/components_Connected.raw";
  QStringList arguments;
  arguments << this->VolumeProcessedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ
            << this->ComponentsConnectedFile << this->ConnectedComponentsSize;
  
  this->OutputDisplay->append("Executing ConnCompntwFldfill\n");
  this->Time.restart();
  this->ConnCompntwFldfill->start("./ConnCompntwFldfill", arguments);
}

void MDLGUI::RunAnisoDiffuse()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "ConnCompntwFldfill completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->AnisoDiffusedFile =
    this->BackboneFile.dir().absolutePath() + "/Aniso_Diffused.raw";

  QStringList arguments;
  arguments << this->ComponentsConnectedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ << this->AnisoDiffusedFile;

  this->OutputDisplay->append("Executing AnisoDiffuse\n");
  this->Time.restart();
  this->AnisoDiffuse->start("./AnisoDiffuse", arguments);
}

void MDLGUI::RunGradientVecField()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "AnisoDiffuse completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->VectorFile = this->BackboneFile.dir().absolutePath() + "/out.vec";

  QStringList arguments;
  arguments << this->AnisoDiffusedFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << this->VectorFile;

  this->OutputDisplay->append("Executing GradientVecField\n");
  this->Time.restart();
  this->GradientVecField->start("./GradientVecField", arguments);
}

void MDLGUI::RunIntegratedskel()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "GradientVecField completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->SeedFile = this->BackboneFile.dir().absolutePath() + "/out.seed";
  this->SkeletonFile = this->BackboneFile.dir().absolutePath() + "/out.skel";
  this->VectorMagnitude = this->VectorMagnitudeInput->text();

  QStringList arguments;
  arguments << this->VectorFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << this->VectorMagnitude << this->SeedFile
            << this->SkeletonFile;

  this->OutputDisplay->append("Executing Integratedskel\n");
  this->Time.restart();
  this->Integratedskel->start("./Integratedskel", arguments);
}

void MDLGUI::RunBackboneExtract1()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Integratedskel completed in " + QString::number(secondsElapsed) + " seconds.\n");
  this->BackboneCandidateFile =  this->BackboneFile.dir().absolutePath() +
    "/BackboneCandidate.vtk";
  this->EdgeRange = this->EdgeRangeInput->text();
  this->MorphStrength = this->MorphStrengthInput->text();

  QStringList arguments;
  arguments << this->DataDir << "out.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->MorphStrength
            << this->BackboneCandidateFile << "1";

  this->OutputDisplay->append("Executing BackboneExtract (first pass)\n");
  this->Time.restart();
  this->BackboneExtract1->start("./BackboneExtract", arguments);
}

void MDLGUI::RunMDABasedSpineExtraction1()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "BackboneExtract (first pass) completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->MDLFeatureFile = this->BackboneFile.dir().absolutePath() +
    "/MDLFeature.txt";
  this->SpineCandidateFile = this->BackboneFile.dir().absolutePath() +
    "/SpineCandidate.vtk";
  this->GraphPruneSize = this->GraphPruneSizeInput->text();
  this->MorphStrength = this->MorphStrengthInput->text();
  this->WeightFactor = this->WeightFactorInput->text();

  QStringList arguments;
  arguments << this->DataDir << "out.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->GraphPruneSize << this->MorphStrength
            << this->WeightFactor <<  this->AnisoDiffusedFile
            << this->MDLFeatureFile << this->SpineCandidateFile;
  this->OutputDisplay->append("Executing MDABasedSpineExtraction (first pass)\n");
  this->Time.restart();
  this->MDABasedSpineExtraction1->start("./MDABasedSpineExtraction", arguments);
}

void MDLGUI::RunBSplineFitting()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "MDABasedSpineExtraction (first pass) completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->SmoothBackboneFile = this->BackboneFile.dir().absolutePath() +
    "/SmoothBackbone.vtk";
  this->ExtraSpineFile = this->BackboneFile.dir().absolutePath() +
    "/ExtraSpine.vtk";

  QStringList arguments;
  arguments << this->ComponentsConnectedFile << this->BackboneCandidateFile
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->SmoothBackboneFile << this->ExtraSpineFile;

  this->OutputDisplay->append("Executing BSplineFitting\n");
  this->Time.restart();
  this->BSplineFitting->start("./BSplineFitting", arguments);
}

void MDLGUI::RunRefiningSkeleton1()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "BSplineFitting completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->RefinedSkeletonFile = this->BackboneFile.dir().absolutePath() +
    "/RefinedSkel.skel";

  QStringList arguments;
  arguments << this->DataDir << "SmoothBackbone.vtk" << "SpineCandidate.vtk"
            << "ExtraSpine.vtk" << "RefinedSkel.skel" << "1";
  this->OutputDisplay->append("Executing RefiningSkeleton (first pass)\n");
  this->Time.restart();
  this->RefiningSkeleton1->start("./RefiningSkeleton", arguments);
}

void MDLGUI::RunBackboneExtract2()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "RefiningSkeleton (first pass) completed in " + QString::number(secondsElapsed) + " seconds.\n");

  QStringList arguments;
  arguments << this->DataDir << "RefinedSkel.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->MorphStrength
            << this->BackboneFile.absoluteFilePath() << "0";
  this->OutputDisplay->append("Executing BackboneExtract (second pass)\n");
  this->Time.restart();
  this->BackboneExtract2->start("./BackboneExtract", arguments);
}

void MDLGUI::RunRefiningSkeleton2()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "BackboneExtract (second pass) completed in " + QString::number(secondsElapsed) + " seconds.\n");

  QStringList arguments;
  arguments << this->DataDir << this->BackboneFile.fileName() 
            << "SpineCandidate.vtk" <<  "ExtraSpine.vtk" << "RefinedSkel.skel"
            << "0";
  this->OutputDisplay->append("Executing RefiningSkeleton (second pass)\n");
  this->Time.restart();
  this->RefiningSkeleton2->start("./RefiningSkeleton", arguments);
}

void MDLGUI::RunMDABasedSpineExtraction2()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "RefiningSkeleton (second pass) completed in " + QString::number(secondsElapsed) + " seconds.\n");

  QStringList arguments;
  arguments << this->DataDir << "RefinedSkel.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->GraphPruneSize << this->MorphStrength
            << this->WeightFactor <<  this->AnisoDiffusedFile
            << this->MDLFeatureFile << this->SpinesFile.absoluteFilePath();
  this->OutputDisplay->append("Executing MDABasedSpineExtraction (second pass)\n");
  this->Time.restart();
  this->MDABasedSpineExtraction2->start("./MDABasedSpineExtraction", arguments);
}

void MDLGUI::DeleteIntermediaryFiles()
{
  QFile f1(this->VolumeProcessedFile);
  if(f1.exists())
    {
    f1.remove();
    }
  QFile f2(this->ComponentsConnectedFile);
  if(f2.exists())
    {
    f2.remove();
    }
  QFile f3(this->AnisoDiffusedFile);
  if(f3.exists())
    {
    f3.remove();
    }
  QFile f4(this->VectorFile);
  if(f4.exists())
    {
    f4.remove();
    }
  QFile f5(this->SeedFile);
  if(f5.exists())
    {
    f5.remove();
    }
  QFile f6(this->SkeletonFile);
  if(f6.exists())
    {
    f6.remove();
    }
  QFile f7(this->BackboneCandidateFile);
  if(f7.exists())
    {
    f7.remove();
    }
  QFile f8(this->SpineCandidateFile);
  if(f8.exists())
    {
    f8.remove();
    }
  QFile f9(this->SmoothBackboneFile);
  if(f9.exists())
    {
    f9.remove();
    }
  QFile f10(this->ExtraSpineFile);
  if(f10.exists())
    {
    f10.remove();
    }
  QFile f11(this->RefinedSkeletonFile);
  if(f11.exists())
    {
    f11.remove();
    }
}

void MDLGUI::FinishedRunningSkeletonization()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Final step completed in " + QString::number(secondsElapsed) + " seconds.\n");
  this->OutputDisplay->append("Skeletonization complete.\n\n");
  this->RunSkeletonizationButton->setEnabled(true);
  if(this->DeleteIntermediaryFilesCheckBox->isChecked())
    {
    this->DeleteIntermediaryFiles();
    }
}

void MDLGUI::closeEvent(QCloseEvent *event)
{
  event->accept();
}

