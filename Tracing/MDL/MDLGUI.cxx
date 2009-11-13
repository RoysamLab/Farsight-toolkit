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
  this->AnisoDiffuse = new QProcess(this);
  this->BSplineFitting = new QProcess(this);
  this->ConnCompntwFldfill = new QProcess(this);
  this->GradientVecField = new QProcess(this);
  this->Integratedskel = new QProcess(this);
  this->MinSpanTree = new QProcess(this);
  this->volumeProcess = new QProcess(this);
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
  this->connect(this->SelectOutputButton, SIGNAL(pressed()), this,
    SLOT(SelectOutputImage()));
  this->connect(this->SelectSpinesButton, SIGNAL(pressed()), this,
    SLOT(SelectSpinesFile()));
  this->connect(this->ComponentsSizeInput, SIGNAL(textChanged(const QString &)),
    this, SLOT(CheckInputs()));
  this->connect(this->RunSkeletonizationButton, SIGNAL(pressed()), this,
    SLOT(RunSkeletonization()));
  this->connect(this->CombineOutputsButton, SIGNAL(toggled(bool)), this,
    SLOT(ToggleSpinesInputButton()));
  
  //set up the QSignalMapper so that each process can append to the
  //DisplayOutput field
  this->Mapper->setMapping(this->AnisoDiffuse, this->AnisoDiffuse);
  this->Mapper->setMapping(this->BSplineFitting, this->BSplineFitting);
  this->Mapper->setMapping(this->ConnCompntwFldfill, this->ConnCompntwFldfill);
  this->Mapper->setMapping(this->GradientVecField, this->GradientVecField);
  this->Mapper->setMapping(this->Integratedskel, this->Integratedskel);
  this->Mapper->setMapping(this->MinSpanTree, this->MinSpanTree);
  this->Mapper->setMapping(this->volumeProcess, this->volumeProcess);
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
  connect(this->MinSpanTree, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
  connect(this->volumeProcess, SIGNAL(readyReadStandardOutput()), this->Mapper,
    SLOT(map()));
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
                this, SLOT(RunMinSpanTree()));
  this->connect(this->MinSpanTree, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RunBSplineFitting()));
  this->connect(this->BSplineFitting, SIGNAL(finished(int, QProcess::ExitStatus)),
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

void MDLGUI::SelectOutputImage()
{
   QString filePath = QFileDialog::getSaveFileName(this,
     tr("Select Output File"), "", tr("VTK files (*.vtk)"));
   this->OutputFile = QFileInfo(filePath);
   this->OutputImageLabel->setText(this->OutputFile.fileName());
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

void MDLGUI::ToggleSpinesInputButton()
{
  if(this->CombineOutputsButton->isChecked())
    {
    this->SelectSpinesButton->setEnabled(false);
    this->SpinesLabel->hide();
    }
  else
    {
    this->SelectSpinesButton->setEnabled(true);
    this->SpinesLabel->show();
    }
  this->CheckInputs();
}

void MDLGUI::CheckInputs()
{
  if(!this->CombineOutputsButton->isChecked() &&
     this->SpinesFile.filePath() == "")
    {
    this->RunSkeletonizationButton->setEnabled(false);
    return;
    }

  if(this->InputFile.filePath() != "" && this->OutputFile.filePath() != "" &&
     this->ComponentsSizeInput->text() != "")
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

void MDLGUI:: RunvolumeProcess()
{
  //construct the list of arguments
  QStringList arguments;
  this->VolumeProcessedFile =
    this->OutputFile.dir().absolutePath() + "/volume_Processed.raw";
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
  this->OutputDisplay->append("Step #1, executing volumeProcess\n");

  //start a timer
  this->Time.start();

  cout << "*** volume ***" << endl;
  cout << arguments.size() << endl;
  for (int i = 0; i < arguments.size(); ++i)
      cout << arguments.at(i).toLocal8Bit().constData() << endl;

  //run the executable
  this->volumeProcess->start("./volumeProcess", arguments);
}

void MDLGUI:: RunConnCompntwFldfill()
{
  //provide information about how long the last step took.
  int secondsElapsed = this->Time.elapsed() / 1000.0;

  this->OutputDisplay->append(
    "Step #1 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  cout << "*** Compnt ***" << endl;
  if(!(this->RawInput))
    {
    cout << "I go here" << endl;
    //also have to grab image size out of the output
    QRegExp rx("input image size: \\[(\\d+), (\\d+), (\\d+)\\]"); 
    if(this->OutputDisplay->document()->toPlainText().indexOf(rx) != -1)
      {
      QStringList list = rx.capturedTexts();
      this->ImageSizeX = list[1];
      this->ImageSizeY = list[2];
      this->ImageSizeZ = list[3];
      cout << " X == " << this->ImageSizeX.toStdString() << endl;
      cout << " Y == " << this->ImageSizeY.toStdString() << endl;
      cout << " Z == " << this->ImageSizeZ.toStdString() << endl;
      }
    else
      {
      cout << "fux0red..." << endl;
      QMessageBox msgBox;
      msgBox.setText("BUG: Didn't find the image size in output...");
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
      }
    }

  //construct a list of arguments
  this->ComponentsConnectedFile = 
    this->OutputFile.dir().absolutePath() + "/components_Connected.raw";
  QStringList arguments;
  arguments << this->VolumeProcessedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ
            << this->ComponentsConnectedFile << this->ConnectedComponentsSize;
  
  this->OutputDisplay->append("Step #2, executing ConnCompntwFldfill\n");
  this->Time.restart();
  this->ConnCompntwFldfill->start("./ConnCompntwFldfill", arguments);
}

void MDLGUI:: RunAnisoDiffuse()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Step #2 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->AnisoDiffusedFile =
    this->OutputFile.dir().absolutePath() + "/Aniso_Diffused.raw";

  QStringList arguments;
  arguments << this->ComponentsConnectedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ << this->AnisoDiffusedFile;

  this->OutputDisplay->append("Step #3, executing AnisoDiffuse\n");
  this->Time.restart();
  this->AnisoDiffuse->start("./AnisoDiffuse", arguments);
}

void MDLGUI:: RunGradientVecField()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Step #3 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->VectorFile = this->OutputFile.dir().absolutePath() + "/out.vec";

  QStringList arguments;
  arguments << this->AnisoDiffusedFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << this->VectorFile;

  this->OutputDisplay->append("Step #4, executing GradientVecField\n");
  this->Time.restart();
  this->GradientVecField->start("./GradientVecField", arguments);
}

void MDLGUI:: RunIntegratedskel()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Step #4 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  this->SeedFile = this->OutputFile.dir().absolutePath() + "/out.seed";
  this->SkeletonFile = this->OutputFile.dir().absolutePath() + "/out.skel";

  QStringList arguments;
  arguments << this->VectorFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << "0.05" << this->SeedFile
            << this->SkeletonFile;

  this->OutputDisplay->append("Step #5, executing Integratedskel\n");
  this->Time.restart();
  this->Integratedskel->start("./Integratedskel", arguments);
}

void MDLGUI:: RunMinSpanTree()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Step #5 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  QString dataDir = this->OutputFile.dir().absolutePath() + "/";
  this->SeedFile = this->OutputFile.dir().absolutePath() + "/out.seed";
  this->SkeletonFile = this->OutputFile.dir().absolutePath() + "/out.skel";

  QStringList arguments;
  arguments << dataDir << "out.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ
            << "8" << "4" << "170" << "0.7" << this->AnisoDiffusedFile;
  if(this->SmoothBackbonesButton->isChecked())
    {
    this->UnsmoothedBackbonesFile = this->OutputFile.filePath();
    this->UnsmoothedBackbonesFile.replace(QString(".vtk"),
                                         QString("-unsmoothed.vtk"));
    arguments << this->UnsmoothedBackbonesFile;
    if(this->CombineOutputsButton->isChecked())
      {
      arguments << "--combine-outputs";
      }
    else
      {
      this->UnsmoothedSpinesFile = this->SpinesFile.filePath();
      this->UnsmoothedSpinesFile.replace(QString(".vtk"),
                                         QString("-unsmoothed.vtk"));
      arguments << this->UnsmoothedSpinesFile;
      }
    }
  else
    {
    arguments << this->OutputFile.filePath();
    if(this->CombineOutputsButton->isChecked())
      {
      arguments << "--combine-outputs";
      }
    else
      {
      arguments << this->SpinesFile.filePath();
      }
    }
  cout << "*** tree ***" << endl;
  for (int i = 0; i < arguments.size(); ++i)
      cout << arguments.at(i).toLocal8Bit().constData() << endl;

  this->OutputDisplay->append("Step #6, executing MinSpanTree\n");
  this->Time.restart();
  this->MinSpanTree->start("./MinSpanTree", arguments);
}

void MDLGUI:: RunBSplineFitting()
{
  if(!this->SmoothBackbonesButton->isChecked())
    {
    this->FinishedRunningSkeletonization();
    return;
    }

  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Step #6 completed in " + QString::number(secondsElapsed) + " seconds.\n");

  QString dataDir = this->OutputFile.dir().absolutePath() + "/";
  QStringList arguments;
  arguments << this->InputFile.filePath() << this->UnsmoothedBackbonesFile
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->OutputFile.filePath() << this->SpinesFile.filePath();

  this->OutputDisplay->append("Step #7, executing BSplineFitting\n");
  this->Time.restart();
  this->BSplineFitting->start("./BSplineFitting", arguments);
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
  QFile f7(this->BackboneFile);
  if(f7.exists())
    {
    f7.remove();
    }
  QFile f8(this->UnsmoothedBackbonesFile);
  if(f8.exists())
    {
    f8.remove();
    }
  QFile f9(this->UnsmoothedSpinesFile);
  if(f9.exists())
    {
    f9.remove();
    }
}

void MDLGUI::FinishedRunningSkeletonization()
{
  int secondsElapsed = this->Time.elapsed() / 1000.0;
  this->OutputDisplay->append(
    "Last step completed in " + QString::number(secondsElapsed) + " seconds.\n");
  this->RunSkeletonizationButton->setEnabled(true);
  if(this->DeleteIntermediaryFilesCheckBox->isChecked())
    {
    this->DeleteIntermediaryFiles();
    }
  //also need to reset the processes and clear the inputs
}

void MDLGUI::closeEvent(QCloseEvent *event)
{
  event->accept();
}

