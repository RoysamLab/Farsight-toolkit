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

#include "vtkColorTransferFunction.h"
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include "vtkVolumeProperty.h"

#include "MDLWizard.h"

//this has to be after MDLWizard.h for some reason...
#include <QVTKWidget.h>

//-----------------------------------------------------------------------------
MDLWizard::MDLWizard()
{
  //create QProcesses which do the actual work 
  this->VolumeProcess = new QProcess(this);
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

  //these executables live in a different place depending on OS
  //find their location now
  this->ExecutablePath = ".";
  QString exeName; 
  #ifdef WIN32
    exeName = "/volumeProcess.exe";
  #else
    exeName = "/volumeProcess";
  #endif
  if(!QFileInfo(this->ExecutablePath + exeName).exists())
    {
    this->ExecutablePath = QDir::currentPath() + "/MDLWizard.app/Contents/bin";
    if(!QFileInfo(this->ExecutablePath + exeName).exists())
      {
      this->ExecutablePath = "../../..";
      if(!QFileInfo(this->ExecutablePath + exeName).exists())
        {
        QMessageBox msgBox;
        msgBox.setText(
          "BUG: can't find required executables.  This isn't going to work...");
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        }
      }
    }

  //don't run automatically by default
  this->InteractiveExecution = true;

  //use a signal mapper to connect the many processes to a single output
  //window
  this->Mapper = new QSignalMapper(this);

  //set the Wizard's name and background image
  this->setObjectName("SkeletonizationWizard");
  this->setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/banner.png"));
  this->setWindowTitle(tr("Skeletonization Wizard"));
  
  //initialize this object
  this->Initialize();
}

//-----------------------------------------------------------------------------
MDLWizard::~MDLWizard()
{
  delete this->OutputWindow;
  this->OutputWindow = 0;
  delete this->HelpWindow;
  this->HelpWindow = 0;
  delete this->RenderWidget;
  this->RenderWidget = 0;
}

//-----------------------------------------------------------------------------
void MDLWizard::Initialize()
{
  //setup the user interface
  this->setupUi(this);
  this->setAttribute(Qt::WA_DeleteOnClose);
  this->move(40, 59);
  this->setWindowTitle(tr("Skeletonization Wizard"));

  //add a help button to display the help window
  this->setOption(QWizard::HaveHelpButton, true);
  this->connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
  this->HelpWindow = new QTextEdit();
  this->HelpWindow->resize(400, 400);
  this->HelpWindow->move(839, 504);
  this->HelpWindow->setWindowTitle("Help");
  this->HelpWindow->hide();

  //set up the output display: a text window that shows the outputs & errors
  //from the various running processes
  this->OutputWindow = new QTextEdit();
  this->OutputWindow->resize(400, 400);
  this->OutputWindow->move(39, 504);
  this->OutputWindow->setWindowTitle("Output Display");
  this->OutputWindow->show();


  //set up the render widget: a QVTKWidget to show intermediary results
  this->RenderWidget = new QVTKWidget();
  this->Renderer = vtkSmartPointer<vtkRenderer>::New();
  this->Renderer->SetBackground(0.33, 0.35, 0.43);
  this->RenderWidget->GetRenderWindow()->AddRenderer(this->Renderer);
  this->RenderWidget->resize(400, 400);
  this->RenderWidget->move(439, 504);
  this->RenderWidget->setWindowTitle("Render Window");
  this->RenderWidget->show();

  //initialize properties used to control the flow through this wizard 
  this->PreprocessingDone = false;
  this->PhaseOneDone = false;
  this->PhaseTwoDone = false;

  //connect the signals and slots of the various QWidgets used by this class
  this->SetupSignalsAndSlots();
}

//-----------------------------------------------------------------------------
void MDLWizard::SetupSignalsAndSlots()
{
  this->connect(this, SIGNAL(currentIdChanged(int)),
                this, SLOT(UpdateHelpWindow()));
  this->connect(this->SelectInputButton, SIGNAL(clicked()),
                this, SLOT(SelectInputImage()));
  this->connect(this->SelectBackboneButton, SIGNAL(clicked()),
                this, SLOT(SelectBackboneFile()));
  this->connect(this->SelectSpinesButton, SIGNAL(clicked()),
                this, SLOT(SelectSpinesFile()));
  this->connect(this->DeleteFilesButton, SIGNAL(clicked()),
                this, SLOT(DeleteIntermediaryFiles()));

  //all these connections are for updating GUI elements that appear on more
  //than one page.
  this->connect(this->EdgeRangeInput1, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateEdgeRangeFrom1()));
  this->connect(this->EdgeRangeInput2, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateEdgeRangeFrom2()));
  this->connect(this->MorphStrengthInput1, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateMorphStrengthFrom1()));
  this->connect(this->MorphStrengthInput2, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateMorphStrengthFrom2()));
  this->connect(this->WeightFactorInput1, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateWeightFactorFrom1()));
  this->connect(this->WeightFactorInput2, SIGNAL(textChanged(const QString &)),
                this, SLOT(UpdateWeightFactorFrom2()));
 
  //connect each process to the output display using a signal mapper.
  connect(this->Mapper, SIGNAL(mapped(QObject *)),
          this, SLOT(AppendOutputToDisplay(QObject *)));
  connect(this->VolumeProcess, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->VolumeProcess, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->VolumeProcess, this->VolumeProcess);
  connect(this->ConnCompntwFldfill, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->ConnCompntwFldfill, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->ConnCompntwFldfill, this->ConnCompntwFldfill);
  connect(this->AnisoDiffuse, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->AnisoDiffuse, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->AnisoDiffuse, this->AnisoDiffuse);
  connect(this->GradientVecField, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->GradientVecField, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->GradientVecField, this->GradientVecField);
  connect(this->Integratedskel, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->Integratedskel, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->Integratedskel, this->Integratedskel);
  connect(this->BackboneExtract1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->BackboneExtract1, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->BackboneExtract1, this->BackboneExtract1);
  connect(this->MDABasedSpineExtraction1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->MDABasedSpineExtraction1, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->MDABasedSpineExtraction1,
                           this->MDABasedSpineExtraction1);
  connect(this->BSplineFitting, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->BSplineFitting, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->BSplineFitting, this->BSplineFitting);
  connect(this->RefiningSkeleton1, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->RefiningSkeleton1, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->RefiningSkeleton1, this->RefiningSkeleton1);
  connect(this->BackboneExtract2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->BackboneExtract2, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->BackboneExtract2, this->BackboneExtract2);
  connect(this->RefiningSkeleton2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->RefiningSkeleton2, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->RefiningSkeleton2, this->RefiningSkeleton2);
  connect(this->MDABasedSpineExtraction2, SIGNAL(readyReadStandardOutput()),
          this->Mapper, SLOT(map()));
  connect(this->MDABasedSpineExtraction2, SIGNAL(readyReadStandardError()),
          this->Mapper, SLOT(map()));
  this->Mapper->setMapping(this->MDABasedSpineExtraction2,
                           this->MDABasedSpineExtraction2);

  //connect the processes to their "run" buttons
  this->connect(this->VolumeProcessButton, SIGNAL(clicked()),
                this, SLOT(RunVolumeProcess()));
  this->connect(this->ConnCompntButton, SIGNAL(clicked()),
                this, SLOT(RunConnCompntwFldfill()));
  this->connect(this->AnisoDiffuseButton, SIGNAL(clicked()),
                this, SLOT(RunAnisoDiffuse()));
  this->connect(this->GradientVecFieldButton, SIGNAL(clicked()),
                this, SLOT(RunGradientVecField()));
  this->connect(this->IntegratedskelButton, SIGNAL(clicked()),
                this, SLOT(RunIntegratedskel()));
  this->connect(this->BackboneExtractButton1, SIGNAL(clicked()),
                this, SLOT(RunBackboneExtract1()));
  this->connect(this->SpineExtractionButton1, SIGNAL(clicked()),
                this, SLOT(RunMDABasedSpineExtraction1()));
  this->connect(this->BSplineFittingButton, SIGNAL(clicked()),
                this, SLOT(RunBSplineFitting()));
  this->connect(this->RefiningSkeletonButton1, SIGNAL(clicked()),
                this, SLOT(RunRefiningSkeleton1()));
  this->connect(this->BackboneExtractButton2, SIGNAL(clicked()),
                this, SLOT(RunBackboneExtract2()));
  this->connect(this->RefiningSkeletonButton2, SIGNAL(clicked()),
                this, SLOT(RunRefiningSkeleton2()));
  this->connect(this->SpineExtractionButton2, SIGNAL(clicked()),
                this, SLOT(RunMDABasedSpineExtraction2()));

  //when a process finishes, run a method that re-enables the appropriate "run"
  //button and displays the results
  this->connect(this->VolumeProcess,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(VolumeProcessFinished()));
  this->connect(this->ConnCompntwFldfill,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(ConnCompntwFldfillFinished()));
  this->connect(this->AnisoDiffuse,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(AnisoDiffuseFinished()));
  this->connect(this->GradientVecField,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(GradientVecFieldFinished()));
  this->connect(this->Integratedskel,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(IntegratedskelFinished()));
  this->connect(this->BackboneExtract1,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(BackboneExtract1Finished()));
  this->connect(this->MDABasedSpineExtraction1,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(MDABasedSpineExtraction1Finished()));
  this->connect(this->BSplineFitting,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(BSplineFittingFinished()));
  this->connect(this->RefiningSkeleton1,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RefiningSkeleton1Finished()));
  this->connect(this->BackboneExtract2,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(BackboneExtract2Finished()));
  this->connect(this->RefiningSkeleton2,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(RefiningSkeleton2Finished()));
  this->connect(this->MDABasedSpineExtraction2,
                SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(MDABasedSpineExtraction2Finished()));
}

//-----------------------------------------------------------------------------
void MDLWizard::SelectInputImage()
{
  this->InputFileName = QFileDialog::getOpenFileName(this,
    tr("Select Input Image"), "", tr("Image Files (*.*)"));
  this->InputFile = QFileInfo(this->InputFileName);
  this->InputImageLabel->setText(this->InputFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::SelectBackboneFile()
{
  this->BackboneFileName = QFileDialog::getSaveFileName(this,
    tr("Select Backbone File"), "", tr("VTK files (*.vtk)"));
  this->BackboneFile = QFileInfo(this->BackboneFileName);
  this->DataDir = this->BackboneFile.dir().absolutePath() + "/";
  this->BackboneOutputLabel->setText(this->BackboneFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::SelectSpinesFile()
{
  this->SpinesFileName = QFileDialog::getSaveFileName(this,
    tr("Select Spines File"), "", tr("VTK files (*.vtk)"));
  this->SpinesFile = QFileInfo(this->SpinesFileName);
  this->SpinesOutputLabel->setText(this->SpinesFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::SetInputImage(const char *filename)
{
  this->InputFileName = filename;
  this->InputFile = QFileInfo(this->InputFileName);
  if(!this->InputFile.exists())
    {
    cerr << "ERROR: input file " << filename << " doesn't exist!" << endl;
    }
  this->InputImageLabel->setText(this->InputFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::SetBackboneFile(const char *filename)
{
  this->BackboneFileName = filename;
  this->BackboneFile = QFileInfo(this->BackboneFileName);
  this->DataDir = this->BackboneFile.dir().absolutePath() + "/";
  this->BackboneOutputLabel->setText(this->BackboneFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::SetSpinesFile(const char *filename)
{
  this->SpinesFileName = filename;
  this->SpinesFile = QFileInfo(this->SpinesFileName);
  this->SpinesOutputLabel->setText(this->SpinesFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateEdgeRangeFrom1()
{
  cout << "I got called" << endl;
  if(this->EdgeRangeInput1->text() != "")
    {
    this->EdgeRange = this->EdgeRangeInput1->text();
    this->EdgeRangeInput2->setText(this->EdgeRange);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateEdgeRangeFrom2()
{
  if(this->EdgeRangeInput2->text() != "")
    {
    this->EdgeRange = this->EdgeRangeInput2->text();
    this->EdgeRangeInput1->setText(this->EdgeRange);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateMorphStrengthFrom1()
{
  if(this->MorphStrengthInput1->text() != "")
    {
    this->MorphStrength = this->MorphStrengthInput1->text();
    this->MorphStrengthInput2->setText(this->MorphStrength);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateMorphStrengthFrom2()
{
  if(this->MorphStrengthInput2->text() != "")
    {
    this->MorphStrength = this->MorphStrengthInput2->text();
    this->MorphStrengthInput1->setText(this->MorphStrength);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateWeightFactorFrom1()
{
  if(this->WeightFactorInput1->text() != "")
    {
    this->WeightFactor = this->WeightFactorInput1->text();
    this->WeightFactorInput2->setText(this->WeightFactor);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateWeightFactorFrom2()
{
  if(this->WeightFactorInput2->text() != "")
    {
    this->WeightFactor = this->WeightFactorInput2->text();
    this->WeightFactorInput1->setText(this->WeightFactor);
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::AppendOutputToDisplay(QObject *mapped)
{
  QProcess *process = static_cast<QProcess *>(mapped);
  if(this->OutputWindow->isHidden())
    {
    this->OutputWindow->show();
    }
  this->OutputWindow->append(QString(process->readAllStandardError()));
  this->OutputWindow->append(QString(process->readAllStandardOutput()));
}

//-----------------------------------------------------------------------------
void MDLWizard::RenderVolume(QFileInfo volumeFile)
{
  //remove all old actors
  this->Renderer->RemoveAllViewProps();
  
  //change the window title so the user knows what file is being displayed
  this->RenderWidget->setWindowTitle(volumeFile.fileName());

  vtkSmartPointer<vtkVolume> volume = 
    this->ConvertRawToVolume(volumeFile.filePath().toStdString().c_str());
  this->Renderer->AddVolume(volume);
  this->Renderer->ResetCamera();
  this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
void MDLWizard::RenderPolyData(QFileInfo polyDataFile)
{
  //remove all old actors
  this->Renderer->RemoveAllViewProps();

  //change the window title so the user knows what file is being displayed
  this->RenderWidget->setWindowTitle(polyDataFile.fileName());

  vtkSmartPointer<vtkActor> actor =
    this->CreateActorFromPolyDataFile(
      polyDataFile.filePath().toStdString().c_str());
  this->Renderer->AddActor(actor);
  this->Renderer->ResetCamera();
  this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
void MDLWizard::RenderPoints(QString pointsFileName)
{
  //remove all old actors
  this->Renderer->RemoveAllViewProps();

  //change the window title so the user knows what file is being displayed
  QFileInfo info(pointsFileName);
  this->RenderWidget->setWindowTitle(info.fileName());

  vtkSmartPointer<vtkActor> actor =
    this->CreateActorFromPointsFile(pointsFileName);
  this->Renderer->AddActor(actor);

  //also render AnisoDiffused.raw so the points have some context
  vtkSmartPointer<vtkVolume> volume = 
    this->ConvertRawToVolume(this->AnisoDiffusedFile.toStdString().c_str());
  this->Renderer->AddVolume(volume);

  this->Renderer->ResetCamera();
  this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkActor> MDLWizard::CreateActorFromPointsFile(
  QString pointsFileName)
{
  //declare our return value
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();

  //open up the file and see how many lines (points) there are
  int numPoints = 0;
  QFile pointsFile(pointsFileName);
  if (!pointsFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
    cerr << "Error opening " << pointsFileName.toStdString() << endl;
    return actor;
    }
  while (!pointsFile.atEnd())
    {
    QByteArray notUsed = pointsFile.readLine();
    notUsed.squeeze();
    numPoints++;
    }

  cout << pointsFileName.toStdString() << " has " << numPoints << " lines." << endl;

  if(!pointsFile.seek(0))
    {
    cerr << "Error rewinding " << pointsFileName.toStdString() << endl;
    return actor;
    }

  //initialize vtkPoints object for visualization
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->Allocate(numPoints);

  //read the points from the file and insert them into the vtkPoints object
  for(int i = 0; i < numPoints; i++)
    {
    QString line(pointsFile.readLine());
    QStringList list = line.split(" ");
    double x = list[0].toDouble();
    double y = list[1].toDouble();
    double z = list[2].toDouble();
    points->InsertPoint(i, x, y, z);
    }

  //finish setting up the actor
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);

  vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
  glyph->SetInput(polyData);
  vtkSmartPointer<vtkGlyphSource2D> glyphSource =
    vtkSmartPointer<vtkGlyphSource2D>::New();
  glyphSource->SetGlyphTypeToCross();
  glyphSource->SetScale(5.0);
  glyphSource->Update();
  glyph->SetInputConnection(1, glyphSource->GetOutputPort());

  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(glyph->GetOutput());

  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(1.0, 1.0, 0.0);

  return actor;
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkVolume> MDLWizard::ConvertRawToVolume(const char *filename)
{
  int sizeX = this->ImageSizeX.toInt();
  int sizeY = this->ImageSizeY.toInt();
  int sizeZ = this->ImageSizeZ.toInt();

  unsigned char *buf =
    new unsigned char[sizeX * sizeY * sizeZ];
  unsigned char *p = buf;
  FILE *fp = fopen(filename, "r"); 
  fread(buf, sizeX * sizeY, sizeZ, fp);
  fclose(fp);

  vtkSmartPointer<vtkImageData> imageData =
    vtkSmartPointer<vtkImageData>::New();
  imageData->SetDimensions(sizeX, sizeY, sizeZ);
  imageData->SetSpacing(1.0,1.0,1.0);
  imageData->SetOrigin(0.0,0.0,0.0);
  imageData->SetScalarTypeToUnsignedChar(); 
  imageData->SetNumberOfScalarComponents(1); 
  imageData->AllocateScalars(); 

  int* dims = imageData->GetDimensions();
  unsigned char *a = (unsigned char*)(imageData->GetScalarPointer());
  for (int z=0; z<dims[2]; z++)
    {
    for (int y=0; y<dims[1]; y++)
      {
      for (int x=0; x<dims[0]; x++)
        {
        *a = *p;
        a++;
        p++;
        }
      }
    }

  vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> mapper =
    vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
  mapper->SetInput(imageData);

  vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
  volume->SetMapper(mapper);

  double range[2];
  imageData->GetScalarRange(range);

  vtkSmartPointer<vtkVolumeProperty> volumeProperty =
    this->NewRGBVolumeProperty(range);
  volume->SetProperty(volumeProperty);

  delete [] buf;
  return volume;
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkActor> MDLWizard::CreateActorFromPolyDataFile(
  const char *filename)
{
  vtkSmartPointer<vtkPolyDataReader> reader =
    vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(filename);

  vtkPolyData *polyData = reader->GetOutput();

  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(polyData);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  return actor;
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkVolumeProperty> MDLWizard::NewRGBVolumeProperty(
  const double range[])
{
  // Create transfer mapping scalar value to opacity.
  vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction =
    vtkSmartPointer<vtkPiecewiseFunction>::New();
  opacityTransferFunction->AddPoint( range[0], 0.0);
  opacityTransferFunction->AddPoint(range[1], 0.5);

  // Create transfer mapping scalar value to color.
  vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction =
    vtkSmartPointer<vtkColorTransferFunction>::New();
  colorTransferFunction->SetColorSpaceToRGB();
  //these hardcoded values aren't gonna work in the long run
  //talk to the ParaView folks about how they manage this problem...
  colorTransferFunction->AddRGBPoint(range[0], 0.0, 0.0, 1.0);
  colorTransferFunction->AddRGBPoint(range[1], 1.0, 0.0, 0.0);

  vtkSmartPointer<vtkVolumeProperty> volumeProperty =
    vtkSmartPointer<vtkVolumeProperty>::New();
  volumeProperty->SetColor(colorTransferFunction);
  volumeProperty->SetScalarOpacity(opacityTransferFunction);
  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationTypeToLinear();
  volumeProperty->SetScalarOpacityUnitDistance(0.75);

  return volumeProperty;
}


//-----------------------------------------------------------------------------
void MDLWizard::RunVolumeProcess()
{
  this->VolumeProcessButton->setEnabled(false);
  this->OutputWindow->clear();

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

  //run the executable
  this->VolumeProcess->start(this->ExecutablePath + "/volumeProcess",
                             arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunConnCompntwFldfill()
{
  this->VolumeProcessButton->setEnabled(false);
  this->ConnCompntButton->setEnabled(false);
  this->ConnectedComponentsSize = this->ComponentsSizeInput->text();
  this->OutputWindow->clear();

  //construct a list of arguments
  this->ComponentsConnectedFile = 
    this->BackboneFile.dir().absolutePath() + "/components_Connected.raw";
  QStringList arguments;
  arguments << this->VolumeProcessedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ
            << this->ComponentsConnectedFile << this->ConnectedComponentsSize;
  
  this->ConnCompntwFldfill->start(this->ExecutablePath + "/ConnCompntwFldfill",
                                  arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunAnisoDiffuse()
{
  this->VolumeProcessButton->setEnabled(false);
  this->ConnCompntButton->setEnabled(false);
  this->AnisoDiffuseButton->setEnabled(false);
  this->OutputWindow->clear();

  this->AnisoDiffusedFile =
    this->BackboneFile.dir().absolutePath() + "/Aniso_Diffused.raw";

  QStringList arguments;
  arguments << this->ComponentsConnectedFile << this->ImageSizeX
            << this->ImageSizeY << this->ImageSizeZ << this->AnisoDiffusedFile;

  this->AnisoDiffuse->start(this->ExecutablePath + "/AnisoDiffuse", arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunGradientVecField()
{
  this->VolumeProcessButton->setEnabled(false);
  this->ConnCompntButton->setEnabled(false);
  this->AnisoDiffuseButton->setEnabled(false);
  this->GradientVecFieldButton->setEnabled(false);
  this->OutputWindow->clear();

  this->VectorFile = this->BackboneFile.dir().absolutePath() + "/out.vec";

  QStringList arguments;
  arguments << this->AnisoDiffusedFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << this->VectorFile;

  this->GradientVecField->start(this->ExecutablePath + "/GradientVecField",
                                arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunIntegratedskel()
{
  this->IntegratedskelButton->setEnabled(false);
  this->OutputWindow->clear();

  this->SeedFile = this->BackboneFile.dir().absolutePath() + "/out.seed";
  this->SkeletonFile = this->BackboneFile.dir().absolutePath() + "/out.skel";
  this->VectorMagnitude = this->VectorMagnitudeInput->text();

  QStringList arguments;
  arguments << this->VectorFile << this->ImageSizeX << this->ImageSizeY
            << this->ImageSizeZ << this->VectorMagnitude << this->SeedFile
            << this->SkeletonFile;

  this->Integratedskel->start(this->ExecutablePath + "/Integratedskel",
                              arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunBackboneExtract1()
{
  this->IntegratedskelButton->setEnabled(false);
  this->BackboneExtractButton1->setEnabled(false);
  this->OutputWindow->clear();

  this->BackboneCandidateFile =  this->BackboneFile.dir().absolutePath() +
    "/BackboneCandidate.vtk";
  this->EdgeRange = this->EdgeRangeInput1->text();
  this->MorphStrength = this->MorphStrengthInput1->text();

  QStringList arguments;
  arguments << this->DataDir << "out.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->MorphStrength
            << this->BackboneCandidateFile << "1";

  this->BackboneExtract1->start(this->ExecutablePath + "/BackboneExtract",
                                arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunMDABasedSpineExtraction1()
{
  this->IntegratedskelButton->setEnabled(false);
  this->BackboneExtractButton1->setEnabled(false);
  this->SpineExtractionButton1->setEnabled(false);
  this->OutputWindow->clear();

  this->MDLFeatureFile = this->BackboneFile.dir().absolutePath() +
    "/MDLFeature.txt";
  this->SpineCandidateFile = this->BackboneFile.dir().absolutePath() +
    "/SpineCandidate.vtk";
  this->GraphPruneSize = this->GraphPruneSizeInput->text();
  this->MorphStrength = this->MorphStrengthInput2->text();
  this->EdgeRange = this->EdgeRangeInput2->text();
  this->WeightFactor = this->WeightFactorInput1->text();

  QStringList arguments;
  arguments << this->DataDir << "out.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->GraphPruneSize << this->MorphStrength
            << this->WeightFactor <<  this->AnisoDiffusedFile
            << this->MDLFeatureFile << this->SpineCandidateFile;
  
  this->MDABasedSpineExtraction1->start(
    this->ExecutablePath + "/MDABasedSpineExtraction", arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunBSplineFitting()
{
  this->IntegratedskelButton->setEnabled(false);
  this->BackboneExtractButton1->setEnabled(false);
  this->SpineExtractionButton1->setEnabled(false);
  this->BSplineFittingButton->setEnabled(false);
  this->OutputWindow->clear();

  this->SmoothBackboneFile = this->BackboneFile.dir().absolutePath() +
    "/SmoothBackbone.vtk";
  this->ExtraSpineFile = this->BackboneFile.dir().absolutePath() +
    "/ExtraSpine.vtk";

  this->BSplineOrder = this->OrderInput->text();
  this->BSplineLevels = this->LevelsInput->text();

  QStringList arguments;
  arguments << this->ComponentsConnectedFile << this->BackboneCandidateFile
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->BSplineOrder << this->BSplineLevels
            << this->SmoothBackboneFile << this->ExtraSpineFile;

  this->BSplineFitting->start(this->ExecutablePath + "/BSplineFitting",
                              arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunRefiningSkeleton1()
{
  this->IntegratedskelButton->setEnabled(false);
  this->BackboneExtractButton1->setEnabled(false);
  this->SpineExtractionButton1->setEnabled(false);
  this->BSplineFittingButton->setEnabled(false);
  this->RefiningSkeletonButton1->setEnabled(false);
  this->OutputWindow->clear();

  this->RefinedSeedFile = this->BackboneFile.dir().absolutePath() +
    "/RefinedSeed.seed";
  this->RefinedSkeletonFile = this->BackboneFile.dir().absolutePath() +
    "/RefinedSkel.skel";

  QStringList arguments;
  arguments << this->DataDir << "SmoothBackbone.vtk" << "SpineCandidate.vtk"
            << "ExtraSpine.vtk" << "RefinedSeed.seed" << "RefinedSkel.skel"
            << "1";
  
  this->RefiningSkeleton1->start(this->ExecutablePath + "/RefiningSkeleton",
                                 arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunBackboneExtract2()
{
  this->BackboneExtractButton2->setEnabled(false);
  this->OutputWindow->clear();

  this->MorphStrength = this->MorphStrengthInput2->text();
  this->EdgeRange = this->EdgeRangeInput2->text();

  QStringList arguments;
  arguments << this->DataDir << "RefinedSkel.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->MorphStrength
            << this->BackboneFile.absoluteFilePath() << "0";
  
  this->BackboneExtract2->start(this->ExecutablePath + "/BackboneExtract",
                                arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunRefiningSkeleton2()
{
  this->BackboneExtractButton2->setEnabled(false);
  this->RefiningSkeletonButton2->setEnabled(false);
  this->OutputWindow->clear();

  QStringList arguments;
  arguments << this->DataDir << this->BackboneFile.fileName() 
            << "SpineCandidate.vtk" << "ExtraSpine.vtk" << "RefinedSeed.seed"
            << "RefinedSkel.skel" << "0";
 
  this->RefiningSkeleton2->start(this->ExecutablePath + "/RefiningSkeleton",
                                 arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::RunMDABasedSpineExtraction2()
{
  this->BackboneExtractButton2->setEnabled(false);
  this->RefiningSkeletonButton2->setEnabled(false);
  this->SpineExtractionButton2->setEnabled(false);
  this->OutputWindow->clear();
  
  this->MorphStrength = this->MorphStrengthInput2->text();
  this->EdgeRange = this->EdgeRangeInput2->text();
  this->GraphPruneSize = this->GraphPruneSizeInput->text();
  this->WeightFactor = this->WeightFactorInput2->text();

  QStringList arguments;
  arguments << this->DataDir << "RefinedSkel.skel" << "components_Connected.raw"
            << this->ImageSizeX << this->ImageSizeY << this->ImageSizeZ 
            << this->EdgeRange << this->GraphPruneSize << this->MorphStrength
            << this->WeightFactor <<  this->AnisoDiffusedFile
            << this->MDLFeatureFile << this->SpinesFile.absoluteFilePath();
  
  this->MDABasedSpineExtraction2->start(
    this->ExecutablePath + "/MDABasedSpineExtraction", arguments);
}

//-----------------------------------------------------------------------------
void MDLWizard::VolumeProcessFinished()
{
  //we have to grab the image size out of the output if our input wasn't a
  //raw image.
  if(!(this->RawInput))
    {
    QRegExp rx("input image size: \\[(\\d+), (\\d+), (\\d+)\\]"); 
    if(this->OutputWindow->document()->toPlainText().indexOf(rx) != -1)
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
  this->VolumeProcessButton->setEnabled(true);
  this->ConnCompntButton->setEnabled(true);
  this->ConnCompntButton->setFocus();
  this->RenderVolume(QFileInfo(this->VolumeProcessedFile));
  if(!this->InteractiveExecution)
    {
    this->RunConnCompntwFldfill();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::ConnCompntwFldfillFinished()
{
  this->VolumeProcessButton->setEnabled(true);
  this->ConnCompntButton->setEnabled(true);
  this->AnisoDiffuseButton->setEnabled(true);
  this->AnisoDiffuseButton->setFocus();
  this->RenderVolume(QFileInfo(this->ComponentsConnectedFile));
  if(!this->InteractiveExecution)
    {
    this->RunAnisoDiffuse();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::AnisoDiffuseFinished()
{
  this->VolumeProcessButton->setEnabled(true);
  this->ConnCompntButton->setEnabled(true);
  this->AnisoDiffuseButton->setEnabled(true);
  this->GradientVecFieldButton->setEnabled(true);
  this->GradientVecFieldButton->setFocus();
  this->RenderVolume(QFileInfo(this->AnisoDiffusedFile));
  if(!this->InteractiveExecution)
    {
    this->RunGradientVecField();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::GradientVecFieldFinished()
{
  this->VolumeProcessButton->setEnabled(true);
  this->ConnCompntButton->setEnabled(true);
  this->AnisoDiffuseButton->setEnabled(true);
  this->GradientVecFieldButton->setEnabled(true);
  this->PreprocessingDone = true;
  this->PreprocessingPage->CheckIfComplete();
  if(!this->InteractiveExecution)
    {
    this->next();
    this->RunIntegratedskel();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::IntegratedskelFinished()
{
  this->IntegratedskelButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->BackboneExtractButton1->setFocus();
  this->RenderPoints(this->SeedFile);
  if(!this->InteractiveExecution)
    {
    this->RunBackboneExtract1();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::BackboneExtract1Finished()
{
  this->IntegratedskelButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->SpineExtractionButton1->setEnabled(true);
  this->SpineExtractionButton1->setFocus();
  this->RenderPolyData(QFileInfo(this->BackboneCandidateFile));
  if(!this->InteractiveExecution)
    {
    this->RunMDABasedSpineExtraction1();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::MDABasedSpineExtraction1Finished()
{
  this->IntegratedskelButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->SpineExtractionButton1->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->BSplineFittingButton->setFocus();
  this->RenderPolyData(QFileInfo(this->SpineCandidateFile));
  if(!this->InteractiveExecution)
    {
    this->RunBSplineFitting();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::BSplineFittingFinished()
{
  this->IntegratedskelButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->SpineExtractionButton1->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->RefiningSkeletonButton1->setEnabled(true);
  this->RefiningSkeletonButton1->setFocus();
  this->RenderPolyData(QFileInfo(this->SmoothBackboneFile));
  if(!this->InteractiveExecution)
    {
    this->RunRefiningSkeleton1();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::RefiningSkeleton1Finished()
{
  this->IntegratedskelButton->setEnabled(true);
  this->IntegratedskelButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->SpineExtractionButton1->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->BackboneExtractButton1->setEnabled(true);
  this->SpineExtractionButton1->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->RefiningSkeletonButton1->setEnabled(true);
  this->PhaseOneDone = true;
  this->PhaseOnePage->CheckIfComplete();
  this->RenderPoints(this->RefinedSeedFile);
  if(!this->InteractiveExecution)
    {
    this->next();
    this->RunBackboneExtract2();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::BackboneExtract2Finished()
{
  this->BackboneExtractButton2->setEnabled(true);
  this->RefiningSkeletonButton2->setEnabled(true);
  this->RefiningSkeletonButton2->setFocus();
  this->RenderPolyData(this->BackboneFile);
  if(!this->InteractiveExecution)
    {
    this->RunRefiningSkeleton2();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::RefiningSkeleton2Finished()
{
  this->BackboneExtractButton2->setEnabled(true);
  this->RefiningSkeletonButton2->setEnabled(true);
  this->SpineExtractionButton2->setEnabled(true);
  this->SpineExtractionButton2->setFocus();
  this->RenderPoints(this->RefinedSeedFile);
  if(!this->InteractiveExecution)
    {
    this->RunMDABasedSpineExtraction2();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::MDABasedSpineExtraction2Finished()
{
  this->BackboneExtractButton2->setEnabled(true);
  this->RefiningSkeletonButton2->setEnabled(true);
  this->SpineExtractionButton2->setEnabled(true);
  this->DeleteFilesButton->setEnabled(true);
  this->PhaseTwoDone = true;
  this->PhaseTwoPage->CheckIfComplete();
  this->RenderPolyData(this->SpinesFile);
  if(!this->InteractiveExecution)
    {
    this->DeleteIntermediaryFiles();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::DeleteIntermediaryFiles()
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
  QFile f11(this->RefinedSeedFile);
  if(f11.exists())
    {
    f11.remove();
    }
  QFile f12(this->RefinedSkeletonFile);
  if(f12.exists())
    {
    f12.remove();
    }
  this->OutputWindow->append(QString("Intermediary files deleted."));
  if(!this->InteractiveExecution)
    {
    this->close();
    }
}

//-----------------------------------------------------------------------------
void MDLWizard::UpdateHelpWindow()
{
  if(this->HelpWindow->isVisible())
    {
    this->showHelp();
    }
}


//-----------------------------------------------------------------------------
void MDLWizard::SetConnectedComponentsSize(const char *s)
{
  this->ConnectedComponentsSize = s;
  this->ComponentsSizeInput->setText(this->ConnectedComponentsSize);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetVectorMagnitude(const char *s)
{
  this->VectorMagnitude = s;
  this->VectorMagnitudeInput->setText(this->VectorMagnitude);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetEdgeRange(const char *s)
{
  this->EdgeRange = s;
  this->EdgeRangeInput1->setText(this->EdgeRange);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetGraphPruneSize(const char *s)
{
  this->GraphPruneSize = s;
  this->GraphPruneSizeInput->setText(this->GraphPruneSize);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetMorphStrength(const char *s)
{
  this->MorphStrength = s;
  this->MorphStrengthInput1->setText(this->MorphStrength);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetWeightFactor(const char *s)
{
  this->WeightFactor = s;
  this->WeightFactorInput1->setText(this->WeightFactor);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetBSplineOrder(const char *s)
{
  this->BSplineOrder = s;
  this->OrderInput->setText(this->BSplineOrder);
}

//-----------------------------------------------------------------------------
void MDLWizard::SetBSplineLevels(const char *s)
{
  this->BSplineLevels = s;
  this->LevelsInput->setText(this->BSplineLevels);
}

//-----------------------------------------------------------------------------
void MDLWizard::showHelp()
{
   static QString lastHelpMessage;
   QString message;
   QFile file;

   switch (currentId())
     {
     enum { Page_Intro, Page_Preprocessing, Page_PhaseOne, Page_PhaseTwo };
     case Page_Intro:
       file.setFileName(":/html/intro.html");
       if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
         {
         return;
         }
       this->HelpWindow->setHtml(file.readAll());
       break;
     case Page_Preprocessing:
       file.setFileName(":/html/preprocessing.html");
       if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
         {
         return;
         }
       this->HelpWindow->setHtml(file.readAll());
     break;
     case Page_PhaseOne:
       file.setFileName(":/html/phaseone.html");
       if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
         {
         return;
         }
       this->HelpWindow->setHtml(file.readAll());
     break;
     case Page_PhaseTwo:
       file.setFileName(":/html/phasetwo.html");
       if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
         {
         return;
         }
       this->HelpWindow->setHtml(file.readAll());
     break;
   default:
     this->HelpWindow->setText("This help is likely not to be of any help.");
   }
   this->HelpWindow->show();
}

//-----------------------------------------------------------------------------
void MDLWizard::closeEvent(QCloseEvent *event)
{
  event->accept();
  this->OutputWindow->close();
  this->HelpWindow->close();
}

