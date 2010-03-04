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

#include "MDL2WizardHelper.h"
#include "MDL2Wizard.h"

#include <QVTKWidget.h>

//-----------------------------------------------------------------------------
MDL2Wizard::MDL2Wizard()
{
  //create a helper, which does most of the actual work
  this->Helper = new MDL2WizardHelper();
  this->Helper->start();

  //don't run automatically by default
  this->InteractiveExecution = true;

  //set the Wizard's name and background image
  this->setObjectName("SkeletonizationWizard");
  this->setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/banner.png"));
  this->setWindowTitle(tr("Skeletonization Wizard"));
  
  //initialize this object
  this->Initialize();
}

//-----------------------------------------------------------------------------
MDL2Wizard::~MDL2Wizard()
{
  delete this->OutputWindow;
  this->OutputWindow = 0;
  delete this->HelpWindow;
  this->HelpWindow = 0;
  delete this->RenderWidget;
  this->RenderWidget = 0;
  this->Helper->exit();
  delete this->Helper;
  this->Helper = 0;
}

//-----------------------------------------------------------------------------
void MDL2Wizard::Initialize()
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
void MDL2Wizard::SetupSignalsAndSlots()
{
  this->connect(this, SIGNAL(currentIdChanged(int)),
                this, SLOT(UpdateHelpWindow()));
  this->connect(this->SelectInputButton, SIGNAL(clicked()),
                this, SLOT(SelectInputImage()));
  this->connect(this->SelectBackboneButton, SIGNAL(clicked()),
                this, SLOT(SelectBackboneFile()));
  this->connect(this->SelectSkeletonButton, SIGNAL(clicked()),
                this, SLOT(SelectSkeletonFile()));

  //connect the "run" buttons to their associated functions here
  this->connect(this->MaskUsingGraphCutsButton, SIGNAL(clicked()),
                this, SLOT(MaskUsingGraphCuts()));
  this->connect(this->MaskSmallConnCompButton, SIGNAL(clicked()),
                this, SLOT(MaskSmallConnComp()));
  this->connect(this->IntegratedskelButton, SIGNAL(clicked()),
                this, SLOT(Integratedskel()));
  this->connect(this->CreateGraphAndMSTButton1, SIGNAL(clicked()),
                this, SLOT(CreateGraphAndMST1()));
  this->connect(this->ErodeAndDilateButton1, SIGNAL(clicked()),
                this, SLOT(ErodeAndDilateNodeDegree1()));
  this->connect(this->BackboneExtractButton, SIGNAL(clicked()),
                this, SLOT(BackboneExtract()));
  this->connect(this->BSplineFittingButton, SIGNAL(clicked()),
                this, SLOT(BSplineFitting()));
  this->connect(this->CreateGraphAndMSTButton2, SIGNAL(clicked()),
                this, SLOT(CreateGraphAndMST2()));
  this->connect(this->ErodeAndDilateButton2, SIGNAL(clicked()),
                this, SLOT(ErodeAndDilateNodeDegree2()));
  this->connect(this->SaveOutputButton, SIGNAL(clicked()),
                this, SLOT(SaveOutput()));

  //this object emits signals that tell the helper to do the actual work
  this->connect(this, SIGNAL(InputChanged(std::string)),
                this->Helper, SLOT(ReadImage(std::string)));
  this->connect(this, SIGNAL(StartMaskUsingGraphCuts()),
                this->Helper, SLOT(RunMaskUsingGraphCuts()), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartMaskSmallConnComp(int)),
                this->Helper, SLOT(RunMaskSmallConnComp(int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartIntegratedskel(double)),
                this->Helper, SLOT(RunIntegratedskel(double)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartCreateGraphAndMST1(int)),
                this->Helper, SLOT(RunCreateGraphAndMST1(int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartErodeAndDilateNodeDegree1(int)),
                this->Helper, SLOT(RunErodeAndDilateNodeDegree1(int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartBackboneExtract()),
                this->Helper, SLOT(RunBackboneExtract()), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartBSplineFitting(unsigned int, unsigned int)),
                this->Helper, SLOT(RunBSplineFitting(unsigned int, unsigned int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartCreateGraphAndMST2(int)),
                this->Helper, SLOT(RunCreateGraphAndMST2(int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(StartErodeAndDilateNodeDegree2(int)),
                this->Helper, SLOT(RunErodeAndDilateNodeDegree2(int)), Qt::QueuedConnection);
  this->connect(this, SIGNAL(ReadyToSaveBackbone(const char *)),
                this->Helper, SLOT(WriteBackbone(const char *)));
  this->connect(this, SIGNAL(ReadyToSaveSkeleton(const char *)),
                this->Helper, SLOT(WriteSkeleton(const char *)));

  //when a process finishes, run a method that re-enables the appropriate 
  //buttons and displays the results
  this->connect(this->Helper, SIGNAL(MaskUsingGraphCutsFinished()), 
                this, SLOT(DisplayMaskUsingGraphCutsResults()));
  this->connect(this->Helper, SIGNAL(MaskSmallConnCompFinished()), 
                this, SLOT(DisplayMaskSmallConnCompResults()));
  this->connect(this->Helper, SIGNAL(IntegratedskelFinished()), 
                this, SLOT(DisplayIntegratedskelResults()));
  this->connect(this->Helper, SIGNAL(CreateGraphAndMST1Finished()), 
                this, SLOT(DisplayCreateGraphAndMST1Results()));
  this->connect(this->Helper, SIGNAL(ErodeAndDilateNodeDegree1Finished()), 
                this, SLOT(DisplayErodeAndDilateNodeDegree1Results()));
  this->connect(this->Helper, SIGNAL(BackboneExtractFinished()), 
                this, SLOT(DisplayBackboneExtractResults()));
  this->connect(this->Helper, SIGNAL(BSplineFittingFinished()), 
                this, SLOT(DisplayBSplineFittingResults()));
  this->connect(this->Helper, SIGNAL(CreateGraphAndMST2Finished()), 
                this, SLOT(DisplayCreateGraphAndMST2Results()));
  this->connect(this->Helper, SIGNAL(ErodeAndDilateNodeDegree2Finished()), 
                this, SLOT(DisplayErodeAndDilateNodeDegree2Results()));
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SelectInputImage()
{
  this->InputFileName = QFileDialog::getOpenFileName(this,
    tr("Select Input Image"), "", tr("Image Files (*.*)"));
  this->InputFile = QFileInfo(this->InputFileName);
  this->InputImageLabel->setText(this->InputFile.filePath());

  emit InputChanged(this->InputFileName.toStdString());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SelectBackboneFile()
{
  this->BackboneFileName = QFileDialog::getSaveFileName(this,
    tr("Select Backbone File"), "", tr("VTK files (*.vtk)"));
  this->BackboneFile = QFileInfo(this->BackboneFileName);
  //this->DataDir = this->BackboneFile.dir().absolutePath() + "/";
  this->BackboneOutputLabel->setText(this->BackboneFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SelectSkeletonFile()
{
  this->SkeletonFileName = QFileDialog::getSaveFileName(this,
    tr("Select Skeleton File"), "", tr("VTK files (*.vtk)"));
  this->SkeletonFile = QFileInfo(this->SkeletonFileName);
  this->SkeletonOutputLabel->setText(this->SkeletonFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetInputImage(const char *filename)
{
  this->InputFileName = filename;
  this->InputFile = QFileInfo(this->InputFileName);
  if(!this->InputFile.exists())
    {
    cerr << "ERROR: input file " << filename << " doesn't exist!" << endl;
    }
  this->InputImageLabel->setText(this->InputFile.filePath());
  emit InputChanged(this->InputFileName.toStdString());
  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetBackboneFile(const char *filename)
{
  this->BackboneFileName = filename;
  this->BackboneFile = QFileInfo(this->BackboneFileName);
  //this->DataDir = this->BackboneFile.dir().absolutePath() + "/";
  this->BackboneOutputLabel->setText(this->BackboneFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetSkeletonFile(const char *filename)
{
  this->SkeletonFileName = filename;
  this->SkeletonFile = QFileInfo(this->SkeletonFileName);
  this->SkeletonOutputLabel->setText(this->SkeletonFile.filePath());

  this->IntroPage->CheckIfComplete();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::RenderVolume(QFileInfo volumeFile)
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
void MDL2Wizard::RenderImage(vtkImageData *image, const char *imageName)
{
  //remove all old actors
  this->Renderer->RemoveAllViewProps();
  
  //change the window title so the user knows what file is being displayed
  this->RenderWidget->setWindowTitle(imageName);

  //create a volume from the specified image
  vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> mapper =
    vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
  mapper->SetInput(image);

  vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
  volume->SetMapper(mapper);

  double range[2];
  image->GetScalarRange(range);

  vtkSmartPointer<vtkVolumeProperty> volumeProperty =
    this->NewRGBVolumeProperty(range);
  volume->SetProperty(volumeProperty);

  this->Renderer->AddVolume(volume);
  this->Renderer->ResetCamera();
  this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::RenderPolyData(QFileInfo polyDataFile)
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
void MDL2Wizard::RenderPoints(QString pointsFileName)
{
  //remove all old actors
  this->Renderer->RemoveAllViewProps();

  //render the clean image so that the points have some context
  this->RenderImage(this->Helper->GetImageData(), "anisotropic diffusion");

  //change the window title so the user knows what file is being displayed
  QFileInfo info(pointsFileName);
  this->RenderWidget->setWindowTitle(info.fileName());

  vtkSmartPointer<vtkActor> actor =
    this->CreateActorFromPointsFile(pointsFileName);
  this->Renderer->AddActor(actor);

  this->Renderer->ResetCamera();
  this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkActor> MDL2Wizard::CreateActorFromPointsFile(
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
vtkSmartPointer<vtkVolume> MDL2Wizard::ConvertRawToVolume(const char *filename)
{
/*
  int sizeX = this->ImageSizeX.toInt();
  int sizeY = this->ImageSizeY.toInt();
  int sizeZ = this->ImageSizeZ.toInt();
*/
  int sizeX = 0;
  int sizeY = 0;
  int sizeZ = 0;

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
vtkSmartPointer<vtkActor> MDL2Wizard::CreateActorFromPolyDataFile(
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
vtkSmartPointer<vtkVolumeProperty> MDL2Wizard::NewRGBVolumeProperty(
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
void MDL2Wizard::MaskUsingGraphCuts()
{
  this->MaskUsingGraphCutsButton->setEnabled(false);
  this->OutputWindow->clear();

  emit StartMaskUsingGraphCuts();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::MaskSmallConnComp()
{
  this->MaskUsingGraphCutsButton->setEnabled(false);
  this->MaskSmallConnCompButton->setEnabled(false);
  this->OutputWindow->clear();

  int componentsSize = this->ComponentsSizeInput->text().toInt();
  emit StartMaskSmallConnComp(componentsSize);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::Integratedskel()
{
  this->IntegratedskelButton->setEnabled(false);
  this->OutputWindow->clear();
 
  double vectorMagnitude = this->VectorMagnitudeInput->text().toDouble();
  emit StartIntegratedskel(vectorMagnitude);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::CreateGraphAndMST1()
{
  this->IntegratedskelButton->setEnabled(false);
  this->CreateGraphAndMSTButton1->setEnabled(false);
  this->OutputWindow->clear();
 
  int edgeRange = this->EdgeRangeInput1->text().toInt();
  emit StartCreateGraphAndMST1(edgeRange);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::ErodeAndDilateNodeDegree1()
{
  this->IntegratedskelButton->setEnabled(false);
  this->CreateGraphAndMSTButton1->setEnabled(false);
  this->ErodeAndDilateButton1->setEnabled(false);
  this->OutputWindow->clear();
 
  int morphStrength = this->MorphStrengthInput1->text().toInt();
  emit StartErodeAndDilateNodeDegree1(morphStrength);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::BackboneExtract()
{
  this->IntegratedskelButton->setEnabled(false);
  this->CreateGraphAndMSTButton1->setEnabled(false);
  this->ErodeAndDilateButton1->setEnabled(false);
  this->BackboneExtractButton->setEnabled(false);
  this->OutputWindow->clear();

  emit StartBackboneExtract();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::BSplineFitting()
{
  this->IntegratedskelButton->setEnabled(false);
  this->CreateGraphAndMSTButton1->setEnabled(false);
  this->ErodeAndDilateButton1->setEnabled(false);
  this->BackboneExtractButton->setEnabled(false);
  this->BSplineFittingButton->setEnabled(false);
  this->OutputWindow->clear();

  unsigned int order = this->OrderInput->text().toUInt();
  unsigned int levels = this->LevelsInput->text().toUInt();
  emit StartBSplineFitting(order, levels);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::CreateGraphAndMST2()
{
  this->IntegratedskelButton->setEnabled(false);
  this->CreateGraphAndMSTButton2->setEnabled(false);
  this->OutputWindow->clear();
  
  int edgeRange = this->EdgeRangeInput2->text().toInt();
  emit StartCreateGraphAndMST2(edgeRange);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::ErodeAndDilateNodeDegree2()
{
  this->CreateGraphAndMSTButton2->setEnabled(false);
  this->ErodeAndDilateButton2->setEnabled(false);
  this->OutputWindow->clear();
  
  int morphStrength = this->MorphStrengthInput2->text().toInt();
  emit StartErodeAndDilateNodeDegree2(morphStrength);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayMaskUsingGraphCutsResults()
{
  this->MaskUsingGraphCutsButton->setEnabled(true);
  this->MaskSmallConnCompButton->setEnabled(true);
  this->MaskSmallConnCompButton->setFocus();
  //not sure what results to display here...
  if(!this->InteractiveExecution)
    {
    this->MaskSmallConnComp();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayMaskSmallConnCompResults()
{
  this->MaskUsingGraphCutsButton->setEnabled(true);
  this->MaskSmallConnCompButton->setEnabled(true);
  this->MaskSmallConnCompButton->setFocus();
  this->RenderImage(this->Helper->GetImageData(), "connected components");
  this->PreprocessingDone = true;
  this->PreprocessingPage->CheckIfComplete();
  if(!this->InteractiveExecution)
    {
    this->next();
    this->Integratedskel();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayIntegratedskelResults()
{
  this->IntegratedskelButton->setEnabled(true);
  this->CreateGraphAndMSTButton1->setEnabled(true);
  this->CreateGraphAndMSTButton1->setFocus();
  this->RenderPoints("out.seed");
  if(!this->InteractiveExecution)
    {
    this->CreateGraphAndMST1();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayCreateGraphAndMST1Results()
{
  this->IntegratedskelButton->setEnabled(true);
  this->CreateGraphAndMSTButton1->setEnabled(true);
  this->ErodeAndDilateButton1->setEnabled(true);
  this->ErodeAndDilateButton1->setFocus();
  this->RenderPolyData(QFileInfo("InitialMST.vtk"));
  if(!this->InteractiveExecution)
    {
    this->ErodeAndDilateNodeDegree1();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayErodeAndDilateNodeDegree1Results()
{
  this->IntegratedskelButton->setEnabled(true);
  this->CreateGraphAndMSTButton1->setEnabled(true);
  this->ErodeAndDilateButton1->setEnabled(true);
  this->BackboneExtractButton->setEnabled(true);
  this->BackboneExtractButton->setFocus();
  //not sure how to display degree.txt...
  if(!this->InteractiveExecution)
    {
    this->BackboneExtract();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayBackboneExtractResults()
{
  this->IntegratedskelButton->setEnabled(true);
  this->CreateGraphAndMSTButton1->setEnabled(true);
  this->ErodeAndDilateButton1->setEnabled(true);
  this->BackboneExtractButton->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->BSplineFittingButton->setFocus();
  this->RenderPolyData(QFileInfo("BackboneCandidate.vtk"));
  if(!this->InteractiveExecution)
    {
    this->BSplineFitting();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayBSplineFittingResults()
{
  this->IntegratedskelButton->setEnabled(true);
  this->CreateGraphAndMSTButton1->setEnabled(true);
  this->ErodeAndDilateButton1->setEnabled(true);
  this->BackboneExtractButton->setEnabled(true);
  this->BSplineFittingButton->setEnabled(true);
  this->RenderPolyData(QFileInfo("SmoothBackbone.vtk"));
  this->PhaseOneDone = true;
  this->PhaseOnePage->CheckIfComplete();
  if(!this->InteractiveExecution)
    {
    this->next();
    this->CreateGraphAndMST2();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayCreateGraphAndMST2Results()
{
  this->CreateGraphAndMSTButton2->setEnabled(true);
  this->ErodeAndDilateButton2->setEnabled(true);
  this->ErodeAndDilateButton2->setFocus();
  this->RenderPolyData(QFileInfo("InitialMST.vtk"));
  if(!this->InteractiveExecution)
    {
    this->ErodeAndDilateNodeDegree2();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::DisplayErodeAndDilateNodeDegree2Results()
{
  this->CreateGraphAndMSTButton2->setEnabled(true);
  this->ErodeAndDilateButton2->setEnabled(true);
  this->SaveOutputButton->setEnabled(true);
  this->SaveOutputButton->setFocus();
  //not sure how to display degree.txt...
  if(!this->InteractiveExecution)
    {
    this->SaveOutput();
    }
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SaveOutput()
{
  emit this->ReadyToSaveBackbone(
    this->BackboneFile.absoluteFilePath().toStdString().c_str());
  emit this->ReadyToSaveSkeleton(
    this->SkeletonFile.absoluteFilePath().toStdString().c_str());
  this->OutputWindow->append("Output saved!");
  this->SaveParameters();
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SaveParameters()
{
  QFile paramFile("MDLParameters.xml");
  if(!paramFile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
    return;
    }
  QTextStream out(&paramFile);
  out << "<?xml version='1.0'?>\n";
  out << "<parameters>\n";
  out << "<parameter name='ComponentsSize' value='";
  out << this->ComponentsSizeInput->text();
  out << "'></parameter>\n";
  out << "<parameter name='VectorMagnitude' value='";
  out << this->VectorMagnitudeInput->text();
  out << "'></parameter>\n";
  out << "<parameter name='EdgeRange1' value='";
  out << this->EdgeRangeInput1->text();
  out << "'></parameter>\n";
  out << "<parameter name='MorphStrength1' value='";
  out << this->MorphStrengthInput1->text();
  out << "'></parameter>\n";
  out << "<parameter name='Order' value='";
  out << this->OrderInput->text();
  out << "'></parameter>\n";
  out << "<parameter name='Levels' value='";
  out << this->LevelsInput->text();
  out << "'></parameter>\n";
  out << "<parameter name='EdgeRange2' value='";
  out << this->EdgeRangeInput2->text();
  out << "'></parameter>\n";
  out << "<parameter name='MorphStrength2' value='";
  out << this->MorphStrengthInput2->text();
  out << "'></parameter>\n";
  out << "</parameters>\n";
  paramFile.close();

  this->PhaseTwoDone = true;
  this->PhaseTwoPage->CheckIfComplete();
  if(!this->InteractiveExecution)
    {
    this->close();
    }
}


//-----------------------------------------------------------------------------
void MDL2Wizard::UpdateHelpWindow()
{
  if(this->HelpWindow->isVisible())
    {
    this->showHelp();
    }
}


/*
//-----------------------------------------------------------------------------
void MDL2Wizard::SetConnectedComponentsSize(const char *s)
{
  this->ConnectedComponentsSize = s;
  this->ComponentsSizeInput->setText(this->ConnectedComponentsSize);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetVectorMagnitude(const char *s)
{
  this->VectorMagnitude = s;
  this->VectorMagnitudeInput->setText(this->VectorMagnitude);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetEdgeRange(const char *s)
{
  this->EdgeRange = s;
  this->EdgeRangeInput1->setText(this->EdgeRange);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetMorphStrength(const char *s)
{
  this->MorphStrength = s;
  this->MorphStrengthInput1->setText(this->MorphStrength);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetBSplineOrder(const char *s)
{
  this->BSplineOrder = s;
  this->OrderInput->setText(this->BSplineOrder);
}

//-----------------------------------------------------------------------------
void MDL2Wizard::SetBSplineLevels(const char *s)
{
  this->BSplineLevels = s;
  this->LevelsInput->setText(this->BSplineLevels);
}
*/

//-----------------------------------------------------------------------------
void MDL2Wizard::showHelp()
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
void MDL2Wizard::closeEvent(QCloseEvent *event)
{
  event->accept();
  this->OutputWindow->close();
  this->HelpWindow->close();
}

