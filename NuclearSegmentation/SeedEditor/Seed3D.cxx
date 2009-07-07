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

/*
  Farsight ToolKit 3D Viewer:
  v1: implements render window functionality
    render window creation in "Seed3D::RenderWin"
    line objects mapped and actors created in "Seed3D::LineAct"
    actor added and window created in "Seed3D::addAct"
  v2: add picking functionality:
    access to observing a pick event is in "Seed3D::interact"
    the interacter calls for "Seed3D::PickPoint"
  v3: #include "Seed3D.h"
#include "Seed3DHelperClasses.h"include contourFilter and rayCast renderers
  v4: converted to Qt, member names changed to fit "VTK style" more closely.
*/
#include "Seed3D.h"
#include "Seed3DHelperClasses.h"

Seed3D::Seed3D(int argc, char **argv)
{

this->Initialize();
this->rayCast(argv[1],argv[2]);
this->AddVolumeSliders();
this->QVTK = 0;

}

Seed3D::~Seed3D()
{
  if(this->QVTK)
    {
    delete this->QVTK;
    }
  //delete this->tobj;
}

void Seed3D::Initialize()
{
  this->CreateGUIObjects();
  this->CreateLayout();
  this->CreateInteractorStyle();
  //this->CreateActors();
  this->resize(640, 480);
  this->setWindowTitle(tr("Seed Viewer"));
  this->QVTK->GetRenderWindow()->Render();
}

/*  set up the components of the interface */
void Seed3D::CreateGUIObjects()
{
  //Setup a QVTK Widget for embedding a VTK render window in Qt.
this->QVTK = new QVTKWidget(this);
this->Renderer = vtkRenderer::New();
this->QVTK->GetRenderWindow()->AddRenderer(this->Renderer);

}

void Seed3D::CreateLayout()
{
  //layout for the main window
  QGridLayout *viewerLayout = new QGridLayout(this);
  viewerLayout->addWidget(this->QVTK, 0, 0);
}




void Seed3D::CreateInteractorStyle()
{
  this->Interactor = this->QVTK->GetRenderWindow()->GetInteractor();
  //use trackball control for mouse commands
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkInteractorStyleTrackballCamera::New();
}


void Seed3D::rayCast(char *raySource, char *pointSource)
{
  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  ReaderType::Pointer i2spReader = ReaderType::New();
  i2spReader->SetFileName( raySource );
  try
  {
    i2spReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    //return EXIT_FAILURE;
    }
  std::cout << "Image Read " << std::endl;
//itk-vtk connector
  typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer connector= ConnectorType::New();
  connector->SetInput( i2spReader->GetOutput() );
  vtkImageToStructuredPoints *i2sp = vtkImageToStructuredPoints::New();
  i2sp->SetInput(connector->GetOutput());



  ImageType::SizeType size = i2spReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  vtkSmartPointer<vtkImageData> vtkim = vtkSmartPointer<vtkImageData>::New();
  vtkim->SetScalarTypeToUnsignedChar();
  vtkim->SetDimensions(size[0],size[1],size[2]);
  vtkim->SetNumberOfScalarComponents(1);
  vtkim->AllocateScalars();

  memcpy(vtkim->GetScalarPointer(),i2spReader->GetOutput()->GetBufferPointer(),size[0]*size[1]*size[2]*sizeof(unsigned char));

// Create transfer mapping scalar value to opacity
  vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();

  opacityTransferFunction->AddPoint(2,0.0);
  opacityTransferFunction->AddPoint(20,0.1);
  opacityTransferFunction->AddPoint(40,0.1);
  // Create transfer mapping scalar value to color
  // Play around with the values in the following lines to better vizualize data
  vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
  colorTransferFunction->AddRGBPoint(50.0,0.5,0.5,0);

  // The property describes how the data will look
  vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
  //  volumeProperty->ShadeOn();
    volumeProperty->SetInterpolationTypeToLinear();

  vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
  volumeMapper->SetSampleDistance(0.75);
  volumeMapper->SetInput(vtkim);

  // The volume holds the mapper and the property and
  // can be used to position/orient the volume
    vtkVolume *volume = vtkVolume::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
    volume->SetPickable(0);


 // Read the co-ordinates from a file:

point p;
std::vector<point> spPoint;
spPoint.clear();
ifstream inputstream(pointSource);

if (inputstream.fail())
{
cout << "ERROR: Incorrect File Type" << endl;
//return 1;
}
else
{
cout << "OK: File Imported" << endl;
}

inputstream >> p.x >> p.y>>p.z;
while(!inputstream.eof())
{
spPoint.push_back(p);
inputstream >> p.x >> p.y>>p.z;
}


// Create a sphere only at these particular coordinates only.

  // Create a float array which represents the points.
  vtkFloatArray* pcoords = vtkFloatArray::New();
   // Note that by default, an array has 1 component.
  // We have to change it to 3 for points

    pcoords->SetNumberOfComponents(3);
    pcoords->SetNumberOfTuples(spPoint.size());

  // Assign each tuple. There are 5 specialized versions of SetTuple:
  // SetTuple1 SetTuple2 SetTuple3 SetTuple4 SetTuple9
  // These take 1, 2, 3, 4 and 9 components respectively.

  for (int j=0; j<spPoint.size(); j++)
    {
    float pts[3] = {spPoint[j].x, spPoint[j].y , spPoint[j].z };
    pcoords->SetTuple(j, pts);
    }

  // Create vtkPoints and assign pcoords as the internal data array.
  vtkPoints* points = vtkPoints::New();
  points->SetData(pcoords);
// Create the dataset. In this case, we create a vtkPolyData
  vtkPolyData* polydata = vtkPolyData::New();
  // Assign points and cells
  polydata->SetPoints(points);


cout << "CS2?" << endl;

// create the spikes by glyphing the sphere with a cone.  Create the mapper
// and actor for the glyphs.
    vtkSphereSource *sphere = vtkSphereSource::New();
    vtkGlyph3D *glyph = vtkGlyph3D::New();
    glyph->SetInput(polydata);
    glyph->SetSource(sphere->GetOutput());
    glyph->SetVectorModeToUseNormal();
    glyph->SetScaleModeToScaleByVector();
    glyph->SetScaleFactor(10);
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInput(glyph->GetOutput());
    vtkLODActor *sphereActor = vtkLODActor::New();
    sphereActor->SetMapper(sphereMapper);
    Renderer->AddVolume(volume);
    Renderer->AddActor(sphereActor);

    vtkCamera *cam1 = Renderer->GetActiveCamera();
    Renderer->ResetCamera();


    this->Volume = volume;
    this->QVTK->GetRenderWindow()->Render();
    std::cout << "RayCast rendered \n";
}


void Seed3D::AddVolumeSliders()
{
  vtkSliderRepresentation2D *sliderRep = vtkSliderRepresentation2D::New();
  sliderRep->SetValue(0.8);
  sliderRep->SetTitleText("Opacity");
  sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep->GetPoint1Coordinate()->SetValue(0.2,0.1);
  sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep->GetPoint2Coordinate()->SetValue(0.8,0.1);
  sliderRep->SetSliderLength(0.02);
  sliderRep->SetSliderWidth(0.03);
  sliderRep->SetEndCapLength(0.01);
  sliderRep->SetEndCapWidth(0.03);
  sliderRep->SetTubeWidth(0.005);
  sliderRep->SetMinimumValue(0.0);
  sliderRep->SetMaximumValue(1.0);

  vtkSliderWidget *sliderWidget = vtkSliderWidget::New();
  sliderWidget->SetInteractor(Interactor);
  sliderWidget->SetRepresentation(sliderRep);
  sliderWidget->SetAnimationModeToAnimate();

  vtkSlider2DCallbackBrightness *callback_brightness = vtkSlider2DCallbackBrightness::New();
  callback_brightness->volume = this->Volume;
  sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback_brightness);
  sliderWidget->EnabledOn();


// slider 2

  vtkSliderRepresentation2D *sliderRep2 = vtkSliderRepresentation2D::New();
  sliderRep2->SetValue(0.8);
  sliderRep2->SetTitleText("Brightness");
  sliderRep2->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint1Coordinate()->SetValue(0.2,0.9);
  sliderRep2->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint2Coordinate()->SetValue(0.8,0.9);
  sliderRep2->SetSliderLength(0.02);
  sliderRep2->SetSliderWidth(0.03);
  sliderRep2->SetEndCapLength(0.01);
  sliderRep2->SetEndCapWidth(0.03);
  sliderRep2->SetTubeWidth(0.005);
  sliderRep2->SetMinimumValue(0.0);
  sliderRep2->SetMaximumValue(1.0);

  vtkSliderWidget *sliderWidget2 = vtkSliderWidget::New();
  sliderWidget2->SetInteractor(Interactor);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContrast *callback_contrast = vtkSlider2DCallbackContrast::New();
  callback_contrast->volume = this->Volume;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contrast);
  sliderWidget2->EnabledOn();

}













/*void Seed3D::closeEvent(QCloseEvent *event)
{
  event->accept();
}  this->Renderer = vtkRenderer::New();*/
