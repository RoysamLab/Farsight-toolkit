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
this->mode = 0;
this->counter =0;
this->stateAdd = 0;
this->stateDelete = 0;
this->stateSplit = 0;
this->stateMerge = 0;
this->flag =1;
this->QVTK = 0;
this->Initialize();
this->Renderer->AddActor(ConeActor);
this->rayCast(argv[1],argv[2]);
this->AddVolumeSliders();
}

Seed3D::~Seed3D()
{
  if(this->QVTK)
    {
    delete this->QVTK;
    }
}

void Seed3D::Initialize()
{
  this->CreateGUIObjects();
  this->CreateLayout();
  this->CreateInteractorStyle();
  this->resize(800, 600);
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
  this->QVTK1 = new QVTKWidget(this);
  this->Renderer1 = vtkRenderer::New();
  this->QVTK1->GetRenderWindow()->AddRenderer(this->Renderer1);
  this->QVTK2 = new QVTKWidget(this);
  this->Renderer2 = vtkRenderer::New();
  this->QVTK2->GetRenderWindow()->AddRenderer(this->Renderer2);

  //Setup the buttons that the user will use to interact with this program. 
  this->AddBox = new QCheckBox("&Add Seeds", this);
  //this->ClearButton = new QPushButton("&clear selection", this);
  this->DeleteBox = new QCheckBox("&Delete Seeds", this);
  this->MergeBox = new QCheckBox("&Merge Nuclei", this);
  this->SplitBox = new QCheckBox("&Split Nuclei", this);
  this->UndoDelBox = new QCheckBox("&Undo Selections for Delete/Merge", this);	
  this->PlaceButton = new QPushButton("&Place Seed", this);
  this->ApplyButton = new QPushButton("&Apply", this);std::cout << "RayCast rendered \n";  
  connect(this->AddBox, SIGNAL(stateChanged(int)), this, SLOT(AddSeed()));
  connect(this->DeleteBox, SIGNAL(stateChanged(int)), this, SLOT(DeleteSeed()));
  connect(this->SplitBox, SIGNAL(stateChanged(int)), this, SLOT(SplitSeeds()));
  connect(this->MergeBox, SIGNAL(stateChanged(int)), this, SLOT(MergeSeeds()));
  connect(this->UndoDelBox, SIGNAL(stateChanged(int)), this, SLOT(UndoDeleteSeeds()));
  connect(this->PlaceButton, SIGNAL(clicked()), this, SLOT(PlaceSeed()));
  connect(this->ApplyButton, SIGNAL(clicked()), this, SLOT(Apply()));
}

void Seed3D::CreateLayout()
{
 //layout for the main window
  QGridLayout *buttonLayout = new QGridLayout();
  buttonLayout->addWidget(this->AddBox, 0, 2);
  buttonLayout->addWidget(this->DeleteBox, 0, 3);
  buttonLayout->addWidget(this->SplitBox, 0, 4);
  buttonLayout->addWidget(this->MergeBox, 0, 5);
  buttonLayout->addWidget(this->UndoDelBox, 0, 6);
  buttonLayout->addWidget(this->PlaceButton, 1, 3);
  buttonLayout->addWidget(this->ApplyButton, 1, 4);
  QGridLayout *viewerLayout = new QGridLayout(this);
  viewerLayout->addWidget(this->QVTK, 0, 0,1,2);
  viewerLayout->addWidget(this->QVTK1, 1, 0);
  viewerLayout->addWidget(this->QVTK2, 1, 1); 	
  viewerLayout->addLayout(buttonLayout, 2, 0); 	
}




void Seed3D::CreateInteractorStyle()
{
  this->Interactor = this->QVTK->GetRenderWindow()->GetInteractor();
  this->Interactor1 = this->QVTK1->GetRenderWindow()->GetInteractor();
  this->Interactor2 = this->QVTK2->GetRenderWindow()->GetInteractor();
  //use trackball control for mouse commands
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkInteractorStyleTrackballCamera::New();

 vtkSmartPointer<vtkInteractorStyleRubberBand2D> style1 =
    vtkInteractorStyleRubberBand2D::New();
 
  this->Interactor->SetInteractorStyle(style);
  this->PointPicker = vtkPointPicker::New();
  this->PointPicker->SetTolerance(0.004);
  this->JustPicker = vtkPropPicker::New();
  this->Interactor->SetPicker(this->JustPicker);
  this->Interactor->SetPicker(this->PointPicker);
  this->isPicked = vtkCallbackCommand::New();
  this->isPicked->SetCallback(PickCell);

  //isPicked caller allows observer to intepret click 
  this->isPicked->SetClientData(this);            
  this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,isPicked); 
  this->sphereWidget = vtkSphereWidget::New();;
  this->sphereWidget->PlaceWidget();
  this->sphereWidget->SetInteractor(this->Interactor);
  this->Interactor1->SetInteractorStyle(style1);
  this->Interactor2->SetInteractorStyle(style1);
  //this->Interactor2 = vtkRenderWindowInteractor::New();
  //this->ImageViewer->SetupInteractor(Interactor2);

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

   this->VTKim = vtkim;

   memcpy(vtkim->GetScalarPointer(),i2spReader->GetOutput()->GetBufferPointer(),size[0]*size[1]*size[2]*sizeof(unsigned char));
   vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
   vtkInteractorStyleTrackballCamera::New();

// Create transfer AddObservermapping scalar value to opacity
   vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();
   opacityTransferFunction->AddPoint(2,0.0);
  // opacityTransferFunction->AddPoint(20,0.1);
   opacityTransferFunction->AddPoint(100,0.5);

  // Create transfer mapping scalar value to color
  // Play around with the values in the following lines to better vizualize data
   vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
   colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
   colorTransferFunction->AddRGBPoint(50.0,0.5,0.5,0.5);

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
    this->Volume = volume;
    

 // Create a sphere only at these particular coordinates only.

    std::vector<point> spPoint1;
    std::vector<point> MarkedPoints;
    std::vector<point> MarkedPoints2add;
	std::vector<point> dup_points; //used while deleteing seeds
    
	spPoint1 = ReadPoints(pointSource);
    this->p = spPoint1;
    this->dup_points = spPoint1;
		
    // Create a float array which represents the points.
    this->pcoords = vtkFloatArray::New();
    this->delpcoords = vtkFloatArray::New();	
     // Note that by default, an array has 1 component.
    // We have to change it to 3 for points
    this->pcoords->SetNumberOfComponents(3);
    this->pcoords->SetNumberOfTuples(spPoint1.size());
    this->delpcoords->SetNumberOfComponents(3);
    this->delpcoords->SetNumberOfTuples(1);		

  // Assign each tuple. There are 5 specialized versions of SetTuple:
  // SetTuple1 SetTuple2 SetTuple3 SetTuple4 SetTuple9
  // These take 1, 2, 3, 4 and 9 components respectively.

  	
  for (int j=0; j<spPoint1.size(); j++)
    {
    float pts[3] = {spPoint1[j].x,spPoint1[j].y , spPoint1[j].z };	
    this->pcoords->SetTuple(j, pts);
    }

    float pts[3] = {0.0,0.0,0.0};
    this->delpcoords->SetTuple(0,pts);	
     
// Create vtkPoints and assign pcoords as the internal data array.
  this->point1 = vtkPoints::New();
  this->point2 = vtkPoints::New(); //Seeds marked for deletion
  this->point3 = vtkPoints::New(); //Seeds marked for addition
  this->allpoints = vtkPoints::New();//Yousefseg downstream
  
  this->point1->SetData(this->pcoords);	  
  this->point2->SetData(this->delpcoords);
  this->point3->SetData(this->delpcoords);
 
  // Assign points and cells
  this->polydata1 = vtkPolyData::New();
  this->polydata1->SetPoints(this->point1);
  this->polydata2 = vtkPolyData::New();	
  this->polydata2->SetPoints(this->point2);
  
  this->polydata3 = vtkPolyData::New();	
  this->polydata3->SetPoints(this->point3);

   
// create the spikes by glyphing the sphere with a cone.  Create the mapper
// and actor for the glyphs.

    vtkSphereSource *sphere = vtkSphereSource::New();
    vtkGlyph3D *glyph = vtkGlyph3D::New();
    vtkGlyph3D *delglyph = vtkGlyph3D::New();
    vtkGlyph3D *addglyph = vtkGlyph3D::New();
    glyph->SetInput(this->polydata1);
    delglyph->SetInput(this->polydata2);
    addglyph->SetInput(this->polydata3);
    glyph->SetSource(sphere->GetOutput());
    delglyph->SetSource(sphere->GetOutput());
    addglyph->SetSource(sphere->GetOutput());
    glyph->SetVectorModeToUseNormal();
    delglyph->SetVectorModeToUseNormal();
    addglyph->SetVectorModeToUseNormal();
    glyph->SetScaleModeToScaleByVector();
    delglyph->SetScaleModeToScaleByVector();
    addglyph->SetScaleModeToScaleByVector();
    glyph->SetScaleFactor(1);
    delglyph->SetScaleFactor(1);
    addglyph->SetScaleFactor(1);
    glyph->GeneratePointIdsOn();
    delglyph->GeneratePointIdsOn();
    addglyph->GeneratePointIdsOn();
    this->Glyph = glyph;
    this->delglyph = delglyph;
    this->addglyph = addglyph;
    this->imActor = imActor;
    
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper *delsphereMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper *addsphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInput(this->Glyph->GetOutput());
    delsphereMapper->SetInput(this->delglyph->GetOutput());
    addsphereMapper->SetInput(this->addglyph->GetOutput());
    vtkLODActor *sphereActor = vtkLODActor::New();
    vtkLODActor *delsphereActor = vtkLODActor::New();
    vtkLODActor *addsphereActor = vtkLODActor::New();
    sphereActor->SetMapper(sphereMapper);
    delsphereActor->SetMapper(delsphereMapper);
    delsphereActor->GetProperty()->SetColor(0.5,0.0,0.0); 	
    addsphereActor->SetMapper(addsphereMapper);
    addsphereActor->GetProperty()->SetColor(0.0,0.5,0.0);
    Renderer->AddVolume(volume);
    Renderer->AddActor(sphereActor);
    Renderer->AddActor(delsphereActor);
    Renderer->AddActor(addsphereActor);	
    this->SphereMapper = sphereMapper;
    this->SphereActor = sphereActor;
    this->DelSphereMapper = delsphereMapper;
    this->DelSphereActor = delsphereActor;		
    this->AddSphereMapper = addsphereMapper;
    this->AddSphereActor = addsphereActor;

    vtkCamera *cam1 = Renderer->GetActiveCamera();
    cam1->SetViewUp (0, 1, 0);
    cam1->SetPosition (0, 0, 50);
    Renderer->ResetCamera();	
    this->Volume = volume;
    std::cout << "RayCast rendered \n";  
 
    
    //Remove the initial red and green glyphs. 
    vtkDataArray* delpoint = this->point2->GetData();    
    delpoint->RemoveTuple((vtkIdType)0);
    this->point2->SetData(delpoint);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);
    
    vtkDataArray* addpoint = this->point3->GetData();    
    addpoint->RemoveTuple((vtkIdType)0);
    this->point3->SetData(addpoint);
    this->addglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	


	
  // Create the widget and its representation
  vtkPointHandleRepresentation3D *handle = vtkPointHandleRepresentation3D::New();
  vtkHandleWidget *widget = vtkHandleWidget::New();
  widget->SetInteractor(this->Interactor);
  widget->SetRepresentation(handle);
  this->widget = widget;
  this->handle = handle; 
  //this->rep = rep;
 
 
  vtkSeedCallback *scbk = vtkSeedCallback::New();
  widget->AddObserver(vtkCommand::EndInteractionEvent,scbk);
  widget->SetEnabled(true);

 static double axialElements[16] = {
           1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1, 0,
           0, 0, 0, 1 };
	
  vtkMatrix4x4 *resliceAxes = vtkMatrix4x4::New();
  resliceAxes->DeepCopy(axialElements);
  vtkImageReslice *reslice = vtkImageReslice::New();
  reslice->SetInput(vtkim);
  reslice->SetOutputDimensionality(2);
  reslice->SetResliceAxes(resliceAxes);
  reslice->SetInterpolationModeToLinear();
  
  double* origincalc = this->Volume->GetBounds();
  reslice->SetResliceAxesOrigin((origincalc[1]-origincalc[0])/2.0,(origincalc[3]-origincalc[2])/2.0,(origincalc[5]-origincalc[4])/2.0);
  this->Reslice = reslice;


  vtkDataSetMapper *im_mapper = vtkDataSetMapper::New();
  im_mapper->SetInput(reslice->GetOutput());
  vtkActor *imActor = vtkActor::New();
  imActor->SetMapper(im_mapper);
  
  this->Renderer1->AddActor(imActor);
  this->QVTK1->GetRenderWindow()->Render();
  this->Renderer1->ResetCamera();	
  this->imActor = imActor;
	

  static double sagittalElements[16] = { 0, 0,-1, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1 }; 
  vtkMatrix4x4 *resliceAxes1 = vtkMatrix4x4::New();
  resliceAxes->DeepCopy(sagittalElements);

  vtkImageReslice *reslice1 = vtkImageReslice::New();
  reslice1->SetInput(vtkim);
  reslice1->SetOutputDimensionality(2);
  reslice1->SetResliceAxes(resliceAxes1);
  reslice1->SetInterpolationModeToLinear();
  
  origincalc = this->Volume->GetBounds();
  reslice1->SetResliceAxesOrigin((origincalc[1]-origincalc[0])/2.0,(origincalc[3]-origincalc[2])/2.0,(origincalc[5]-origincalc[4])/2.0);
  this->Reslice1 = reslice1;

  vtkDataSetMapper *im_mapper1 = vtkDataSetMapper::New();
  im_mapper1->SetInput(reslice1->GetOutput());
  vtkActor *imActor1 = vtkActor::New();
  imActor1->SetMapper(im_mapper1);
  
  this->Renderer2->AddActor(imActor1);
  this->QVTK2->GetRenderWindow()->Render();
  this->Renderer2->ResetCamera();
  this->imActor1 = imActor1;

  
  scbk->handle = this->handle;
  scbk->handle1 = this->handle1;
  scbk->handle2 = this->handle2;
  scbk->reslice = this->Reslice;
  scbk->reslice1 = this->Reslice1;
  scbk->imActor1 = this->imActor1; 
  scbk->imActor =  this->imActor;
  scbk->QVTK1 = this->QVTK1;
  scbk->QVTK2 = this->QVTK2;
  scbk->vol = this->Volume;  

  vtkPointHandleRepresentation2D *handle1 = vtkPointHandleRepresentation2D::New();
  vtkHandleWidget *widget1 = vtkHandleWidget::New();
  widget1->SetInteractor(this->Interactor1);
  widget1->SetRepresentation(handle1);
  this->widget1 = widget1;
  this->handle1 = handle1; 
	  
  vtkSeedCallback1 *scbk1 = vtkSeedCallback1::New();
  widget1->AddObserver(vtkCommand::EndInteractionEvent,scbk1);
  widget1->SetEnabled(true);
  

  vtkPointHandleRepresentation2D *handle2 = vtkPointHandleRepresentation2D::New();
  vtkHandleWidget *widget2 = vtkHandleWidget::New();
  widget2->SetInteractor(this->Interactor2);
  widget2->SetRepresentation(handle2);
  this->widget2 = widget2;
  this->handle2 = handle2; 
	  
  vtkSeedCallback2 *scbk2 = vtkSeedCallback2::New();
  widget2->AddObserver(vtkCommand::EndInteractionEvent,scbk2);
  widget2->SetEnabled(true);
  scbk->handle1 = this->handle1;
  scbk->handle2 = this->handle2;
  scbk1->handle1= this->handle1;	
  scbk1->handle2= this->handle2;	
  scbk1->QVTK2 = this->QVTK2;
  scbk1->QVTK = this->QVTK;	
  scbk1->vol = this->Volume;
  scbk1->handle = this->handle; 
  scbk2->handle1= this->handle1;	
  scbk2->handle2= this->handle2;
  scbk2->QVTK1 = this->QVTK1;	
  scbk2->vol = this->Volume;
  scbk2->QVTK = this->QVTK;
  scbk2->handle = this->handle; 
  std::cout << "RayCast rendered \n";  
  this->QVTK->GetRenderWindow()->Render();    

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



//seed size

  vtkSliderRepresentation2D *sliderRep3 = vtkSliderRepresentation2D::New();
  sliderRep3->SetValue(0.8);
  sliderRep3->SetTitleText("Seed Size");
  sliderRep3->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep3->GetPoint1Coordinate()->SetValue(0.1,0.2);
  sliderRep3->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep3->GetPoint2Coordinate()->SetValue(0.1,0.9);
  sliderRep3->SetSliderLength(0.01);
  sliderRep3->SetSliderWidth(0.03);
  sliderRep3->SetEndCapLength(0.01);
  sliderRep3->SetEndCapWidth(0.03);
  sliderRep3->SetTubeWidth(0.005);
  sliderRep3->SetMinimumValue(1.0);
  sliderRep3->SetMaximumValue(15.0);

  vtkSliderWidget *sliderWidget3 = vtkSliderWidget::New();
  sliderWidget3->SetInteractor(Interactor);
  sliderWidget3->SetRepresentation(sliderRep3);
  sliderWidget3->SetAnimationModeToAnimate();
  //sliderWidget3->EnabledOn();
  vtkSlider2DCallbackSeedSize *callback_seedsize = vtkSlider2DCallbackSeedSize::New();
  callback_seedsize->Glyph = this->Glyph;
  callback_seedsize->addglyph = this->addglyph;
  callback_seedsize->delglyph = this->delglyph;	
  sliderWidget3->AddObserver(vtkCommand::InteractionEvent,callback_seedsize);
  sliderWidget3->EnabledOn();
}



void   Seed3D::PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{ 
   Seed3D* seed = (Seed3D*)clientdata;
 /*  PickPoint allows fot the point id and coordinates to be returned 
  as well as adding a marker on the last picked point
  R_click to select point on line  */

    int *pos = seed->Interactor->GetEventPosition();
    seed->Interactor->GetPicker()->Pick(pos[0],pos[1],0.0,seed->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
    vtkPointPicker *point_picker = (vtkPointPicker *)seed->Interactor->GetPicker();
    double pickPos[3];
    seed->PointPicker->GetPickPosition(pickPos);    //this is the coordinates of the pick  
    cout<<"Point Selected " << pickPos[0]<<"-"<<pickPos[1]<<"-"<<pickPos[2]<<endl;

    //Useful to check if clicked on a seed !   
//If I am in Delete mode 

if((seed->mode == 2||seed->mode == 4) && seed->flag == 1){

    vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
    int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 
    cout<<pID<<" is the iD value"<<endl;		

if(pID<=seed->dup_points.size())    //The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    float dist1;
    int index;
    float finpt[3];
    for (int j=0; j<seed->dup_points.size(); j++)
    {
	
    	float p1[3] = {seed->dup_points[j].x, seed->dup_points[j].y ,seed->dup_points[j].z };	
    	float dist1= sqrt(pow((p1[0]-pickPos[0]),2) + pow((p1[1]-pickPos[1]),2) + pow((p1[2]-pickPos[2]),2));   
        if (dist1<dist)
       		{
		dist = dist1;
        	finpt[0] = p1[0];
		finpt[1] = p1[1];
		finpt[2] = p1[2];
                index = j;
       		}    

}
	cout<<index<<" is the inD value"<<endl;
	seed->dup_points.erase(seed->dup_points.begin()+index);

       //Remove the glyph		
        vtkDataArray* points2del = seed->point1->GetData();    
        vtkDataArray* points2delred;
        points2del->RemoveTuple((vtkIdType)index);
        seed->point1->SetData(points2del);
	seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);

      // Add the new red glyph 
        seed->point2->InsertNextPoint(finpt);
	    seed->polydata2->SetPoints(seed->point2);
    	seed->delglyph->SetInput(seed->polydata2);
        seed->DelSphereMapper->SetInput(seed->delglyph->GetOutput());	      
        seed->delglyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        seed->QVTK->GetRenderWindow()->Render();
       

        

//Keep Track of seeds marked for deletion/merge in a vector 
//Useful while undoing it 

vtkIdType Id;
double* pointz;
point p;
Id = seed->point2->GetNumberOfPoints();

    pointz = seed->point2->GetPoint(Id-1);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    seed->MarkedPoints.push_back(p);

}
}


//If all points are marked, don't allow to delete anymore seeds.
//Set the flag so that it does not enter into delete mode functionality 

if(seed->MarkedPoints.size()==seed->p.size() || seed->dup_points.size()==0) 
{
seed->flag =0; 	 
}
else{
seed->flag =1; 	 
   }


if(seed->mode == 5){

    vtkDataArray* pointIds = seed->delglyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
    int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 
    cout<<pID<<" is the iD value"<<endl;		

if(pID<=seed->MarkedPoints.size())    //The ids of non-seed points is much greater than the ids of the seed points 

{    
   float dist =1000.00;
   float dist1;
   int index;
   float finpt[3];
   for (int j=0; j<seed->MarkedPoints.size(); j++)
    {
    	float p1[3] = {seed->MarkedPoints[j].x, seed->MarkedPoints[j].y ,seed->MarkedPoints[j].z };	
    	float dist1= sqrt(pow((p1[0]-pickPos[0]),2) + pow((p1[1]-pickPos[1]),2) + pow((p1[2]-pickPos[2]),2));   
        if (dist1<dist)
       		{
		dist = dist1;
        finpt[0] = p1[0];
		finpt[1] = p1[1];
		finpt[2] = p1[2];
        index = j;
       		}
    }
	seed->MarkedPoints.erase(seed->MarkedPoints.begin()+index);
     
     //Remove the red glyph		
        vtkDataArray* points2put = seed->point2->GetData();    
        points2put->RemoveTuple((vtkIdType)index);
        seed->point2->SetData(points2put);
	seed->delglyph->SetScaleFactor(seed->delglyph->GetScaleFactor()+0.0001);
	cout<<"MarkedPoints " <<seed->MarkedPoints.size()<<endl;
        cout<<"Number of Points in point2 "<<seed->point2->GetNumberOfPoints()<<endl;
      
     // Add the new silver glyph 
        seed->point1->InsertNextPoint(finpt);
	seed->polydata1->SetPoints(seed->point1);
    	seed->Glyph->SetInput(seed->polydata1);
        seed->SphereMapper->SetInput(seed->Glyph->GetOutput());	      
        seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        seed->QVTK->GetRenderWindow()->Render();
	
//Update the duplicate Seed list vector

point addPoint;
addPoint.x = finpt[0];
addPoint.y = finpt[1];
addPoint.z = finpt[2];
seed->dup_points.push_back(addPoint);
cout<<"Total seeds "<<seed->dup_points.size()<<endl;
		}
        
}
/*
//Set the state for undoDelbox depending on whether there are marked seeds
if(seed->point2->GetNumberOfPoints()==0)
{
seed->UndoDelBox->setCheckState((Qt::CheckState)0);
seed->UndoDelBox->setCheckable(false);
}
else
{
seed->UndoDelBox->setCheckable(true);
}
*/
seed->Check();	
}



std::vector<point> Seed3D::ReadPoints(char* pointzource)

{
 // Read the co-ordinates from a file: 
point p;
std::vector<point> spPoint;
spPoint.clear();
ifstream inputstream(pointzource);

if (inputstream.fail())
{
cout << "ERROR: Incorrect File Type" << endl;
//return 1;
}
else
{
cout << "OK: File Imported" << endl;
}
// Create a sphere only at these particular coordinates only.

  // Create a float array which represents the points.
  this->pcoords = vtkFloatArray::New();
  // Note that by default, an array has 1 component.
  // We have to change it to 3 for points

    this->pcoords->SetNumberOfComponents(3);
    this->pcoords->SetNumberOfTuples(spPoint.size());	
    inputstream >> p.x >> p.y>>p.z;
    while(!inputstream.eof())
  {
    spPoint.push_back(p);
    inputstream >> p.x >> p.y>>p.z;
  }
    return(spPoint);
}




void Seed3D::PlaceSeed()
{

	if(this->mode==1 ||this->mode==3){	
        std::cout << "Place Seed Activated" << std::endl;	
        double* p1 = this->handle1->GetWorldPosition();
        double* p2 = this->handle2->GetWorldPosition();       
	    double* bounds = this->Volume->GetBounds();
	    float placePoint[3] = {(bounds[1]/2.0)+p2[0],(bounds[3]/2.0)+p2[1],(bounds[5]/2.0)+p1[1]};	
        //this->point1->InsertNextPoint(placePoint);
        this->point3->InsertNextPoint(placePoint);
        //this->allpoints->InsertNextPoint(placePoint);//downstream to yousef core
        //this->polydata1->SetPoints(point1);
    	this->polydata3->SetPoints(this->point3);
        //this->Glyph->SetInput(this->polydata1);
        //this->Glyph->SetInput(this->polydata1);
        this->addglyph->SetInput(this->polydata3);
        this->addglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        //this->SphereMapper->SetInput(this->Glyph->GetOutput());	
        this->AddSphereMapper->SetInput(this->addglyph->GetOutput());
        this->Glyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        this->QVTK->GetRenderWindow()->Render();



		
//Keep Track of seeds marked for addition/split in a vector 
//Useful while undoing it 


vtkIdType Id;
double* pointz;
point p;
Id = this->point3->GetNumberOfPoints();

    pointz = this->point3->GetPoint(Id-1);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    this->MarkedPoints2add.push_back(p);



	if(this->stateDelete){ this->mode = 2;} 
	if(this->stateSplit) { this->mode = 3;}
	if(this->stateMerge) { this->mode = 4;}
	if(this->stateAdd)   { this->mode = 1;}

                      }
        
}



void Seed3D::AddSeed()
{
this->mode = 1;
this->stateAdd = this->AddBox->checkState();
if(this->MarkedPoints.size()!=0 && this->stateAdd){
this->msgBox = new QMessageBox(this);
//if(this->stateDelete)
this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes?");
if(this->stateMerge)
this->msgBox->setText("You have marked seeds for merging. Do you want to apply the changes?"); 
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
  // if(this->stateDelete)
   this->mode =2;
   if(this->stateMerge)
   this->mode =4;
   Apply();
   if(this->stateDelete )
   this->DeleteBox->setCheckState((Qt::CheckState)0); 
   if(this->stateMerge)
   this->MergeBox->setCheckState((Qt::CheckState)0);
   this->mode =1;
   this->AddBox->setCheckState((Qt::CheckState)2);
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //DeleteSeed(); 
	   if(this->stateDelete){
	   this->AddBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->DeleteBox->setCheckState((Qt::CheckState)2);
	   this->mode =2;}



	   if(this->stateMerge){
	   this->AddBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);
           this->MergeBox->setCheckState((Qt::CheckState)2);
	              this->mode =4;}


          if(this->stateUndoDel){
	    this->AddBox->setCheckState((Qt::CheckState)0);
	    this->UndoDelBox->setCheckState((Qt::CheckState)2);	
	    this->mode =5;
	   }




       break;
                            }
   default:
       // should never be reached
       break;
 }

}


else if(this->MarkedPoints2add.size()!=0 && this->stateAdd){
	if(stateSplit){
this->msgBox = new QMessageBox(this);
this->msgBox->setText("You have marked seeds for splitting nuclei. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   this->mode =3;
   Apply();
   this->SplitBox->setCheckState((Qt::CheckState)0); 
   this->mode =1;
   this->AddBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //SplitSeeds(); 
       this->AddBox->setCheckState((Qt::CheckState)0);
       this->UndoDelBox->setCheckState((Qt::CheckState)0);
       this->SplitBox->setCheckState((Qt::CheckState)2);
       this->mode =3;
       break;
                            }
   default:
       // should never be reached
       break;
 }

}
}

else{

this->AddBox->setCheckState((Qt::CheckState)this->stateAdd);
if(this->stateAdd){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateSplit) { this->SplitBox->setCheckState((Qt::CheckState)0); }
	if(this->stateMerge) { this->MergeBox->setCheckState((Qt::CheckState)0); }
  	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
        this->mode = 1;
  	cout<<"Add"<<endl; 
		 }
else{
this->mode=0;
}

}
}

void Seed3D::DeleteSeed()
{
if(this->stateUndoDel) {this->UndoDelBox->setCheckState((Qt::CheckState)0);}
this->mode = 2;
this->stateDelete = this->DeleteBox->checkState();
if(this->MarkedPoints2add.size()!=0 && this->stateDelete){
this->msgBox = new QMessageBox(this);
if(this->stateAdd)
this->msgBox->setText("You have marked seeds for addition. Do you want to apply the changes?");
if(this->stateSplit)
this->msgBox->setText("You have marked seeds for splitting. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   if(this->stateAdd)
   this->mode =1;
   if(this->stateSplit)
   this->mode =3;
   Apply();
   if(this->stateAdd)
   this->AddBox->setCheckState((Qt::CheckState)0);
   if(this->stateSplit)
   this->SplitBox->setCheckState((Qt::CheckState)0); 
   this->mode =2;
   this->DeleteBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //AddSeed(); 
       
	   if(this->stateAdd){
	   this->DeleteBox->setCheckState((Qt::CheckState)0);
           this->UndoDelBox->setCheckState((Qt::CheckState)0);	
	   this->AddBox->setCheckState((Qt::CheckState)2);
	   this->mode =1;}
	   if(this->stateSplit){
	   this->DeleteBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->SplitBox->setCheckState((Qt::CheckState)2);
	   this->mode =3;}
	   
       break;
                            }
   default:
       // should never be reached
       break;
 }

}

else if(this->MarkedPoints.size()!=0 && this->stateDelete){
	if(stateMerge){
this->msgBox = new QMessageBox(this);
this->msgBox->setText("You have marked seeds for merging. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   this->mode =4;
   Apply();
   this->MergeBox->setCheckState((Qt::CheckState)0); 
   this->mode =2;
   this->DeleteBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //AddSeed(); 
	   this->DeleteBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->MergeBox->setCheckState((Qt::CheckState)2);
           this->mode =4;
	   break;
                            }
  default:
       // should never be reached
       break;
 }

}
}
else{
this->DeleteBox->setCheckState((Qt::CheckState)this->stateDelete);
if(this->stateDelete){
if(this->stateSplit) { this->SplitBox->setCheckState((Qt::CheckState)0); }
if(this->stateMerge) { this->MergeBox->setCheckState((Qt::CheckState)0); }
if(this->stateAdd) { cout<<"hi"<<endl;this->AddBox->setCheckState((Qt::CheckState)0); }
if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
this->mode = 2;
cout<<"Delete"<<endl; 
}

else{
this->mode=0;
}
}
}

void Seed3D::SplitSeeds()
{

this->mode =3;
this->stateSplit = this->SplitBox->checkState();
if(MarkedPoints.size()!=0 && this->stateSplit){

	this->msgBox = new QMessageBox(this);
	this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes?");
	//if(this->stateDelete)	
	//this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes?");
	if(this->stateMerge)	
	this->msgBox->setText("You have seeds marked for merging. Do you want to apply the changes?");
	this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
	this->msgBox->setDefaultButton(QMessageBox::Yes);
	int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   
   //if(this->stateDelete)
   this->mode =2;
   if(this->stateMerge)
   this->mode =4;
   Apply();
   if(this->stateDelete)	
   this->DeleteBox->setCheckState((Qt::CheckState)0); 
   if(this->stateMerge)	
   this->MergeBox->setCheckState((Qt::CheckState)0); 
   this->mode =3;
   this->SplitBox->setCheckState((Qt::CheckState)2);
   break;
	}

   case QMessageBox::No:{
       // Cancel was clicked
       //DeleteSeed(); 
	   if(this->stateDelete){
           this->SplitBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->DeleteBox->setCheckState((Qt::CheckState)2);
	   this->mode =2;
}
       if(this->stateMerge){
	   this->SplitBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->MergeBox->setCheckState((Qt::CheckState)2);
	   this->mode =4;
	   }
	  

if(this->stateMerge){
	   this->SplitBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
           this->MergeBox->setCheckState((Qt::CheckState)2);
	   this->mode =4;
	   }
	

if(this->stateUndoDel){
	   this->SplitBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)2);	
	   this->mode =5;
	   }
       break;

                            }
   default:
       // should never be reached
       break;
 }

}

else if(this->MarkedPoints2add.size()!=0 && this->stateSplit){
	if(stateAdd){
this->msgBox = new QMessageBox(this);
this->msgBox->setText("You have seeds marked for addition. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   this->mode =1;
   Apply();
   this->AddBox->setCheckState((Qt::CheckState)0);
   this->mode =3;
   this->SplitBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   
   
   // Split or add howdyu decide ??
   
   case QMessageBox::No:{
       // Cancel was clicked
       //AddSeed(); 
       
       this->SplitBox->setCheckState((Qt::CheckState)0);
       this->UndoDelBox->setCheckState((Qt::CheckState)0);	
       this->AddBox->setCheckState((Qt::CheckState)2);
       this->mode =1;
       break;
                            }
   default:
       // should never be reached
       break;
 }
	}
}

else{
this->SplitBox->setCheckState((Qt::CheckState)this->stateSplit);
if(this->stateSplit){
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }
if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
if(this->stateMerge) { this->MergeBox->setCheckState((Qt::CheckState)0); }
if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
this->mode = 3;
		    }

else{
this->mode=0;
}
}
}

//Merge Seeds
void Seed3D::MergeSeeds()
{
if(this->stateUndoDel) {this->UndoDelBox->setCheckState((Qt::CheckState)0);}
this->mode =4;
this->stateMerge = this->MergeBox->checkState();
if(MarkedPoints.size()!=0 && this->stateMerge){
	if(this->stateDelete || this->stateUndoDel){
this->msgBox = new QMessageBox(this);
this->msgBox->setText("You have seeds marked for deletion. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();
switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   this->mode =2;
   Apply();
   this->DeleteBox->setCheckState((Qt::CheckState)0); 
   this->mode =4;
   //If all seeds are deleted no more seeds can be merged
   //Should not allow the user to click merge
   if(this->dup_points.size()!=0)
   this->MergeBox->setCheckState((Qt::CheckState)2);
   else
   this->MergeBox->setCheckState((Qt::CheckState)0); 
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //DeleteSeed(); 
       this->MergeBox->setCheckState((Qt::CheckState)0);
       this->UndoDelBox->setCheckState((Qt::CheckState)0);
       this->DeleteBox->setCheckState((Qt::CheckState)2);
              this->mode =2;
	   break;
                            }
   default:
       // should never be reached
       break;
 }
	}
}

else if(this->MarkedPoints2add.size()!=0 && this->stateMerge){
this->msgBox = new QMessageBox(this);
if(this->stateAdd)
this->msgBox->setText("You have seeds marked for addition. Do you want to apply the changes?");
if(this->stateSplit)
this->msgBox->setText("You have seeds marked for splitting. Do you want to apply the changes?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   if(this->stateAdd)
   this->mode =1;
   if(this->stateSplit)
   this->mode =3;
   Apply();
   if(this->stateAdd)
   this->AddBox->setCheckState((Qt::CheckState)0);
   if(this->stateSplit)
   this->SplitBox->setCheckState((Qt::CheckState)0); 
   this->mode =4;
   this->MergeBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       //AddSeed(); 
       
	   if(this->stateAdd){
	   
	   this->MergeBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);
           this->AddBox->setCheckState((Qt::CheckState)2);
	   this->mode =1;}
	   if(this->stateSplit){
	   this->MergeBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);
           this->SplitBox->setCheckState((Qt::CheckState)2);
	   this->mode =3;
	   }
	   
       break;
                            }
   default:
       // should never be reached
       break;
 }

}

else{
this->MergeBox->setCheckState((Qt::CheckState)this->stateMerge);

if(this->stateMerge){
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }
if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
if(this->stateSplit) { this->SplitBox->setCheckState((Qt::CheckState)0); }
this->UndoDelBox->setCheckState((Qt::CheckState)0); 
this->mode = 4;
		            }
else{
this->mode=0;
    }
}
}

void Seed3D::UndoDeleteSeeds()
{
this->stateUndoDel = this->UndoDelBox->checkState();
this->UndoDelBox->setCheckState((Qt::CheckState)this->stateUndoDel);
if(this->stateUndoDel){
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }
if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
if(this->stateSplit) { this->SplitBox->setCheckState((Qt::CheckState)0); }
if(this->stateMerge) { this->MergeBox->setCheckState((Qt::CheckState)0); }
this->UndoDelBox->setCheckState((Qt::CheckState)2);
this->stateUndoDel = this->UndoDelBox->checkState();
this->mode = 5;
}
else{
this->mode=0;
}

}


void Seed3D::Apply()
{
//Put added seeds in a vector for Yousef
//Also add it to the main vector in case we want 
//to delete the added seeds
 
if(this->mode==1)
	{
	vtkIdType Id;
	double* pointz;
	point p;
	Id = this->point3->GetNumberOfPoints();
	std::vector<point> tobeAdded;
	
	for(int i =0 ; i<Id;i++)
		{
    			pointz = this->point3->GetPoint(i);
				this->point1->InsertNextPoint(pointz);
				this->polydata1->SetPoints(this->point1);
				this->Glyph->SetInput(this->polydata1);	
				p.x = pointz[0];
    			p.y = pointz[1];
    			p.z = pointz[2];	
    			tobeAdded.push_back(p);
				this->dup_points.push_back(p);
		}


//Remove the green glyph		
        vtkDataArray* points2add = this->point3->GetData();    
        points2add->Reset();
        points2add->Squeeze();
        this->point3->SetData(points2add);
	    this->addglyph->SetScaleFactor(this->addglyph->GetScaleFactor()+0.0001);
        this->SphereMapper->SetInput(this->Glyph->GetOutput());	      
        this->Glyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        this->QVTK->GetRenderWindow()->Render();
		this->MarkedPoints2add.erase(this->MarkedPoints2add.begin(),this->MarkedPoints2add.end());
 
}


//Delete the seeds from the screen
if(this->mode==2)
        { 
this->Glyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
vtkIdType Id;
double* pointz;
point p;
Id = this->point2->GetNumberOfPoints();
std::vector<point> tobeDeleted;
cout<<Id<<endl;
for(int i =0 ;i<Id;i++)
{
    pointz = this->point2->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    tobeDeleted.push_back(p);
    
}

vtkDataArray* points2del = this->point2->GetData();    
//Remove the glyph		

for(int i =0 ;i<Id;i++)
{		
    points2del->RemoveTuple((vtkIdType)0);
    cout<<"deleted red"<<endl;		
}
    this->point2->SetData(points2del);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);    
    this->QVTK->GetRenderWindow()->Render();
    this->MarkedPoints.erase(this->MarkedPoints.begin(),this->MarkedPoints.end());
    //this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	
}

//Split Nuclei
if(this->mode==3)
	{
	vtkIdType Id;
	double* pointz;
	point p;
	Id = this->point3->GetNumberOfPoints();
	std::vector<point> tobeSplit;
	
	for(int i =0 ; i<Id;i++)
		{
    			pointz = this->point3->GetPoint(i);
				this->point1->InsertNextPoint(pointz);
				this->polydata1->SetPoints(this->point1);
				this->Glyph->SetInput(this->polydata1);	
				p.x = pointz[0];
    			p.y = pointz[1];
    			p.z = pointz[2];	
    			tobeSplit.push_back(p);
				this->dup_points.push_back(p);
		}


//Remove the green glyph		
        vtkDataArray* points2add = this->point3->GetData();    
        points2add->Reset();
        points2add->Squeeze();
        this->point3->SetData(points2add);
	    this->addglyph->SetScaleFactor(this->addglyph->GetScaleFactor()+0.0001);
        this->SphereMapper->SetInput(this->Glyph->GetOutput());	      
        this->Glyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
        this->QVTK->GetRenderWindow()->Render();
		this->MarkedPoints2add.erase(this->MarkedPoints2add.begin(),this->MarkedPoints2add.end());
 
}




//Delete the seeds from the screen
if(this->mode==4)
        { 
this->Glyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
vtkIdType Id;
double* pointz;
point p;
Id = this->point2->GetNumberOfPoints();
std::vector<point> tobeMerged;
for(int i =0 ;i<Id;i++)
{
    pointz = this->point2->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    tobeMerged.push_back(p);
    
}

vtkDataArray* points2merge = this->point2->GetData();    
//Remove the glyph		

for(int i =0 ;i<Id;i++)
{		
    points2merge->RemoveTuple((vtkIdType)0);
    cout<<"merge red"<<endl;		
}
    this->point2->SetData(points2merge);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);    
    this->QVTK->GetRenderWindow()->Render();
    this->MarkedPoints.erase(this->MarkedPoints.begin(),this->MarkedPoints.end());
    //this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	
}

Check();	

}

void Seed3D::Check()
{
if(this->dup_points.size()==0)
{

//this->DeleteBox->setCheckState((Qt::CheckState)0);
this->SplitBox->setCheckState((Qt::CheckState)0);
//this->MergeBox->setCheckState((Qt::CheckState)0);
//this->DeleteBox->setCheckable(false);
this->SplitBox->setCheckable(false);
//this->MergeBox->setCheckable(false);

}

else
{
//this->DeleteBox->setCheckable(true);
this->SplitBox->setCheckable(true);
//this->MergeBox->setCheckable(true);

}


//Atleast one marked point required to undo selection for delete
if(this->MarkedPoints.size()==0){
this->UndoDelBox->setCheckState((Qt::CheckState)0);
this->UndoDelBox->setCheckable(false);
}
else{
this->UndoDelBox->setCheckable(true);
}
}




