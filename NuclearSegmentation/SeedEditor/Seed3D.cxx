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

#include "Seed3D.h"
#include "Seed3DHelperClasses.h"

//*******************************************************************************
// SeedEditor
// This widget is a main window for editing nuclei.  
//********************************************************************************

Seed3D::Seed3D(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{

	//Create the GUI
	createMenus();
	createStatusBar();
	setWindowTitle(tr("Seed Editor"));
	lastPath = ".";

//Variable instantiation
	this->mode = 0;
	this->counter =0;
	this->stateAdd = 0;
	this->stateDelete = 0;
	this->stateTP = 0;
	this->stateCP = 0;
	this->stateFP = 0;
	this->stateFN = 0;
	
	this->flag =1;
	
	
	// QVTK  -> QVTK Widget for the main window displaying the volume.
	// QVTK1 -> QVTK1 Widget for the window displaying the YZ plane.
	// QVTK2 -> QVTK2 Widget for the window displaying the XY plane.
	
	this->QVTK = 0;
	this->QVTK1 = 0;
	this->QVTK2 = 0;
    
    //This variable is used to indicate if an image has already been loaded
	this->iRender=0;

  //Set up the Widgets & Renderers to display the actors and the Volume	
  this->browse = new QWidget(this);
  this->setCentralWidget(browse);
  this->QVTK = new QVTKWidget(this->browse);
  this->Renderer = vtkRenderer::New();
  this->QVTK->GetRenderWindow()->AddRenderer(this->Renderer);
  this->QVTK1 = new QVTKWidget(this->browse);
  this->Renderer1 = vtkRenderer::New();
  this->QVTK1->GetRenderWindow()->AddRenderer(this->Renderer1);
  this->QVTK2 = new QVTKWidget(this->browse);
  this->Renderer2 = vtkRenderer::New();
  this->QVTK2->GetRenderWindow()->AddRenderer(this->Renderer2);

  
  
  
  
  //Setup the buttons that the user will use to interact with this program. 
  
  this->EditRbutton = new QRadioButton("Edit Mode",this->browse);
  this->ValidateRbutton = new QRadioButton("Validation Mode",this->browse);
  this->AddBox = new QCheckBox("&Add Seeds", this->browse);
  this->DeleteBox = new QCheckBox("&Delete Seeds", this->browse);
  this->UndoDelBox = new QCheckBox("&Undo Selections for Delete", this->browse);
  this->EditRbutton->click();
  
  // Validation checkboxes
  this->TPBox = new QCheckBox("&True Positives", this->browse);
  this->CPBox = new QCheckBox("&Cluster Positives", this->browse);
  this->FPBox = new QCheckBox("&False Positives", this->browse);
  this->FNBox = new QCheckBox("&False Negatives", this->browse);
    
  this->PlaceButton = new QPushButton("&Place Seed", this->browse);
  this->ApplyButton = new QPushButton("&Apply", this->browse);  
  this->ValidateButton = new QPushButton("&Validate", this->browse);  
//********************************************************************************
// Layouts 

  QGridLayout *buttonLayout = new QGridLayout(this);
  
  buttonLayout->addWidget(this->EditRbutton,10,1);
  buttonLayout->addWidget(this->ValidateRbutton,10,2);
  buttonLayout->addWidget(this->AddBox,11, 1);
  buttonLayout->addWidget(this->DeleteBox,11,2);
  buttonLayout->addWidget(this->UndoDelBox,11,3);
  
  buttonLayout->addWidget(this->TPBox,13,1);
  buttonLayout->addWidget(this->CPBox,13,2);
  buttonLayout->addWidget(this->FPBox,13,3);
  buttonLayout->addWidget(this->FNBox,13,4);
  	
  buttonLayout->addWidget(this->PlaceButton,12, 1);
  buttonLayout->addWidget(this->ApplyButton,12, 2);
  buttonLayout->addWidget(this->ValidateButton,14, 1);	


  QGridLayout *viewerLayout = new QGridLayout(this->browse);
  viewerLayout->addWidget(this->QVTK, 0,0,1,2);
  viewerLayout->addWidget(this->QVTK1,1,0);
  viewerLayout->addWidget(this->QVTK2,1,1); 	
  viewerLayout->addLayout(buttonLayout,2,0); 

  //********************************************************************************
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
  this->Interactor->SetPicker(this->PointPicker);
  this->isPicked = vtkCallbackCommand::New();
  this->isPicked->SetCallback(PickCell);

  //isPicked caller allows observer to intepret click 
  this->isPicked->SetClientData(this);            
  this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,isPicked); 
  this->Interactor1->SetInteractorStyle(style1);
  this->Interactor2->SetInteractorStyle(style1);
  

  //What happens when the user clicks the buttons or the Checkboxes is determined by the "SLOTS"
  connect(this->PlaceButton, SIGNAL(clicked()), this, SLOT(PlaceSeed()));
  connect(this->ValidateButton, SIGNAL(clicked()), this, SLOT(ValidateSeeds()));
  connect(this->AddBox, SIGNAL(stateChanged(int)), this, SLOT(AddSeed()));
  connect(this->DeleteBox, SIGNAL(stateChanged(int)), this, SLOT(DeleteSeed()));
  connect(this->ApplyButton, SIGNAL(clicked()), this, SLOT(Apply()));		
  connect(this->UndoDelBox, SIGNAL(stateChanged(int)), this, SLOT(UndoDeleteSeeds()));
  connect(this->TPBox, SIGNAL(stateChanged(int)), this, SLOT(TruePositives()));
  connect(this->CPBox, SIGNAL(stateChanged(int)), this, SLOT(ClusterPositives()));
  connect(this->FPBox, SIGNAL(stateChanged(int)), this, SLOT(FalsePositives()));
  connect(this->FNBox, SIGNAL(stateChanged(int)), this, SLOT(FalseNegatives()));

//Resize the Window 
  this->resize(800,800);
}


void Seed3D::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));
	loadAction = new QAction(tr("Load Image and Seed File..."), this);
	loadAction->setStatusTip(tr("Load an image into the browser"));
	connect(loadAction, SIGNAL(triggered()), this, SLOT(loadImage()));
	fileMenu->addAction(loadAction); 
	fileMenu->addSeparator();

	saveAction = new QAction(tr("Save Seeds"), this);
	saveAction->setStatusTip(tr("Save Changes (Edits, etc)"));
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(saveResult()));
	fileMenu->addAction(saveAction);
	fileMenu->addSeparator();

    exitAction = new QAction(tr("Exit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
    fileMenu->addAction(exitAction);

}

void Seed3D::createStatusBar()
{
    statusLabel = new QLabel(tr(" Ready"));
    statusBar()->addWidget(statusLabel, 1);
}



// The core logic is present in this function
void Seed3D::loadImage()
{  
   // Browser to open the Image
   this->fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							    "All Files (*.*)"));
  if(this->fileName == "")
	return;
  this->lastPath = QFileInfo(this->fileName).absolutePath();
  

//Seeds
this->fileNameSeed = QFileDialog::getOpenFileName(
                             this, "Select seed file to open", lastPath,
                             tr("Text files (*.txt);; \n"
							    "All Files (*.*)"));
  if(this->fileNameSeed == "")
	return;
  QString name = QFileInfo(this->fileNameSeed).baseName();	
	
// if image already exists delete the objects

  if(this->iRender==1){
	  DeleteObjects();
					  }

//*******************************************************************************

//Read the image
  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;
  ReaderType::Pointer i2spReader = ReaderType::New();
  i2spReader->SetFileName( fileName.toStdString() );
  try
  {
    i2spReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    }
  
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
    

 // Create seeds (spheres) only at coordinates specified by the text file.
 // I also create Red (Markedpoints) and Green Spheres (Markedpoints2add) initially and delete them.
//  So that the first delete and addition by the user would be 
// "inserts" into the vector !
// Note: 1) All the seeds form a single glyph.
//       2) The green spheres form a single glyph. These are seeds added but not applied permanently.
//		 3) The red spheres form a single glyph.These seeds are marked for Deletion but can be unmarked
	
	std::vector<point> spPoint;
    std::vector<point> MarkedPoints;
    std::vector<point> MarkedPoints2add;
	std::vector<point> dup_points; //used while deleteing seeds
    	
	std::string sfilename = fileNameSeed.toStdString();
	char* seedname;
	seedname = &sfilename[0];

	//spPoint is the vector of seeds 
	spPoint = ReadPoints(seedname);
    this->dup_points = spPoint; //dup_points contains the upto date location of seeds 
		
    // Create a float array which represents the points.
    this->pcoords = vtkFloatArray::New();
    this->delpcoords = vtkFloatArray::New();	
    this->TPpcoords = vtkFloatArray::New();	
	this->CPpcoords = vtkFloatArray::New();	
	this->FPpcoords = vtkFloatArray::New();	
	this->FNpcoords = vtkFloatArray::New();	




	// Note that by default, an array has 1 component.
    // We have to change it to 3 for points
    this->pcoords->SetNumberOfComponents(3);
    this->pcoords->SetNumberOfTuples(spPoint.size());
    this->delpcoords->SetNumberOfComponents(3);
    this->delpcoords->SetNumberOfTuples(1);		
	this->TPpcoords->SetNumberOfComponents(3);
    this->TPpcoords->SetNumberOfTuples(1);
 	this->CPpcoords->SetNumberOfComponents(3);
    this->CPpcoords->SetNumberOfTuples(1);
	this->FPpcoords->SetNumberOfComponents(3);
    this->FPpcoords->SetNumberOfTuples(1);
	this->FNpcoords->SetNumberOfComponents(3);
    this->FNpcoords->SetNumberOfTuples(1);



  for (unsigned int j=0; j<spPoint.size(); j++)
    {
    float pts[3] = {spPoint[j].x,spPoint[j].y , spPoint[j].z };	
    this->pcoords->SetTuple(j, pts);
    }
    float pts[3] = {0.0,0.0,0.0};
	this->delpcoords->SetTuple(0,pts);	
    this->TPpcoords->SetTuple(0,pts);	
	this->CPpcoords->SetTuple(0,pts);	
	this->FPpcoords->SetTuple(0,pts);	
	this->FNpcoords->SetTuple(0,pts);	
// Create vtkPoints and assign pcoords as the internal data array.
  this->point1 = vtkPoints::New(); //Original Points
  this->point2 = vtkPoints::New(); //Seeds marked for deletion
  this->point3 = vtkPoints::New(); //Seeds marked for addition
  this->pointTP = vtkPoints::New();
  this->pointCP = vtkPoints::New();
  this->pointFP = vtkPoints::New();
  this->pointFN = vtkPoints::New();
  this->allpoints = vtkPoints::New();//Yousefseg downstream
  
  this->point1->SetData(this->pcoords);	  
  this->point2->SetData(this->delpcoords);
  this->point3->SetData(this->delpcoords);
  this->pointTP->SetData(this->TPpcoords);	
  this->pointCP->SetData(this->CPpcoords);	
  this->pointFP->SetData(this->FPpcoords);	
  this->pointFN->SetData(this->FNpcoords);	

  // Assign points and cells
  this->polydata1 = vtkPolyData::New();
  this->polydata1->SetPoints(this->point1);
  this->polydata2 = vtkPolyData::New();	
  this->polydata2->SetPoints(this->point2);
  
  this->polydata3 = vtkPolyData::New();	
  this->polydata3->SetPoints(this->point3);

  this->polydataTP = vtkPolyData::New();	
  this->polydataTP->SetPoints(this->pointTP);
  this->polydataCP = vtkPolyData::New();	
  this->polydataCP->SetPoints(this->pointCP);
  this->polydataFP = vtkPolyData::New();	
  this->polydataFP->SetPoints(this->pointFP);
  this->polydataFN = vtkPolyData::New();	
  this->polydataFN->SetPoints(this->pointFN);

// create the glyphing the sphere with a cone.  Create the mapper
// and actor for the glyphs.

    vtkSphereSource *sphere = vtkSphereSource::New();
    vtkGlyph3D *glyph = vtkGlyph3D::New();
    vtkGlyph3D *delglyph = vtkGlyph3D::New();
    vtkGlyph3D *addglyph = vtkGlyph3D::New();
    vtkGlyph3D *TPglyph = vtkGlyph3D::New();
    vtkGlyph3D *CPglyph = vtkGlyph3D::New();
	vtkGlyph3D *FPglyph = vtkGlyph3D::New();
	vtkGlyph3D *FNglyph = vtkGlyph3D::New();




	
	glyph->SetInput(this->polydata1);
    delglyph->SetInput(this->polydata2);
    addglyph->SetInput(this->polydata3);
    TPglyph->SetInput(this->polydataTP);
	CPglyph->SetInput(this->polydataCP);
	FPglyph->SetInput(this->polydataFP);
	FNglyph->SetInput(this->polydataFN);

	glyph->SetSource(sphere->GetOutput());
    delglyph->SetSource(sphere->GetOutput());
    addglyph->SetSource(sphere->GetOutput());
    TPglyph->SetSource(sphere->GetOutput());
	CPglyph->SetSource(sphere->GetOutput());
	FPglyph->SetSource(sphere->GetOutput());
	FNglyph->SetSource(sphere->GetOutput());



	glyph->SetVectorModeToUseNormal();
    delglyph->SetVectorModeToUseNormal();
    addglyph->SetVectorModeToUseNormal();
    TPglyph->SetVectorModeToUseNormal();
	CPglyph->SetVectorModeToUseNormal();
	FPglyph->SetVectorModeToUseNormal();
	FNglyph->SetVectorModeToUseNormal();


	glyph->SetScaleModeToScaleByVector();
    delglyph->SetScaleModeToScaleByVector();
    addglyph->SetScaleModeToScaleByVector();
    TPglyph->SetScaleModeToScaleByVector();
	CPglyph->SetScaleModeToScaleByVector();
	FPglyph->SetScaleModeToScaleByVector();
	FNglyph->SetScaleModeToScaleByVector();



	glyph->SetScaleFactor(1);
    delglyph->SetScaleFactor(1);
    addglyph->SetScaleFactor(1);
    TPglyph->SetScaleFactor(1);
	CPglyph->SetScaleFactor(1);
	FPglyph->SetScaleFactor(1);
	FNglyph->SetScaleFactor(1);



	glyph->GeneratePointIdsOn();
    delglyph->GeneratePointIdsOn();
    addglyph->GeneratePointIdsOn();
    TPglyph->GeneratePointIdsOn();
	CPglyph->GeneratePointIdsOn();
	FPglyph->GeneratePointIdsOn();
	FNglyph->GeneratePointIdsOn();

	this->Glyph = glyph;
    this->delglyph = delglyph;
    this->addglyph = addglyph;
    this->TPglyph = TPglyph;
	this->CPglyph = CPglyph;
	this->FPglyph = FPglyph;
	this->FNglyph = FNglyph;
	
	
	this->imActor = imActor;
    

	// The Pipeline
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper *delsphereMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper *addsphereMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper *TPsphereMapper = vtkPolyDataMapper::New();
	vtkPolyDataMapper *CPsphereMapper = vtkPolyDataMapper::New();
	vtkPolyDataMapper *FPsphereMapper = vtkPolyDataMapper::New();
	vtkPolyDataMapper *FNsphereMapper = vtkPolyDataMapper::New();

	sphereMapper->SetInput(this->Glyph->GetOutput());
    delsphereMapper->SetInput(this->delglyph->GetOutput());
    addsphereMapper->SetInput(this->addglyph->GetOutput());
    TPsphereMapper->SetInput(this->TPglyph->GetOutput());
	CPsphereMapper->SetInput(this->CPglyph->GetOutput());
	FPsphereMapper->SetInput(this->FPglyph->GetOutput());
	FNsphereMapper->SetInput(this->FNglyph->GetOutput());



	vtkLODActor *sphereActor = vtkLODActor::New();
    vtkLODActor *delsphereActor = vtkLODActor::New();
    vtkLODActor *addsphereActor = vtkLODActor::New();
    vtkLODActor *TPsphereActor = vtkLODActor::New();
	vtkLODActor *CPsphereActor = vtkLODActor::New();
	vtkLODActor *FPsphereActor = vtkLODActor::New();
	vtkLODActor *FNsphereActor = vtkLODActor::New();

	sphereActor->SetMapper(sphereMapper);
    delsphereActor->SetMapper(delsphereMapper);
    delsphereActor->GetProperty()->SetColor(0.5,0.0,0.0); 	
    addsphereActor->SetMapper(addsphereMapper);
    addsphereActor->GetProperty()->SetColor(0.0,0.5,0.0);
    
	TPsphereActor->SetMapper(TPsphereMapper);
    TPsphereActor->GetProperty()->SetColor(0.0,0.0,0.6);

	CPsphereActor->SetMapper(CPsphereMapper);
    CPsphereActor->GetProperty()->SetColor(0.0,0.8,0.6);

	FPsphereActor->SetMapper(FPsphereMapper);
    FPsphereActor->GetProperty()->SetColor(0.4,0.0,0.6);

	FNsphereActor->SetMapper(FNsphereMapper);
    FNsphereActor->GetProperty()->SetColor(0.8,0.1,0.0);

	
	Renderer->AddVolume(volume);
    Renderer->AddActor(sphereActor);
    Renderer->AddActor(delsphereActor);
    Renderer->AddActor(addsphereActor);	
	Renderer->AddActor(TPsphereActor);	
	Renderer->AddActor(CPsphereActor);
	Renderer->AddActor(FPsphereActor);
	Renderer->AddActor(FNsphereActor);

    this->SphereMapper = sphereMapper;
    this->SphereActor = sphereActor;
    this->DelSphereMapper = delsphereMapper;
    this->DelSphereActor = delsphereActor;		
    this->AddSphereMapper = addsphereMapper;
    this->AddSphereActor = addsphereActor;
    this->TPSphereActor = TPsphereActor;
	this->TPsphereMapper = TPsphereMapper;
	this->CPSphereActor = CPsphereActor;
	this->CPsphereMapper = CPsphereMapper;
	this->FPSphereActor = FPsphereActor;
	this->FPsphereMapper = FPsphereMapper;
	this->FNSphereActor = FNsphereActor;
	this->FNsphereMapper = FNsphereMapper;



	// The active camera for the main window!
	vtkCamera *cam1 = Renderer->GetActiveCamera();
    cam1->SetViewUp (0, 1, 0);
    cam1->SetPosition (0, 0, 50);
    Renderer->ResetCamera();	
    this->Volume = volume;
   
    //Remove the initial red and green glyphs. 
    vtkDataArray* delpoint = this->point2->GetData();    
    delpoint->RemoveTuple((vtkIdType)0);
    this->point2->SetData(delpoint);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);
    
    vtkDataArray* addpoint = this->point3->GetData();    
    addpoint->RemoveTuple((vtkIdType)0);
    this->point3->SetData(addpoint);
    this->addglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	


	vtkDataArray* TPpoint = this->pointTP->GetData();    
    TPpoint->RemoveTuple((vtkIdType)0);
    this->pointTP->SetData(TPpoint);
    this->TPglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	

	vtkDataArray* CPpoint = this->pointCP->GetData();    
    CPpoint->RemoveTuple((vtkIdType)0);
    this->pointCP->SetData(CPpoint);
    this->TPglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	


	vtkDataArray* FPpoint = this->pointFP->GetData();    
    FPpoint->RemoveTuple((vtkIdType)0);
    this->pointFP->SetData(FPpoint);
    this->FPglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	


	vtkDataArray* FNpoint = this->pointFN->GetData();    
    FNpoint->RemoveTuple((vtkIdType)0);
    this->pointFN->SetData(FNpoint);
    this->FNglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	


//********************************************************************************
// This section of code deals with the creation of the handles and handle widget

  // Create the widget and its representation
  vtkPointHandleRepresentation3D *handle = vtkPointHandleRepresentation3D::New();
  vtkHandleWidget *widget = vtkHandleWidget::New();
  widget->SetInteractor(this->Interactor);
  widget->SetRepresentation(handle);
  this->widget = widget;
  this->handle = handle; 
  vtkSeedCallback *scbk = vtkSeedCallback::New();
  widget->AddObserver(vtkCommand::EndInteractionEvent,scbk);
  widget->SetEnabled(true);

//********************************************************************************
// Create the Imagereslice planes and the actors to the respectove Qwidgets

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
  
  // Set the extent for the Y-Z slice.
  double* origincalc = this->Volume->GetBounds();	
  int origin[6];
  origin[0] = (int)origincalc[2];
  origin[1] = (int)origincalc[3];
  origin[2] = (int)origincalc[4];
  origin[3] = (int)origincalc[5];
  origin[4] = (int)origincalc[0];
  origin[5] = (int)origincalc[1];

  reslice->SetOutputExtent(origin); 	 
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
  this->imActor->RotateWXYZ(180,0,0,1);

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

//********************************************************************************  
// handle  -> handle for the main window
// handle1 -> handle for the Y-Z plane actor (imActor formed by reslice)
// handle2 -> handle for the X-Y plane actor (imActor formed by reslice1)

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
  scbk1->QVTK1 = this->QVTK1;
  scbk1->QVTK = this->QVTK;	
  scbk1->vol = this->Volume;
  scbk1->handle = this->handle; 
  //scbk1->reslice = this->Reslice;
  scbk1->reslice1 = this->Reslice1;
  scbk1->imActor1 = this->imActor1; 
  scbk2->handle1= this->handle1;	
  scbk2->handle2= this->handle2;
  scbk2->QVTK2 = this->QVTK2;	
  scbk2->QVTK1 = this->QVTK1;
  scbk2->QVTK = this->QVTK;	
  scbk2->vol = this->Volume;
  scbk2->handle = this->handle; 
  scbk2->reslice = this->Reslice;
  scbk2->imActor = this->imActor;
  
  
//********************************************************************************  
//	These are the sliders used to control the opacity, brightness and seed size !

//OPACITY
  this->sliderRep = vtkSliderRepresentation2D::New();
  this->sliderRep->SetValue(0.8);
  this->sliderRep->SetTitleText("Opacity");
  this->sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep->GetPoint1Coordinate()->SetValue(0.2,0.1);
  this->sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep->GetPoint2Coordinate()->SetValue(0.8,0.1);
  this->sliderRep->SetSliderLength(0.02);
  this->sliderRep->SetSliderWidth(0.03);
  this->sliderRep->SetEndCapLength(0.01);
  this->sliderRep->SetEndCapWidth(0.03);
  this->sliderRep->SetTubeWidth(0.005);
  this->sliderRep->SetMinimumValue(0.0);
  this->sliderRep->SetMaximumValue(1.0);

  this->sliderWidget = vtkSliderWidget::New();
  this->sliderWidget->SetInteractor(Interactor);
  this->sliderWidget->SetRepresentation(this->sliderRep);
  this->sliderWidget->SetAnimationModeToAnimate();
  vtkSlider2DCallbackBrightness *callback_brightness = vtkSlider2DCallbackBrightness::New();
  callback_brightness->volume = this->Volume;
  this->sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback_brightness);
  this->sliderWidget->EnabledOn();
    
//Brightness
  this->sliderRep2 = vtkSliderRepresentation2D::New();
  this->sliderRep2->SetValue(0.8);
  this->sliderRep2->SetTitleText("Brightness");
  this->sliderRep2->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep2->GetPoint1Coordinate()->SetValue(0.2,0.9);
  this->sliderRep2->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep2->GetPoint2Coordinate()->SetValue(0.8,0.9);
  this->sliderRep2->SetSliderLength(0.02);
  this->sliderRep2->SetSliderWidth(0.03);
  this->sliderRep2->SetEndCapLength(0.01);
  this->sliderRep2->SetEndCapWidth(0.03);
  this->sliderRep2->SetTubeWidth(0.005);
  this->sliderRep2->SetMinimumValue(0.0);
  this->sliderRep2->SetMaximumValue(1.0);

  this->sliderWidget2 = vtkSliderWidget::New();
  this->sliderWidget2->SetInteractor(Interactor);
  this->sliderWidget2->SetRepresentation(this->sliderRep2);				
  this->sliderWidget2->SetAnimationModeToAnimate();
  vtkSlider2DCallbackContrast *callback_contrast = vtkSlider2DCallbackContrast::New();
  callback_contrast->volume = this->Volume;
  this->sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contrast);
  this->sliderWidget2->EnabledOn();



//seed size
  this->sliderRep3 = vtkSliderRepresentation2D::New();
  this->sliderRep3->SetValue(0.8);
  this->sliderRep3->SetTitleText("Seed Size");
  this->sliderRep3->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep3->GetPoint1Coordinate()->SetValue(0.1,0.2);
  this->sliderRep3->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  this->sliderRep3->GetPoint2Coordinate()->SetValue(0.1,0.9);
  this->sliderRep3->SetSliderLength(0.01);
  this->sliderRep3->SetSliderWidth(0.03);
  this->sliderRep3->SetEndCapLength(0.01);
  this->sliderRep3->SetEndCapWidth(0.03);
  this->sliderRep3->SetTubeWidth(0.005);
  this->sliderRep3->SetMinimumValue(1.0);
  this->sliderRep3->SetMaximumValue(15.0);

  this->sliderWidget3 = vtkSliderWidget::New();
  this->sliderWidget3->SetInteractor(Interactor);
  this->sliderWidget3->SetRepresentation(this->sliderRep3);
  this->sliderWidget3->SetAnimationModeToAnimate();
  

  vtkSlider2DCallbackSeedSize *callback_seedsize = vtkSlider2DCallbackSeedSize::New();
  callback_seedsize->Glyph = this->Glyph;
  callback_seedsize->addglyph = this->addglyph;
  callback_seedsize->delglyph = this->delglyph;	
  callback_seedsize->TPglyph = this->TPglyph;	
  callback_seedsize->CPglyph = this->CPglyph;	
  callback_seedsize->FPglyph = this->FPglyph;	
  callback_seedsize->FNglyph = this->FNglyph;	




  this->sliderWidget3->AddObserver(vtkCommand::InteractionEvent,callback_seedsize);
  this->sliderWidget3->EnabledOn();

  this->QVTK->GetRenderWindow()->Render();   
  this->QVTK1->GetRenderWindow()->Render();   
  this->QVTK2->GetRenderWindow()->Render(); 
  this->iRender = 1;	 

}

void   Seed3D::PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{ 
   Seed3D* seed = (Seed3D*)clientdata;
 /*  PickPoint allows fot the point id and coordinates to be returned 
  as well as adding a marker on the last picked point
  R_click to select point on line  */
	
 int *pos = seed->Interactor->GetEventPosition();
 seed->Interactor->GetPicker()->Pick(pos[0],pos[1],0.0,seed->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
 //vtkPointPicker *point_picker = (vtkPointPicker *)seed->Interactor->GetPicker();
 double pickPos[3];
 seed->PointPicker->GetPickPosition(pickPos);    //this is the coordinates of the pick  
 
//Useful to check if clicked on a seed !   
//If I am in Delete mode 

if((seed->mode == 2) && seed->flag == 1){

vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 
	
if((unsigned int)pID<=seed->dup_points.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    //float dist1;
    int index;
    float finpt[3];
    for (unsigned int j=0; j<seed->dup_points.size(); j++)
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
	////cout<<index<<" is the inD value"<<endl;
	seed->dup_points.erase(seed->dup_points.begin()+index);

	//Remove the glyph		
     vtkDataArray* points2del = seed->point1->GetData();    
     //vtkDataArray* points2delred;
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

if(seed->MarkedPoints.size()==seed->dup_points.size() || seed->dup_points.size()==0) 
{
seed->flag =0; 	
}
else{
seed->flag =1; 	
   }



//Useful to check if clicked on a seed !   
//If I am in Validation TP mode 
if(seed->mode == 6) {


vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 

if((unsigned int)pID<=seed->dup_points.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    //float dist1;
    int index;
    float finpt[3];
    for (unsigned int j=0; j<seed->dup_points.size(); j++)

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
	
	seed->dup_points.erase(seed->dup_points.begin()+index);

	//Remove the glyph		
     vtkDataArray* points2del = seed->point1->GetData();    
     //vtkDataArray* points2delred;
     points2del->RemoveTuple((vtkIdType)index);
     seed->point1->SetData(points2del);
     seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);


	  // Add the new blue glyph 
       seed->pointTP->InsertNextPoint(finpt);
	   seed->polydataTP->SetPoints(seed->pointTP);
	   cout<<seed->polydataTP->GetNumberOfPoints()<<"jjj"<<endl;      
	   seed->TPglyph->SetInput(seed->polydataTP);
       seed->TPsphereMapper->SetInput(seed->TPglyph->GetOutput());	      
       seed->TPglyph->SetScaleFactor(seed->TPglyph->GetScaleFactor()+0.0001);//to rerender immediately
       seed->QVTK->GetRenderWindow()->Render();
       cout<<"Add the new blue glyph"<<endl;

}
}



if(seed->mode == 7) {


vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 

if((unsigned int)pID<=seed->dup_points.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    //float dist1;
    int index;
    float finpt[3];
    for (unsigned int j=0; j<seed->dup_points.size(); j++)

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
	
	seed->dup_points.erase(seed->dup_points.begin()+index);

	//Remove the glyph		
     vtkDataArray* points2del = seed->point1->GetData();    
     //vtkDataArray* points2delred;
     points2del->RemoveTuple((vtkIdType)index);
     seed->point1->SetData(points2del);
     seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);



     // Add the new  glyph 
       seed->pointCP->InsertNextPoint(finpt);
	   cout<<seed->pointCP->GetNumberOfPoints()<<"jjj"<<endl; 
	   cout<<seed->pointTP->GetNumberOfPoints()<<"jjj"<<endl; 	
	   seed->polydataCP->SetPoints(seed->pointCP);
       seed->CPglyph->SetInput(seed->polydataCP);
       seed->CPsphereMapper->SetInput(seed->CPglyph->GetOutput());	      
       
	   seed->CPglyph->SetScaleFactor(seed->CPglyph->GetScaleFactor()+0.0001);//to rerender immediately
       seed->QVTK->GetRenderWindow()->Render();
       cout<<seed->polydataTP->GetNumberOfPoints()<<"jjj"<<endl; 

}
}


if(seed->mode == 8) {


vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 

if((unsigned int)pID<=seed->dup_points.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    //float dist1;
    int index;
    float finpt[3];
    for (unsigned int j=0; j<seed->dup_points.size(); j++)

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
	
	seed->dup_points.erase(seed->dup_points.begin()+index);

	//Remove the glyph		
     vtkDataArray* points2del = seed->point1->GetData();    
     //vtkDataArray* points2delred;
     points2del->RemoveTuple((vtkIdType)index);
     seed->point1->SetData(points2del);
     seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);



     // Add the new blue glyph 
       seed->pointFP->InsertNextPoint(finpt);
	   seed->polydataFP->SetPoints(seed->pointFP);
       seed->FPglyph->SetInput(seed->polydataFP);
       seed->FPsphereMapper->SetInput(seed->FPglyph->GetOutput());	      
       seed->FPglyph->SetScaleFactor(seed->FPglyph->GetScaleFactor()+0.0001);//to rerender immediately
       seed->QVTK->GetRenderWindow()->Render();
       cout<<"Add the new glyph"<<endl;
}
}



if(seed->mode == 9) {


vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 

if((unsigned int)pID<=seed->dup_points.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{   				   //Use this to check if clicked on a seed or not		
    float dist =1000.00;
    //float dist1;
    int index;
    float finpt[3];
    for (unsigned int j=0; j<seed->dup_points.size(); j++)

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
	
	seed->dup_points.erase(seed->dup_points.begin()+index);

	//Remove the glyph		
     vtkDataArray* points2del = seed->point1->GetData();    
     //vtkDataArray* points2delred;
     points2del->RemoveTuple((vtkIdType)index);
     seed->point1->SetData(points2del);
     seed->Glyph->SetScaleFactor(seed->Glyph->GetScaleFactor()+0.0001);

     // Add the new blue glyph 
       seed->pointFN->InsertNextPoint(finpt);
	   seed->polydataFN->SetPoints(seed->pointFN);
       seed->FNglyph->SetInput(seed->polydataFN);
       seed->FNsphereMapper->SetInput(seed->FNglyph->GetOutput());	      
       seed->FNglyph->SetScaleFactor(seed->FNglyph->GetScaleFactor()+0.0001);//to rerender immediately
       seed->QVTK->GetRenderWindow()->Render();
       cout<<"Add the new glyph"<<endl;

}
}


// If the Undo Deletebox is checked, mode = 5 ! 
if(seed->mode == 5){

    vtkDataArray* pointIds = seed->delglyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
    int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 

if((unsigned int)pID<=seed->MarkedPoints.size())    
//The ids of non-seed points is much greater than the ids of the seed points 
{    
   float dist =1000.00;
   //float dist1;
   int index;
   float finpt[3];
   for (unsigned int j=0; j<seed->MarkedPoints.size(); j++)
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
		}
}
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

	if(this->mode==1)
    {
	  //double* p1 = this->handle1->GetWorldPosition();
    //double* p2 = this->handle2->GetWorldPosition();       
	  double* p3 = this->handle->GetWorldPosition();       
		//double* bounds = this->Volume->GetBounds();
		float placePoint[3] = {(float)(int)p3[0],(float)(int)p3[1],(float)(int)p3[2]};
        this->point3->InsertNextPoint(placePoint);
    	this->polydata3->SetPoints(this->point3);
        this->addglyph->SetInput(this->polydata3);
        this->addglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);//to rerender immediately
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
	if(this->stateAdd)   { this->mode = 1;}
                      }       
}

void Seed3D::AddSeed()
{
this->mode = 1;
this->stateAdd = this->AddBox->checkState();
if(this->MarkedPoints.size()!=0 && this->stateAdd){
this->msgBox = new QMessageBox(this);
this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes before performing this action?");
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
     this->mode =2;
	this->Apply();
	if(this->stateDelete )
	this->DeleteBox->setCheckState((Qt::CheckState)0); 
	this->mode =1;
	this->AddBox->setCheckState((Qt::CheckState)2);
	break;
	}
   case QMessageBox::No:{
       // Cancel was clicked
       if(this->stateDelete){
	   this->AddBox->setCheckState((Qt::CheckState)0);
	   this->UndoDelBox->setCheckState((Qt::CheckState)0);	
       this->DeleteBox->setCheckState((Qt::CheckState)2);
	   this->mode =2;}
	   break;
                            }
   default:
       // should never be reached
       break;
 }

}


//else{
this->AddBox->setCheckState((Qt::CheckState)this->stateAdd);
if(this->stateAdd){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
        this->mode = 1;
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
////cout<<Id<<endl;
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
}
    this->point2->SetData(points2del);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);    
    this->QVTK->GetRenderWindow()->Render();
    this->MarkedPoints.erase(this->MarkedPoints.begin(),this->MarkedPoints.end());
}

Check();	

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
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   if(this->stateAdd)
   this->mode =1;
   Apply();
   if(this->stateAdd)
   this->AddBox->setCheckState((Qt::CheckState)0);
   this->mode =2;
   this->DeleteBox->setCheckState((Qt::CheckState)2); 
   break;
	}
   case QMessageBox::No:{
  
	   if(this->stateAdd){
	   this->DeleteBox->setCheckState((Qt::CheckState)0);
           this->UndoDelBox->setCheckState((Qt::CheckState)0);	
	   this->AddBox->setCheckState((Qt::CheckState)2);
	   this->mode =1;}  
       break;
                            }
   default:
       // should never be reached
       break;
 }

}

else if(this->MarkedPoints.size()!=0 && this->stateDelete){
}

else{
this->DeleteBox->setCheckState((Qt::CheckState)this->stateDelete);
if(this->stateDelete){
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }////cout<<"hi"<<endl; }
if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
this->mode = 2;
}

else{
this->mode=0;
}
}
}


void Seed3D::saveResult()
{
	QString name = QFileDialog::getSaveFileName(this, tr("Save File As"), this->lastPath, tr("Text Files (*.txt)") );
	if(name == "")
	  {
	  }
	else
	  {
	  QByteArray   bytes  = name.toAscii();
      const char *ptr11    = bytes.data();	
	  std::ofstream seeds;
	  seeds.open(ptr11);
	  
	  if(this->EditRbutton->isChecked()){
			
					for(std::vector<point>::iterator it = dup_points.begin();it!=this->dup_points.end();++it){
							point pt = *it;
							seeds<<pt.x<<" "<<pt.y<<" "<<pt.z<<"\n";		
																											  }
	  }
			
	  if(this->ValidateRbutton->isChecked()){
					seeds<<"True Positives"<<"\n";
					for(std::vector<point>::iterator it = TPVec.begin();it!=this->TPVec.end();++it){
							
							point pt = *it;
							seeds<<pt.x<<" "<<pt.y<<" "<<pt.z<<"\n";		
																											  }
					seeds<<"Cluster Positives"<<"\n";
					for(std::vector<point>::iterator it = CPVec.begin();it!=this->CPVec.end();++it){
					    	
						    point pt = *it;
							seeds<<pt.x<<" "<<pt.y<<" "<<pt.z<<"\n";		
																											  }
					seeds<<"False Positives"<<"\n";
					for(std::vector<point>::iterator it = FPVec.begin();it!=this->FPVec.end();++it){
					    	
						    point pt = *it;
							seeds<<pt.x<<" "<<pt.y<<" "<<pt.z<<"\n";		
																											  }
					seeds<<"False Negatives"<<"\n";
					for(std::vector<point>::iterator it = FNVec.begin();it!=this->FNVec.end();++it){
					    	
						    point pt = *it;
							seeds<<pt.x<<" "<<pt.y<<" "<<pt.z<<"\n";		
																											  }

	  }
	  seeds.close();

		}		
}





void Seed3D::DeleteObjects(){
  //this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->Volume);
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->AddSphereActor);	
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->DelSphereActor);
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->SphereActor);
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->CPSphereActor);
  //this->Volume->Delete();
  this->AddSphereActor->Delete();
  this->DelSphereActor->Delete();	
  this->SphereActor->Delete();

  this->AddSphereMapper->Delete();
  this->DelSphereMapper->Delete();	
  this->SphereMapper->Delete();
  this->CPsphereMapper->Delete();	


  this->polydata1->Delete();
  this->polydata2->Delete();	
  this->polydata3->Delete();
  this->polydataCP->Delete();

 /* this->handle->Delete();
  this->handle1->Delete();
  this->handle2->Delete();
  this->widget->Delete();
  this->widget1->Delete();	
  this->widget2->Delete();*/

  /*this->sliderRep->Delete();
  this->sliderWidget->Delete();
  this->sliderRep2->Delete();	
  this->sliderWidget2->Delete();
  this->sliderRep3->Delete();	
  this->sliderWidget3->Delete();



  this->QVTK->GetRenderWindow()->Render();   
  this->QVTK1->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->imActor);
  this->imActor->Delete();
  this->QVTK1->GetRenderWindow()->Render();   
  this->QVTK2->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->imActor1);
  this->imActor1->Delete();
  this->QVTK2->GetRenderWindow()->Render();   */
  
}


void Seed3D::UndoDeleteSeeds()
{
this->stateUndoDel = this->UndoDelBox->checkState();
this->UndoDelBox->setCheckState((Qt::CheckState)this->stateUndoDel);
if(this->stateUndoDel){
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }
if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
this->UndoDelBox->setCheckState((Qt::CheckState)2);
this->stateUndoDel = this->UndoDelBox->checkState();
this->mode = 5;
}
else{
this->mode=0;
}
}

void Seed3D::Check()
{
//Atleast one marked point required to undo selection for delete
if(this->MarkedPoints.size()==0){
this->UndoDelBox->setCheckState((Qt::CheckState)0);
this->UndoDelBox->setCheckable(false);
}
else{
this->UndoDelBox->setCheckable(true);
}
}




void Seed3D::TruePositives()
{
this->mode = 6;
this->stateTP = this->TPBox->checkState();
this->mode = 6;
std::cout<<this->mode<<"ggg"<<std::endl;
this->TPBox->setCheckState((Qt::CheckState)this->stateTP);
if(this->stateTP){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
    if(this->stateAdd){this->AddBox->setCheckState((Qt::CheckState)0); }
	if(this->stateCP){this->CPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFP){this->FPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFN){this->FNBox->setCheckState((Qt::CheckState)0); }
	this->mode = 6;
				 }
else{
this->mode=0;
}
}


void Seed3D::ClusterPositives()
{
this->mode = 7;
this->stateCP = this->CPBox->checkState();
this->CPBox->setCheckState((Qt::CheckState)this->stateCP);
if(this->stateCP){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
    if(this->stateAdd){this->AddBox->setCheckState((Qt::CheckState)0); }
	if(this->stateTP){this->TPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFP){this->FPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFN){this->FNBox->setCheckState((Qt::CheckState)0); }
	this->mode = 7;
				 }
else{
this->mode=0;
}
}


void Seed3D::FalsePositives()
{
this->mode = 8;
this->stateFP = this->FPBox->checkState();
this->FPBox->setCheckState((Qt::CheckState)this->stateFP);
if(this->stateFP){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
    if(this->stateAdd){this->AddBox->setCheckState((Qt::CheckState)0); }
	if(this->stateCP){this->CPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateTP){this->TPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFN){this->FNBox->setCheckState((Qt::CheckState)0); }
	this->mode = 8;
				 }
else{
this->mode=0;
}
}

void Seed3D::FalseNegatives()
{
this->mode = 9;
this->stateFN = this->FNBox->checkState();
this->FNBox->setCheckState((Qt::CheckState)this->stateFN);
if(this->stateFN){
	if(this->stateDelete){ this->DeleteBox->setCheckState((Qt::CheckState)0);}
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
    if(this->stateAdd){this->AddBox->setCheckState((Qt::CheckState)0); }
	if(this->stateCP){this->CPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateFP){this->FPBox->setCheckState((Qt::CheckState)0); }
	if(this->stateTP){this->TPBox->setCheckState((Qt::CheckState)0); }
	this->mode = 9;
				 }
else{
this->mode=0;
}
}

void Seed3D::ValidateSeeds()
{
vtkIdType Id;
double* pointz;
point p;
Id = this->pointTP->GetNumberOfPoints();
////cout<<Id<<endl;
for(int i =0 ;i<Id;i++)
{
    pointz = this->pointTP->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    this->TPVec.push_back(p);
    
}


Id = this->pointCP->GetNumberOfPoints();
////cout<<Id<<endl;
for(int i =0 ;i<Id;i++)
{
    pointz = this->pointCP->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    this->CPVec.push_back(p);
    
}

Id = this->pointFP->GetNumberOfPoints();

////cout<<Id<<endl;
for(int i =0 ;i<Id;i++)
{
    pointz = this->pointFP->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    this->FPVec.push_back(p);
    
}

Id = this->pointFN->GetNumberOfPoints();

////cout<<Id<<endl;
for(int i =0 ;i<Id;i++)
{
    pointz = this->pointFN->GetPoint(i);
    p.x = pointz[0];
    p.y = pointz[1];
    p.z = pointz[2];	
    this->FNVec.push_back(p);
    
}

}