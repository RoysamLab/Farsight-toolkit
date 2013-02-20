/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


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

	//this->EditRbutton = new QRadioButton("Edit Mode",this->browse);
	//this->ValidateRbutton = new QRadioButton("Validation Mode",this->browse);
	this->AddBtn = new QRadioButton("&Add Seeds", this->browse);
	this->DelBtn = new QRadioButton("&Delete Seeds", this->browse);
	this->UndoDelBtn = new QRadioButton("&Undo Selections for Delete", this->browse);
	DelBtn->setCheckable(false);
	UndoDelBtn->setCheckable(false);

	// Validation checkboxes
	this->CPBtn = new QRadioButton("&Cluster Positives", this->browse);
	this->FPBtn = new QRadioButton("&False Positives", this->browse);
	this->ResetBtn = new QRadioButton("&Reset FP and CP To Original", this->browse);
	ResetBtn->setCheckable(false);

	SetButtonColor(CPBtn, idxCPglyph);
	SetButtonColor(FPBtn, idxFPglyph);
	SetButtonColor(AddBtn, idxAddglyph);
	SetButtonColor(DelBtn, idxdelOrigglyph);

	this->PlaceButton = new QPushButton("&Place Seed", this->browse);
	this->ApplyAddButton = new QPushButton("&Apply Additions", this->browse);  
	this->ApplyDelButton = new QPushButton("&Apply Deletions", this->browse);  

	//********************************************************************************
	// Layouts 

	QGridLayout *buttonLayout = new QGridLayout(this);

	//buttonLayout->addWidget(this->EditRbutton,10,1);
	//buttonLayout->addWidget(this->ValidateRbutton,10,2);
	buttonLayout->addWidget(this->AddBtn,11, 1);
	buttonLayout->addWidget(this->DelBtn,11,2);
	buttonLayout->addWidget(this->UndoDelBtn,11,3);

	//buttonLayout->addWidget(this->TPBox,13,1);
	buttonLayout->addWidget(this->CPBtn,13,2);
	buttonLayout->addWidget(this->FPBtn,13,3);
	buttonLayout->addWidget(this->ResetBtn,13,4);

	buttonLayout->addWidget(this->PlaceButton,12, 1);
	buttonLayout->addWidget(this->ApplyAddButton,12, 2);
	buttonLayout->addWidget(this->ApplyDelButton,12, 3);
	//buttonLayout->addWidget(this->ValidateButton,14, 1);	


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
	this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,isPicked); 
	this->Interactor1->SetInteractorStyle(style1);
	this->Interactor2->SetInteractorStyle(style1);

	//////// make them radio buttons, connect all to one update states func. <<<<====

	//What happens when the user clicks the buttons or the Checkboxes is determined by the "SLOTS"
	connect(this->PlaceButton, SIGNAL(clicked()), this, SLOT(PlaceSeed()));
	//connect(this->ValidateButton, SIGNAL(clicked()), this, SLOT(ValidateSeeds()));
	//connect(this->AddBtn, SIGNAL(stateChanged(int)), this, SLOT(StateChange(AddBox)));
	connect(this->AddBtn, SIGNAL(released()), this, SLOT(AddSeedCheck()));

	connect(this->ApplyAddButton, SIGNAL(clicked()), this, SLOT(ApplyAdd()));
	connect(this->ApplyDelButton, SIGNAL(clicked()), this, SLOT(ApplyDel()));		

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

	loadActionGT = new QAction(tr("Load Seed File to Compare"), this);
	loadActionGT->setStatusTip(tr("Load another set of seeds such as ground truth for comparison..."));
	connect(loadActionGT, SIGNAL(triggered()), this, SLOT(loadSeedsGT()));
	fileMenu->addAction(loadActionGT); 
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


void Seed3D::loadSeedsGT()
{
	GTglyph->ReadFilePts ( SeedFileDialog("Ground Truth") );
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


	// if image already exists delete the objects

	if(this->iRender==1){
		CleanUp();
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
	Volume = vtkVolume::New();
	Volume->SetMapper(volumeMapper);
	Volume->SetProperty(volumeProperty);
	Volume->SetPickable(0);
	//this->Volume = volume;


	// Create seeds (spheres) only at coordinates specified by the text file.
	// I also create Red (Markedpoints) and Green Spheres (Markedpoints2add) initially and delete them.
	//  So that the first delete and addition by the user would be 
	// "inserts" into the vector !
	// Note: 1) All the seeds form a single glyph.
	//       2) The green spheres form a single glyph. These are seeds added but not applied permanently.
	//		 3) The red spheres form a single glyph.These seeds are marked for Deletion but can be unmarked

	//std::vector<point> spPoint;
	std::vector<point> MarkedPoints;
	std::vector<point> MarkedPoints2add;



	// create the glyphing the sphere with a cone.  Create the mapper
	// and actor for the glyphs.

	vtkSmartPointer<vtkSphereSource> sphere = vtkSphereSource::New();

	this->imActor = imActor;
	Renderer->AddVolume(Volume);

	// create myGlyphs with corresponding colors
	for (int i=0; i<GLYPHMAX_IDX; i++)
		myGlyphArray[i] = new myGlyph(this, colorArray[i], sphere);
	// set their alias names for convenience
	delAddedglyph = myGlyphArray[idxdelAddedglyph]; 
	delOrigglyph  = myGlyphArray[idxdelOrigglyph];
	Addglyph      = myGlyphArray[idxAddglyph];
	Origglyph     = myGlyphArray[idxOrigglyph];
	CPglyph       = myGlyphArray[idxCPglyph];
	FPglyph       = myGlyphArray[idxFPglyph];
	GTglyph       = myGlyphArray[idxGTglyph];

	Origglyph->ReadFilePts( SeedFileDialog("test" ) );
	DelBtn->setCheckable(true);
	// The active camera for the main window!
	vtkCamera *cam1 = Renderer->GetActiveCamera();
	cam1->SetViewUp (0, 1, 0);
	cam1->SetPosition (0, 0, 50);
	Renderer->ResetCamera();	
	//this->Volume = volume;

	//********************************************************************************
	// This section of code deals with the creation of the handles and handle widget

	// Create the widget and its representation
	vtkSmartPointer<vtkPointHandleRepresentation3D> handle = vtkPointHandleRepresentation3D::New();
	vtkSmartPointer<vtkHandleWidget> widget = vtkHandleWidget::New();
	widget->SetInteractor(this->Interactor);
	widget->SetRepresentation(handle);
	//this->widget = widget;
	this->handle = handle; 
	vtkSmartPointer<vtkSeedCallback> scbk = vtkSeedCallback::New();
	widget->AddObserver(vtkCommand::EndInteractionEvent,scbk);
	widget->SetEnabled(true);

	//********************************************************************************
	// Create the Imagereslice planes and the actors to the respectove Qwidgets

	static double axialElements[16] = {
		1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 };

		vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkMatrix4x4::New();
		resliceAxes->DeepCopy(axialElements);
		vtkSmartPointer<vtkImageReslice> reslicexy = vtkImageReslice::New();
		reslicexy->SetInput(vtkim);
		reslicexy->SetOutputDimensionality(2);
		reslicexy->SetResliceAxes(resliceAxes);
		reslicexy->SetInterpolationModeToLinear();

		// Set the extent for the Y-Z slice.
		double* origincalc = this->Volume->GetBounds();	
		int origin[6];
		origin[0] = (int)origincalc[2];
		origin[1] = (int)origincalc[3];
		origin[2] = (int)origincalc[4];
		origin[3] = (int)origincalc[5];
		origin[4] = (int)origincalc[0];
		origin[5] = (int)origincalc[1];

		reslicexy->SetOutputExtent(origin); 	 
		reslicexy->SetResliceAxesOrigin((origincalc[1]-origincalc[0])/2.0,(origincalc[3]-origincalc[2])/2.0,(origincalc[5]-origincalc[4])/2.0);
		//this->Reslice = reslice;

		vtkDataSetMapper *im_mapper = vtkDataSetMapper::New();
		im_mapper->SetInput(reslicexy->GetOutput());
		vtkActor *imActor = vtkActor::New();
		imActor->SetMapper(im_mapper);

		this->Renderer1->AddActor(imActor);
		this->QVTK1->GetRenderWindow()->Render();
		this->Renderer1->ResetCamera();	
		this->imActor = imActor;
		this->imActor->RotateWXYZ(180,0,0,1);

		static double sagittalElements[16] = { 0, 0,-1, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1 }; 
		vtkSmartPointer<vtkMatrix4x4> resliceAxes1 = vtkMatrix4x4::New();
		resliceAxes->DeepCopy(sagittalElements);

		vtkSmartPointer<vtkImageReslice> resliceyz = vtkImageReslice::New();
		resliceyz->SetInput(vtkim);
		resliceyz->SetOutputDimensionality(2);
		resliceyz->SetResliceAxes(resliceAxes1);
		resliceyz->SetInterpolationModeToLinear();

		origincalc = this->Volume->GetBounds();
		resliceyz->SetResliceAxesOrigin((origincalc[1]-origincalc[0])/2.0,(origincalc[3]-origincalc[2])/2.0,(origincalc[5]-origincalc[4])/2.0);
		//this->resliceyz = resliceyz;

		vtkDataSetMapper *im_mapper1 = vtkDataSetMapper::New();
		im_mapper1->SetInput(resliceyz->GetOutput());
		vtkActor *imActor1 = vtkActor::New();
		imActor1->SetMapper(im_mapper1);

		this->Renderer2->AddActor(imActor1);
		this->QVTK2->GetRenderWindow()->Render();
		this->Renderer2->ResetCamera();
		this->imActor1 = imActor1;

		//********************************************************************************  
		// handle  -> handle for the main window
		// handle1 -> handle for the Y-Z plane actor (imActor formed by reslice)
		// handle2 -> handle for the X-Y plane actor (imActor formed by resliceyz)

		  scbk->handle = this->handle;
  
		scbk->reslice  = reslicexy;
		scbk->reslice1 = resliceyz;
		scbk->imActor1 = this->imActor1; 
		scbk->imActor =  this->imActor;
		scbk->QVTK1 = this->QVTK1;
		scbk->QVTK2 = this->QVTK2;
		scbk->vol = this->Volume;  

		vtkSmartPointer<vtkPointHandleRepresentation2D> handle1 = vtkPointHandleRepresentation2D::New();
		vtkSmartPointer<vtkHandleWidget> widget1 = vtkHandleWidget::New();
		widget1->SetInteractor(this->Interactor1);
		widget1->SetRepresentation(handle1);
		//this->widget1 = widget1;
		//this->handle1 = handle1; 

		vtkSeedCallback1 *scbk1 = vtkSeedCallback1::New();
		widget1->AddObserver(vtkCommand::EndInteractionEvent,scbk1);
		widget1->SetEnabled(true);

		vtkSmartPointer<vtkPointHandleRepresentation2D> handle2 = vtkPointHandleRepresentation2D::New();
		vtkSmartPointer<vtkHandleWidget> widget2 = vtkHandleWidget::New();
		widget2->SetInteractor(this->Interactor2);
		widget2->SetRepresentation(handle2);
		//this->widget2 = widget2;
		//this->handle2 = handle2; 

		vtkSeedCallback2 *scbk2 = vtkSeedCallback2::New();
		widget2->AddObserver(vtkCommand::EndInteractionEvent,scbk2);
		widget2->SetEnabled(true);
		scbk->handle1 = handle1;
		scbk->handle2 = handle2;
		scbk1->handle1= handle1;	
		scbk1->handle2= handle2;	
		scbk1->QVTK2 = this->QVTK2;
		scbk1->QVTK1 = this->QVTK1;
		scbk1->QVTK = this->QVTK;	
		scbk1->vol = this->Volume;
		scbk1->handle = handle; 
		//scbk1->reslice = this->Reslice;
		scbk1->reslice1 = resliceyz;
		scbk1->imActor1 = this->imActor1; 
		scbk2->handle1= handle1;	
		scbk2->handle2= handle2;
		scbk2->QVTK2 = this->QVTK2;	
		scbk2->QVTK1 = this->QVTK1;
		scbk2->QVTK = this->QVTK;	
		scbk2->vol = this->Volume;
		scbk2->handle = handle; 
		scbk2->reslice = reslicexy;
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
		callback_seedsize->Origglyph     = 	Origglyph->Glyph;
		callback_seedsize->Addglyph      = Addglyph->Glyph;
		callback_seedsize->delOrigglyph  = delOrigglyph->Glyph;	
		callback_seedsize->delAddedglyph = delAddedglyph->Glyph;
		callback_seedsize->CPglyph       = CPglyph->Glyph;
		callback_seedsize->FPglyph       = FPglyph->Glyph;
		callback_seedsize->GTglyph       = GTglyph->Glyph;


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


	if(seed->DelBtn->isChecked() )
	{
		// delete: search closest seed whether original or added
		// and transfer it to the corresponding delete glyph
		seed->PointTransfer(pickPos, seed->Origglyph, seed->Addglyph, seed->delOrigglyph, seed->delAddedglyph);
		if (seed->Origglyph->isEmpty() && seed->Addglyph->isEmpty())
			seed->DelBtn->setCheckable(false);
		if (!seed->delAddedglyph->isEmpty() || !seed->delOrigglyph->isEmpty())
			seed->UndoDelBtn->setCheckable(true);
	}

	

	if(seed->ResetBtn->isChecked()) 
	{
		// delete: search closest seed whether original or added
		// and transfer it to the corresponding delete glyph
		seed->PointTransfer(pickPos, seed->FPglyph, seed->CPglyph, seed->Origglyph, NULL);
		if (seed->FPglyph->isEmpty() && seed->CPglyph->isEmpty() )
			seed->ResetBtn->setCheckable(false);
		if (!seed->Origglyph->isEmpty()  )
			seed->DelBtn->setCheckable(true);
	}
		

	if(seed->CPBtn->isChecked()) 
	{
		// add to CP
		seed->PointTransfer(pickPos, seed->Origglyph, NULL, seed->CPglyph, NULL);
		if ( seed->Origglyph->isEmpty() )
			seed->DelBtn->setCheckable(false);
		if (!seed->CPglyph->isEmpty() )
			seed->ResetBtn->setCheckable(true);
	}


	if(seed->FPBtn->isChecked()) {
		// add to FP
		seed->PointTransfer(pickPos, seed->Origglyph, NULL, seed->FPglyph, NULL);
		if ( seed->Origglyph->isEmpty()  )
			seed->DelBtn->setCheckable(false);
		if (!seed->FPglyph->isEmpty() )
			seed->ResetBtn->setCheckable(true);
	}



	// If the Undo Deletebox is checked, mode = 5 ! 
	if(seed->UndoDelBtn->isChecked()){
		//undo delete
		seed->PointTransfer(pickPos, seed->delAddedglyph, seed->delOrigglyph, seed->Addglyph, seed->Origglyph);
		if (seed->delAddedglyph->isEmpty() && seed->delOrigglyph->isEmpty())
			seed->UndoDelBtn->setCheckable(false);
		if ( !seed->Origglyph->isEmpty() || !seed->Addglyph->isEmpty()  )
			seed->DelBtn->setCheckable(true);
	}
	//seed->Check();
}



void Seed3D::PlaceSeed()
{

	if(this->AddBtn->isChecked())
	{
		double* p3 = this->handle->GetWorldPosition();       

		float placePoint[3] = {(float)(int)p3[0],(float)(int)p3[1],(float)(int)p3[2]};
		Addglyph->AddPoint(placePoint);
		this->QVTK->GetRenderWindow()->Render();
		DelBtn->setCheckable(true);
	}       
}

void Seed3D::AddSeedCheck()
{
	if (!AddBtn->isChecked())
		return;

	if( !delOrigglyph->isEmpty() || !delAddedglyph->isEmpty() ) 
	{
		QMessageBox* msgBox;

		msgBox = new QMessageBox(this);
		msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes before performing this action?");
		msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
		msgBox->setDefaultButton(QMessageBox::Yes);
		int ret = msgBox->exec();

		delete msgBox;

		switch (ret) 
		{
		case QMessageBox::Yes:
			{
				ApplyDel();
				break;
			}
		default:
			break;
		}

	}

}

void myGlyph::Clear()
{
	int numpts = GetNumberOfPoints();
	for(int i =0; i<numpts; i++)
		points->GetData()->RemoveFirstTuple();
	OwnerClass->QVTK->GetRenderWindow()->Render();
	Glyph->SetScaleFactor(Glyph->GetScaleFactor()+0.0001);
}

void Seed3D::ApplyAdd()
{
	int numpts;
	numpts = Addglyph->GetNumberOfPoints();
	for(int i =0 ; i<numpts;i++)
	{
			float pt[3];
			Addglyph->GetPoint(i, pt);

			Addglyph->RemovePoint(i);
			Origglyph->AddPoint(pt);
	}
	this->QVTK->GetRenderWindow()->Render();
	//Check();	

}

void Seed3D::ApplyDel()
{
	delOrigglyph->Clear();
	delAddedglyph->Clear();
	UndoDelBtn->setCheckable(false);
}



void Seed3D::saveResult()
{

	QString name = QFileDialog::getSaveFileName(this, tr("Save File As"), this->lastPath, tr("Text Files (*.txt)") );
	if(!name.isEmpty())
	{
		QByteArray   bytes  = name.toAscii();
		const char *ptr11    = bytes.data();	
		std::ofstream *outfilestrm;
		outfilestrm = new std::ofstream;
		outfilestrm->open(ptr11);
		//seeds<<"Added Seeds (false negatives) " << newseedcounter << "\n";  
		for (int i=0; i<GLYPHMAX_IDX; i++)
		{
			if (!myGlyphArray[i]->isEmpty())
			{
				*outfilestrm << GLYPH_NAMES[i];
				myGlyphArray[i]->SaveToOutfile(	outfilestrm );
			}
		}
		outfilestrm->close();
	}
}





//was DeleteObjects()
void Seed3D::CleanUp(){
	this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->Volume);
	
	delete CPglyph;
	delete FPglyph;
	delete Origglyph;
	delete delAddedglyph;
	delete delOrigglyph;
	delete Addglyph;
	delete GTglyph;
	//this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->AddSphereActor);	
	//this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->DelSphereActor);
	//this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->SphereActor);
	//this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->CPSphereActor);
	////this->Volume->Delete();
	//this->AddSphereActor->Delete();
	//this->DelSphereActor->Delete();	
	//this->SphereActor->Delete();

	//this->AddSphereMapper->Delete();
	//this->DelSphereMapper->Delete();	
	//this->SphereMapper->Delete();
	//this->CPsphereMapper->Delete();	


	//this->polydata1->Delete();
	//this->polydata2->Delete();	
	//this->polydata3->Delete();
	//this->polydataCP->Delete();

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

myGlyph::~myGlyph()
{
	OwnerClass->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(SphereActor);
}

myGlyph::myGlyph(Seed3D *Owner, const float *color, vtkSmartPointer<vtkSphereSource> sphere)
{
	//float pts[3] = {0.0,0.0,0.0};
	OwnerClass = Owner;
	pcoords = vtkFloatArray::New();
	pcoords->SetNumberOfComponents(3);
	pcoords->SetNumberOfTuples(0);
	//pcoords->SetTuple(0,pts);	

	points = vtkPoints::New();
	points->SetData(pcoords);	  
	// Assign points and cells
	polydata = vtkPolyData::New();
	polydata->SetPoints(points);
	Glyph = vtkGlyph3D::New();
	Glyph->SetInput(polydata);

	Glyph->SetSource(sphere->GetOutput());

	Glyph->SetVectorModeToUseNormal();

	Glyph->SetScaleModeToScaleByVector();

	Glyph->SetScaleFactor(1);

	Glyph->GeneratePointIdsOn();

	SphereMapper = vtkPolyDataMapper::New();
	SphereMapper->SetInput(Glyph->GetOutput());

	
	SphereActor = vtkLODActor::New();

	SphereActor->SetMapper(SphereMapper);
	SphereActor->GetProperty()->SetColor(color[0], color[1], color[2]);
	Owner->Renderer->AddActor(SphereActor);

	/*vtkSmartPointer<vtkDataArray> tmppoint = points->GetData();    
	tmppoint->RemoveTuple((vtkIdType)0);
	points->SetData(tmppoint);*/
	Glyph->SetScaleFactor(Glyph->GetScaleFactor()+0.0001);

}

std::string  Seed3D::SeedFileDialog(const char *seedtype)
{
 	char msg[50];
	sprintf(msg, "Select %s seed file to open", seedtype);
	QString fileNameSeed = QFileDialog::getOpenFileName(
		this, msg, lastPath,
		tr("Text files (*.txt);; \n"
		"All Files (*.*)"));
	/*if(fileNameSeed == "")
		return ;*/
	//QString name = QFileInfo(this->fileNameSeed).baseName();	
	std::string sfilename = fileNameSeed.toStdString();
	return sfilename;
}

void myGlyph::ReadFilePts(std::string sfilename)
{
	//char* seedname;
	//seedname = &sfilename[0];
	if (!sfilename.empty() && sfilename.c_str())
	{
		ifstream inputstream(sfilename.c_str());

		if (inputstream.fail())
		{
			cout << "ERROR: Incorrect File Type" << endl;
			//return 1;
		}
		else
		{
			cout << "Reading seed points file ..." << endl;
		}
		// Create a sphere only at these particular coordinates only.


		float p[3];
		int coordx=0;
		//char comments[200];
		char cchk;
		while(!inputstream.eof())
		{
			cchk = inputstream.peek();
			if (cchk>='0' && cchk<='9')
			{
				inputstream >> p[0] >> p[1] >> p[2];
				//std::cout << "reading pt" << p << std::endl;
				AddPoint(p);
				coordx++;	
			}
			else
			{
				inputstream.ignore(2000, '\n');
			}
		}
		std::cout << "read " <<  coordx << " pts"  << std::endl;
		// not sure if this is needed now
		// pcoords->SetNumberOfTuples(coordx);
	}
	else
		std::cout << "No seed file to read" << std::endl;
}	


vtkIdType myGlyph::Search(double *pickPos, float &mindist, float *finpt)
{

	vtkIdType index=-1;
	//if ( pID<=points->GetNumberOfPoints()  )
		//The ids of non-seed points is much greater than the ids of the seed points 
	//{    
		float dist = MAX_SEARCH_DIST;
		//float dist1;
		//float finpt[3];
		for ( int j=0; j<points->GetNumberOfPoints(); j++)
		{

			float p1[3];
			pcoords->GetTupleValue(j, p1);
			//= {seed->MarkedPoints[j].x, seed->MarkedPoints[j].y ,seed->MarkedPoints[j].z };	
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
		mindist = dist;

	//}
	return index;
}

void Seed3D::PointTransfer(double *pickPos, myGlyph *src1, myGlyph *src2, myGlyph *target1, myGlyph *target2)
{
	float dist1 = MAX_SEARCH_DIST, dist2 = MAX_SEARCH_DIST;
	float delpt1[3], delpt2[3];
	vtkIdType delidx1 =-1, delidx2=-1;
	delidx1 = src1->Search(pickPos, dist1, delpt1);
	myGlyph *t2;

	if (src2)
		delidx2 = src2->Search(pickPos, dist2, delpt2);

	if (delidx1 == -1 && delidx2 == -1)
		return;

	if (dist1 < dist2) 
	{
		src1->RemovePoint(delidx1);
		target1->AddPoint(delpt1);
		//target1->points->InsertNextPoint(delpt1);
	}
	else
		if (src2)
		{
			src2->RemovePoint(delidx2);
			if (target2)
				t2=target2;
			else
				t2=target1;
			t2->AddPoint(delpt2);
			//t2->points->InsertNextPoint(delpt2);
		}

	QVTK->GetRenderWindow()->Render();		
}

//void myGlyph::AddPoint(vtkIdType idx, float *p)
//{
//	pcoords->SetTuple(idx, p);
//	Glyph->SetScaleFactor(Glyph->GetScaleFactor()+0.0001);
//	Glyph->Update();
//}

void myGlyph::AddPoint(float *p)
{
	points->InsertNextPoint(p);
	Glyph->SetScaleFactor(Glyph->GetScaleFactor()+0.0001);
	//Glyph->Update();
}

void myGlyph::RemovePoint(vtkIdType idx)
{
	points->GetData()->RemoveTuple(idx);
	//int numpts=points->GetNumberOfPoints();
	Glyph->SetScaleFactor(Glyph->GetScaleFactor()+0.0001);
	//Glyph->Update();
}

void myGlyph::SaveToOutfile(	std::ofstream *outfstrm )
{
	int num = points->GetNumberOfPoints();
 	*outfstrm << "   --- Total= " << num << "\n";
	for( int i=0; i<num; i++)
	{
		float pt[3];
		pcoords->GetTupleValue(i, pt);
		*outfstrm<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<"\n";		
	}
}

void Seed3D::SetButtonColor(QRadioButton *Btn, GlyphIndex idx)
{
	int rgbcol[3];
	char rgbstr[3][4];
	char style[100];
	for (int i=0; i<3; i++)
	{
		rgbcol[i] = (int)255*colorArray[idx][i];
		sprintf(rgbstr[i], "%d", rgbcol[i]);
	}
	sprintf(style, "* {background-color:rgb(%s,%s,%s); padding: 7px ; color:rgb(255,255,255)}", rgbstr[0],rgbstr[1],rgbstr[2]);
	Btn->setStyleSheet(style);
}
