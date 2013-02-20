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

//***********************************************************************************
// TableWindow provides a classic table view into a model
//***********************************************************************************
#include "mdlGUI.h"

//-----------------------------------------------------------------------------
mdlGUI::~mdlGUI()
{
  //delete this->OutputWindow;
  //this->OutputWindow = 0;
  //delete this->HelpWindow;
  //this->HelpWindow = 0;
	delete this->RenderWidget;
	this->RenderWidget = 0;
  //this->Helper->exit();
  //delete this->Helper;
  //this->Helper = 0;

	if(this->Skel)
    {
		delete this->Skel;
		this->Skel = 0;
    }
	if(this->MinSpanTree)
    {
		delete this->MinSpanTree;
		this->MinSpanTree = 0;
    }
}

//Constructors:
mdlGUI::mdlGUI(QWidget * parent)
: QMainWindow(parent)
{

	bannerLabel = new QLabel();
	bannerPixmap = QPixmap(":/images/banner.png");
	bannerLabel->setPixmap(bannerPixmap);
	bannerLabel->setFixedWidth(bannerPixmap.width());
	bannerLabel->setScaledContents(true);

	loadLayout = new QHBoxLayout();
	loadLabel = new QLabel(tr("Input Image: "));
	loadLayout->addWidget(loadLabel);
	imageLabel = new QLabel();
	loadLayout->addWidget(imageLabel);
	loadButton = new QPushButton(tr("Load..."));
	connect(loadButton, SIGNAL(clicked()), this, SLOT(loadImage()));
	loadLayout->addWidget(loadButton);

	prepLayout = new QHBoxLayout();
	prepButton = new QPushButton(tr("PREPROCESS"));
	prepButton->setEnabled(false);
	connect(prepButton, SIGNAL(clicked()), this, SLOT(preprocess()));
	prepLayout->addWidget(prepButton);

	skelBox = new QGroupBox(tr("Skeletonization"));
	skelLayout = new QHBoxLayout(skelBox);
	vectMagLabel = new QLabel(tr("Vector Magnitude: "));
	skelLayout->addWidget(vectMagLabel);
	vectMagEdit = new QLineEdit();
	vectMagEdit->setText("0.04");
	vectMagEdit->setMaximumWidth(30);
	skelLayout->addWidget(vectMagEdit);
	skelLayout->addStretch(10);
	skelCheckbox = new QCheckBox(tr("Use ITK"));
	skelCheckbox->setChecked(false);
	skelLayout->addWidget(skelCheckbox);
	skelButton = new QPushButton(tr("Skeletonize"));
	skelButton->setEnabled(false);
	connect(skelButton, SIGNAL(clicked()), this, SLOT(integratedSkeleton()));
	skelLayout->addWidget(skelButton);

	bbBox = new QGroupBox(tr("Minimum Spanning Tree"));
	bbLayout = new QHBoxLayout(bbBox);
	edgeLabel = new QLabel(tr("Edge Range: "));
	bbLayout->addWidget(edgeLabel);
	edgeRangeEdit = new QLineEdit();
	edgeRangeEdit->setText("5");
	edgeRangeEdit->setMaximumWidth(30);
	bbLayout->addWidget(edgeRangeEdit);
	morphLabel = new QLabel(tr("Prune Iterations: "));
	bbLayout->addWidget(morphLabel);
	morphStrengthEdit = new QLineEdit();
	morphStrengthEdit->setText("10");
	morphStrengthEdit->setMaximumWidth(30);
	bbLayout->addWidget(morphStrengthEdit);
	bbLayout->addStretch(10);
	bbButton = new QPushButton(tr("Backbone Extract"));
	bbButton->setEnabled(false);
	connect(bbButton, SIGNAL(clicked()), this, SLOT(mstBB()));
	bbLayout->addWidget(bbButton);

	spineBox = new QGroupBox(tr("Spines"));
	spineLayout = new QHBoxLayout(spineBox);
	pruneLabel = new QLabel(tr("Prune Threshold: "));
	spineLayout->addWidget(pruneLabel);
	pruneEdit = new QLineEdit();
	pruneEdit->setText("4");
	pruneEdit->setMaximumWidth(30);
	spineLayout->addWidget(pruneEdit);
	alphaLabel = new QLabel(tr("Alpha: "));
	spineLayout->addWidget(alphaLabel);
	alphaEdit = new QLineEdit();
	alphaEdit->setText("0.5");
	alphaEdit->setMaximumWidth(30);
	spineLayout->addWidget(alphaEdit);
	spineLayout->addStretch(10);
	spineButton = new QPushButton(tr("Spine Extract"));
	spineButton->setEnabled(false);
	connect(spineButton, SIGNAL(clicked()), this, SLOT(mstSpine()));
	spineLayout->addWidget(spineButton);

	saveLayout = new QHBoxLayout();
	saveButton = new QPushButton(tr("Save Result..."));
	saveButton->setEnabled(false);
	connect(saveButton, SIGNAL(clicked()), this, SLOT(saveAll()));
	saveLayout->addWidget(saveButton);

	showRenderLayout = new QHBoxLayout();
	showRenderButton = new QPushButton(tr("Show Render Window"));
	connect(showRenderButton, SIGNAL(clicked()), this, SLOT(showRenderWindow()));
	showRenderLayout->addStretch(10);
	showRenderLayout->addWidget(showRenderButton);

	stepsLayout = new QVBoxLayout();
	stepsLayout->addLayout(loadLayout);
	stepsLayout->addLayout(prepLayout);
	stepsLayout->addWidget(skelBox);
	stepsLayout->addWidget(bbBox);
	stepsLayout->addWidget(spineBox);
	stepsLayout->addLayout(saveLayout);
	stepsLayout->addStretch(10);
	stepsLayout->addLayout(showRenderLayout);

	masterLayout = new QHBoxLayout();
	masterLayout->addWidget(bannerLabel);
	masterLayout->addLayout(stepsLayout);

	masterWidget = new QWidget();
	masterWidget->setLayout(masterLayout);
	this->setCentralWidget(masterWidget);

	this->setWindowTitle(tr("MDL GUI"));
	//this->move(10,10);
	
	//set up the render widget: a QVTKWidget to show intermediary results
	this->RenderWidget = new QVTKWidget();
	this->Renderer = vtkSmartPointer<vtkRenderer>::New();
	this->Renderer->SetBackground(0.33, 0.35, 0.43);
	this->RenderWidget->GetRenderWindow()->AddRenderer(this->Renderer);
	this->RenderWidget->resize(500, 500);
	//this->RenderWidget->move( this->frameGeometry().bottomLeft() );
	this->RenderWidget->setWindowTitle("Render Window");
	this->RenderWidget->show();

	OrigImage = NULL;
	this->ITKtoVTK = ITKtoVTKType::New();
	PrepImage = NULL;
	this->Skel = 0;
	this->MinSpanTree = 0;

	lastPath = "";

}

void mdlGUI::loadImage()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Select Input Image"), lastPath, tr("Image Files (*.*)"));
	if(filename == "")
		return;

	lastPath = QFileInfo(filename).absolutePath() + QDir::separator();

	//Put the filename in the label:
	this->imageLabel->setText(filename);

	typedef itk::ImageFileReader< mdl::ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename.toStdString());
	try
	{
		reader->Update();
    }
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "READER FAILED: " << err << std::endl;
		//emit a signal about failing here
		return;
    }

	this->OrigImage = reader->GetOutput();

	//Render the loaded image:
	this->Renderer->RemoveAllViewProps(); //Remove previous image
	this->RenderImage( this->OrigImage );
	this->Renderer->ResetCamera();
	this->RenderWidget->setWindowTitle("Input Image");
	this->RenderWidget->GetRenderWindow()->Render();

	this->prepButton->setEnabled(true);
	this->prepButton->setFocus();
	this->skelButton->setEnabled(false);
	this->bbButton->setEnabled(false);
	this->saveButton->setEnabled(false);
	this->spineButton->setEnabled(false);
}

void mdlGUI::preprocess()
{
	PreprocessDialog * dialog = new PreprocessDialog(lastPath, this);
	dialog->SetImage( this->OrigImage );

	if(!dialog->exec())
	{
		delete dialog;
		return;
	}

	this->PrepImage = dialog->GetImage();
	delete dialog;

	this->Renderer->RemoveAllViewProps(); //Remove previous image
	this->RenderImage( this->PrepImage );
	this->RenderWidget->setWindowTitle("Preprocessed Image");
	this->RenderWidget->GetRenderWindow()->Render();

	this->skelButton->setEnabled(true);
	this->skelButton->setFocus();
	this->bbButton->setEnabled(false);
	this->saveButton->setEnabled(false);
	this->spineButton->setEnabled(false);
}

void mdlGUI::integratedSkeleton()
{
	double vectorMagnitude = this->vectMagEdit->text().toDouble();

	if(this->skelCheckbox->isChecked())
	{
		ftk::Preprocess *pre = new ftk::Preprocess( this->PrepImage );
		pre->BinaryThinning();
		pre->SaveVTKPoints( lastPath.toStdString() + "out.seed" );
		delete pre;

		mdl::vtkFileHandler * FileHandler = new mdl::vtkFileHandler();
		FileHandler->SetNodes( &(this->SkeletonPoints) );
		FileHandler->Read( lastPath.toStdString() + "out.seed" );
		delete FileHandler;
	}
	else
	{
		if(this->Skel != 0)
		{
			delete this->Skel;
		}
		this->Skel = new mdl::IntegratedSkeleton( this->PrepImage );
		this->Skel->SetDebug(false);
		this->Skel->SetUseXiaoLiangMethod(false);
		this->Skel->SetVectorMagnitude(vectorMagnitude);
		this->Skel->Update();
		this->SkeletonPoints = this->Skel->GetOutput();
	}

	//Render the points:
	this->Renderer->RemoveAllViewProps(); //Remove previous image
	this->RenderImage( this->PrepImage );
	this->RenderPoints( this->SkeletonPoints, 1.0, 1.0, 0.0 );
	this->RenderWidget->setWindowTitle("Skeleton Points");
	this->RenderWidget->GetRenderWindow()->Render();

	this->bbButton->setEnabled(true);
	this->bbButton->setFocus();
	this->saveButton->setEnabled(false);
	this->spineButton->setEnabled(false);
}

void mdlGUI::mstBB()
{
	int edgeRange = this->edgeRangeEdit->text().toInt();
	int morphStrength = this->morphStrengthEdit->text().toInt();

	if(this->MinSpanTree != 0)
    {
		delete this->MinSpanTree;
    }
	this->MinSpanTree = new mdl::MST( this->PrepImage );
	this->MinSpanTree->SetDebug(false);
	this->MinSpanTree->SetUseVoxelRounding(false);
	this->MinSpanTree->SetEdgeRange(edgeRange);
	this->MinSpanTree->SetPower(1);
	this->MinSpanTree->SetSkeletonPoints( &(this->SkeletonPoints) );
	this->MinSpanTree->CreateGraphAndMST();
	this->MinSpanTree->ErodeAndDialateNodeDegree(morphStrength);
	this->Nodes = this->MinSpanTree->GetNodes();
	this->BackbonePairs = this->MinSpanTree->BackboneExtract();

	mdl::vtkFileHandler * FileHandler = new mdl::vtkFileHandler();
	FileHandler->SetNodes( &(this->Nodes) );
	FileHandler->SetLines( &(this->BackbonePairs) );
	FileHandler->Write( lastPath.toStdString() + "BackboneCandidate.vtk" );
	delete FileHandler;

	//Render the Backbone:
	this->Renderer->RemoveAllViewProps(); //Remove previous image
	this->RenderImage( this->PrepImage );
	this->RenderPolyData( lastPath.toStdString() + "BackboneCandidate.vtk", 1.0, 1.0, 1.0 );
	this->RenderWidget->setWindowTitle("Backbone");
	this->RenderWidget->GetRenderWindow()->Render();

	this->saveButton->setEnabled(true);
	this->saveButton->setFocus();
	this->spineButton->setEnabled(true);
}

void mdlGUI::mstSpine()
{
	if(!this->MinSpanTree)
		return;

	//double pruneThreshold = this->pruneEdit->text().toDouble();
	//double alpha = this->alphaEdit->text().toDouble();

	this->SpinePairs = this->MinSpanTree->SpineExtract();

	mdl::vtkFileHandler * FileHandler = new mdl::vtkFileHandler();
	FileHandler->SetNodes( &(this->Nodes) );
	FileHandler->SetLines( &(this->SpinePairs) );
	FileHandler->Write( lastPath.toStdString() + "SpineCandidate.vtk" );
	delete FileHandler;

	//Render the Backbone:
	this->Renderer->RemoveAllViewProps(); //Remove previous image
	this->RenderImage( this->PrepImage );
	this->RenderPolyData( lastPath.toStdString() + "BackboneCandidate.vtk", 1.0, 1.0, 1.0 );
	this->RenderPolyData( lastPath.toStdString() + "SpineCandidate.vtk", 0.0, 1.0, 0.0 );
	this->RenderWidget->setWindowTitle("Backbone + Spine");
	this->RenderWidget->GetRenderWindow()->Render();

	this->saveButton->setFocus();

}

void mdlGUI::saveAll()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Results As"), lastPath, tr("VTK Files (*.vtk*)"));
	if(filename == "")
		return;

	lastPath = QFileInfo(filename).absolutePath() + QDir::separator();

	//Combine both Backbone and Spine:
	std::vector<mdl::pairE> allPairs = this->BackbonePairs;
	allPairs.insert( allPairs.end(), SpinePairs.begin(), SpinePairs.end() );

	mdl::vtkFileHandler * FileHandler = new mdl::vtkFileHandler();
	FileHandler->SetNodes( &(this->Nodes) );
	FileHandler->SetLines( &(allPairs) );
	FileHandler->Write(filename.toStdString());
	delete FileHandler;
}

void mdlGUI::showRenderWindow()
{
	this->RenderWidget->show();
}

//-----------------------------------------------------------------------------
void mdlGUI::RenderImage(mdl::ImageType::Pointer image)
{
	this->ITKtoVTK->SetInput( image );
	this->ITKtoVTK->Update();
	vtkSmartPointer<vtkImageData> vtk_image = this->ITKtoVTK->GetOutput();

	//remove all old actors
	//this->Renderer->RemoveAllViewProps();
	//change the window title so the user knows what file is being displayed
	//this->RenderWidget->setWindowTitle(windowName);

	//create a volume from the specified image
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> mapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	mapper->SetInput( vtk_image );

	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(mapper);

	double range[2];
	vtk_image->GetScalarRange(range);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = this->NewRGBVolumeProperty(range);
	volume->SetProperty(volumeProperty);

	this->Renderer->AddVolume(volume);
	//this->Renderer->ResetCamera();
	//this->RenderWidget->GetRenderWindow()->Render();
}


//-----------------------------------------------------------------------------
vtkSmartPointer<vtkVolumeProperty> mdlGUI::NewRGBVolumeProperty(const double range[])
{
	// Create transfer mapping scalar value to opacity.
	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacityTransferFunction->AddPoint( range[0], 0.0);
	opacityTransferFunction->AddPoint( range[1], 0.5);

	// Create transfer mapping scalar value to color.
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	colorTransferFunction->SetColorSpaceToRGB();
	colorTransferFunction->AddRGBPoint(range[0], 0.0, 0.0, 1.0);
	colorTransferFunction->AddRGBPoint(range[1], 1.0, 0.0, 0.0);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetColor(colorTransferFunction);
	volumeProperty->SetScalarOpacity(opacityTransferFunction);
	volumeProperty->ShadeOff();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetScalarOpacityUnitDistance(0.75);

	return volumeProperty;
}

//-----------------------------------------------------------------------------
void mdlGUI::RenderPoints(std::vector<mdl::fPoint3D> points, double r, double g, double b)
{
  //Remove all old actors
  //this->Renderer->RemoveAllViewProps();
  //Render the clean image so that the points have some context
  //this->RenderImage(this->Helper->GetImageData(), "anisotropic diffusion");
  //Change the window title so the user knows what file is being displayed
  //this->RenderWidget->setWindowTitle(windowName);

	vtkSmartPointer<vtkActor> actor = this->CreateActorFromPoints(points);
	actor->GetProperty()->SetColor(r, g, b);
	this->Renderer->AddActor(actor);

  //this->Renderer->ResetCamera();
  //this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkActor> mdlGUI::CreateActorFromPoints(std::vector<mdl::fPoint3D> points)
{
	//declare our return value
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();

	//initialize vtkPoints object for visualization
	vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();

	int numPoints = (int)points.size();
	vtk_points->Allocate(numPoints);
	for(int i = 0; i < numPoints; i++)
    {
		vtk_points->InsertPoint(i, points.at(i).x, points.at(i).y, points.at(i).z);
    }

	//finish setting up the actor
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(vtk_points);

	vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
	glyph->SetInput(polyData);
	vtkSmartPointer<vtkGlyphSource2D> glyphSource = vtkSmartPointer<vtkGlyphSource2D>::New();
	glyphSource->SetGlyphTypeToCross();
	glyphSource->SetScale(1.0);
	glyphSource->Update();
	glyph->SetInputConnection(1, glyphSource->GetOutputPort());

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput(glyph->GetOutput());

	actor->SetMapper(mapper);
	return actor;
}

//-----------------------------------------------------------------------------
void mdlGUI::RenderPolyData(std::string filename, double r, double g, double b)
{
  //remove all old actors
  //this->Renderer->RemoveAllViewProps();
  //change the window title so the user knows what file is being displayed
  //this->RenderWidget->setWindowTitle(windowName);

	vtkSmartPointer<vtkActor> actor = this->CreateActorFromPolyDataFile(filename.c_str());
	actor->GetProperty()->SetColor(r, g, b);
	this->Renderer->AddActor(actor);
  //this->Renderer->ResetCamera();
  //this->RenderWidget->GetRenderWindow()->Render();
}

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkActor> mdlGUI::CreateActorFromPolyDataFile(const char *filename)
{
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename);

	vtkPolyData *polyData = reader->GetOutput();

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput(polyData);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	return actor;
}


//-----------------------------------------------------------------------------
void mdlGUI::closeEvent(QCloseEvent *event)
{
  //this->OutputWindow->close();
  //this->HelpWindow->close();
  this->RenderWidget->close();
  event->accept();
}