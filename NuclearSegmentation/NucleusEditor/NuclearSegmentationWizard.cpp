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

#include "NuclearSegmentationWizard.h"
#include "Seed3DHelperClasses.h"

NuclearSegmentationWizard::NuclearSegmentationWizard(QWidget *parent)
	: QWizard(parent)
{
	this->setPage(Page_Input, new InputPage);
	this->setPage(Page_Parameters, new ParametersPage);
	this->setPage(Page_Binarize, new BinarizePage);
	this->setPage(Page_Seeds, new SeedsPage);
	this->setPage(Page_Cluster, new ClusterPage);
	this->setPage(Page_Finalize, new FinalizePage);

	SavePage *sp = new SavePage;
	connect(sp, SIGNAL(readyToSave(bool)), this, SLOT(updateExeButton(bool)));
	this->setPage(Page_Save, sp);

	this->setStartId(Page_Input);
	//this->setModal(true);

	this->setOption(QWizard::HaveCustomButton1);
	this->setButtonText(QWizard::CustomButton1,"Execute");
	connect(this, SIGNAL(customButtonClicked(int)), this, SLOT(executeNextStep(int)));

	this->setOption(QWizard::NoBackButtonOnStartPage,true);
	//setOption(HaveHelpButton, true);
	//setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	//connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));

	this->setWindowTitle(tr("Nuclear Segmentation Wizard"));

	seg = new ftk::NuclearSegmentation();
	Seeds = new Seed3D(0,0);	
 }

//is called by QWizard to prepare page id just before it is shown as a result of the user clicking Next
void NuclearSegmentationWizard::initializePage(int id)
{
	switch(id)
	{
	case Page_Input:
		button(QWizard::CustomButton1)->setVisible(false);
		break;
	case Page_Parameters:
		button(QWizard::CustomButton1)->setVisible(false);
		break;
	case Page_Binarize:
		if( this->initSegmentation() )
			button(QWizard::CustomButton1)->setVisible(true);
		else
			button(QWizard::CustomButton1)->setVisible(false);
		break;
	case Page_Seeds:
		button(QWizard::CustomButton1)->setVisible(true);
		break;
	case Page_Cluster:
		button(QWizard::CustomButton1)->setVisible(true);
		break;
	case Page_Finalize:
		button(QWizard::CustomButton1)->setVisible(true);
		break;
	case Page_Save:
		button(QWizard::CustomButton1)->setVisible(false);
		break;
	}
}

//is called by QWizard to clean up page id just before the user leaves it by clicking Back
void NuclearSegmentationWizard::cleanupPage(int id)
{
	switch(id)
	{
	case Page_Input:
	case Page_Parameters:
	case Page_Binarize:
		button(QWizard::CustomButton1)->setVisible(false);
		break;
	case Page_Seeds:
	case Page_Cluster:
	case Page_Finalize:
	case Page_Save:
		if( seg )
			button(QWizard::CustomButton1)->setVisible(true);
		else
			button(QWizard::CustomButton1)->setVisible(false);
		break;
	}
}

void NuclearSegmentationWizard::updateExeButton(bool val)
{
	button(QWizard::CustomButton1)->setVisible(val);
}

bool NuclearSegmentationWizard::initSegmentation(void)
{
	//Also can get the segmentation ready
	QString dataFile = field("InputFile").toString();
	QString paramFile = field("ParamFile").toString();

	if(dataFile == "" || paramFile == "")
		return false;

	if(dataFile.toStdString() == seg->GetDataFilename() && paramFile.toStdString() == seg->GetParamFilename())
	{
		return true;
	}
	else
	{
		((BinarizePage*)page(Page_Binarize))->ShowImages( NULL, NULL );
		((SeedsPage*)page(Page_Seeds))->ShowImages( NULL, NULL );
		((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL,NULL );
		((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
		seg->SetInputs( dataFile.toStdString(), paramFile.toStdString() );
		if(seg->LoadData())
			return true;
		else
			return false;
	}

	return true;
}

void NuclearSegmentationWizard::executeNextStep(int whichButton)
{
	if(whichButton == QWizard::CustomButton1)
	{
		if(!seg)
			return;

		//Find out where in the process I am and execute the appropriate step.
		int page_id = this->currentId();

		switch(page_id)
		{
		case Page_Binarize:
			seg->Binarize();
			((BinarizePage*)page(Page_Binarize))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((SeedsPage*)page(Page_Seeds))->ShowImages( NULL, NULL );
			((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL,NULL );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Seeds:
			seg->DetectSeeds();
			((SeedsPage*)page(Page_Seeds))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL,NULL );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Cluster:
			seg->RunClustering();
			((ClusterPage*)page(Page_Cluster))->ShowImages( seg->getDataImage(), seg->getLabelImage(), Seeds );
			//The next two lines are for the Seed Editor			
			//Seeds->GetImage(seg->getDataImage(),seg->getSeeds());
			Seeds->GetImage(seg,seg->getSeeds());
			((ClusterPage*)page(Page_Cluster))->SeedEdit->setEnabled(true);
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Finalize:
			seg->Finalize();
			((FinalizePage*)page(Page_Finalize))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Save:
			seg->SaveLabel();
			if(field("save.xmlRadio").toBool())
			{
				seg->LabelsToObjects();
				seg->WriteToXML( field("save.xmlFile").toString().toStdString() );
			}
			((SavePage*)page(Page_Save))->ImagesSaved(true);
		}
	}
}

int NuclearSegmentationWizard::nextId() const
{
	switch( this->currentId() )
	{
	case Page_Input:
		return Page_Parameters;
		break;
	case Page_Parameters:
		return Page_Binarize;
		break;
	case Page_Binarize:
		if(field("binarize.jump").toBool())
			return Page_Save;
		else
			return Page_Seeds;
		break;
	case Page_Seeds:
		if(field("seeds.jump").toBool())
			return Page_Save;
		else
			return Page_Cluster;
		break;
	case Page_Cluster:
		if(field("cluster.jump").toBool())
			return Page_Save;
		else
			return Page_Finalize;
		break;
	case Page_Finalize:
		return Page_Save;
		break;
	case Page_Save:
		return -1;
		break;
	}
	return -1;
}


//****************************************************************************
// PAGES:
//****************************************************************************

InputPage::InputPage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Input"));
	//setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.png"));
	
	QVBoxLayout *layout = new QVBoxLayout;

	topLabel = new QLabel(tr("Please choose an input image that you would like to segment:"));
	topLabel->setWordWrap(true);
	layout->addWidget(topLabel);

	imageFileCombo = new QComboBox();
	imageFileCombo->addItem(tr(""));
	imageFileCombo->addItem(tr("Browse..."));
	registerField("InputFile",imageFileCombo,"currentText", SIGNAL(currentIndexChanged(QString)));
	connect(imageFileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ImageBrowse(QString)));
	layout->addWidget(imageFileCombo);

	sWin = new SegmentationWindow();
	layout->addWidget(sWin);

	setLayout(layout);
}

bool InputPage::isComplete() const
{
	if(imageFileCombo->currentText() != "")
		return true;
	else
		return false;
}

void InputPage::ImageBrowse(QString comboSelection)
{
	//First check to see if we have selected browse
	if (comboSelection != tr("Browse..."))
		return;

	QString newfilename  = QFileDialog::getOpenFileName(this,"Choose an Image File",lastPath, 
			tr("PIC Files (*.pic)\n"
			   "TIF Files (*.tif *.tiff)\n"  
			   "LSM Files (*.lsm)\n"
			   "All Files (*.*)\n"));

	if (newfilename == "")
	{
		imageFileCombo->setCurrentIndex(0);
		return;
	}

	if( newfilename == imageFileCombo->currentText() )
		return;

	lastPath = QFileInfo(newfilename).absolutePath();
	imageFileCombo->setItemText(0,newfilename);
	imageFileCombo->setCurrentIndex(0);

	ftk::Image::Pointer img = ftk::Image::New();
	img->LoadFile( imageFileCombo->currentText().toStdString() );
	sWin->SetChannelImage( img );

	emit completeChanged();
}


ParametersPage::ParametersPage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Parameters"));

	QVBoxLayout *layout = new QVBoxLayout;

	paramLabel = new QLabel(tr("Please choose a parameters file:"));
	paramLabel->setWordWrap(true);
	layout->addWidget(paramLabel);

	paramFileCombo = new QComboBox();
	paramFileCombo->addItem(tr(""));
	paramFileCombo->addItem(tr("Browse..."));
	registerField("ParamFile",paramFileCombo,"currentText", SIGNAL(currentIndexChanged(QString)));
	connect(paramFileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));
	layout->addWidget(paramFileCombo);

	fileText = new QTextBrowser();
	fileText->setWordWrapMode(QTextOption::NoWrap);
	layout->addWidget(fileText);

	setLayout(layout);
}

bool ParametersPage::isComplete() const
{
	if(paramFileCombo->currentText() != "")
		return true;
	else
		return false;
}

void ParametersPage::ParamBrowse(QString comboSelection)
{
	//First check to see if we have selected browse
	if (comboSelection != tr("Browse..."))
		return;

	QString newfilename  = QFileDialog::getOpenFileName(this,"Choose a Parameters File",lastPath, 
			tr("INI Files (*.ini)\n" 
			   "TXT Files (*.txt)\n" 
			   "All Files (*.*)\n"));

	if (newfilename == "")
	{
		paramFileCombo->setCurrentIndex(0);
		return;
	}

	if( newfilename == paramFileCombo->currentText() )
		return;

	lastPath = QFileInfo(newfilename).absolutePath();
	paramFileCombo->setCurrentIndex(0);
	paramFileCombo->setItemText(0,newfilename);

	//Now read the text from the file:
	QFile file(newfilename);
	if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		return;
	}
	QTextStream in(&file);
	QString all = in.readAll();
	fileText->setText(all);

	emit completeChanged();
}

BinarizePage::BinarizePage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Binarization"));
	QVBoxLayout *layout = new QVBoxLayout;
	sWin = new SegmentationWindow();
	layout->addWidget(sWin);
	jumpBox = new QCheckBox("Skip to Calculating Features and Saving Result");
	registerField("binarize.jump",jumpBox);
	layout->addWidget(jumpBox);
	setLayout(layout);
	hasImage = false;
}

bool BinarizePage::isComplete() const
{
	if(hasImage)
		return true;
	else
		return false;
}

void BinarizePage::ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label)
{
	sWin->SetChannelImage(data);
	sWin->SetLabelImage(label);
	ftk::Image::Pointer null = NULL;
	if(data == null && label == null)
		hasImage = false;
	else
		hasImage = true;
	emit completeChanged();
}

SeedsPage::SeedsPage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Seed Detection"));
	QVBoxLayout *layout = new QVBoxLayout;
	sWin = new SegmentationWindow();
	layout->addWidget(sWin);
	jumpBox = new QCheckBox("Skip to Calculating Features and Saving Result");
	registerField("seeds.jump",jumpBox);
	layout->addWidget(jumpBox);
	setLayout(layout);
	hasImage = false;
}

bool SeedsPage::isComplete() const
{
	if(hasImage)
		return true;
	else
		return false;
}

void SeedsPage::ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label)
{
	sWin->SetChannelImage(data);
	sWin->SetLabelImage(label);
	ftk::Image::Pointer null = NULL;
	if(data == null && label == null)
		hasImage = false;
	else
		hasImage = true;
	emit completeChanged();
}

bool ClusterPage::isComplete() const
{
	if(hasImage)
		return true;
	else
		return false;
}

ClusterPage::ClusterPage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Clustering"));
	QVBoxLayout *layout = new QVBoxLayout;
	sWin = new SegmentationWindow();
	layout->addWidget(sWin);
	jumpBox = new QCheckBox("Skip to Calculating Features and Saving Result");
	SeedEdit = new QPushButton("&View Seeds in 3D and Edit", this);
	SeedEdit->setDisabled(true);
	connect(SeedEdit, SIGNAL(clicked()), this, SLOT(ViewEditor()));
	registerField("cluster.jump",jumpBox);
	layout->addWidget(jumpBox);
	layout->addWidget(SeedEdit);
	setLayout(layout);
	hasImage = false;
}



 void ClusterPage::ViewEditor()
{
 	
	seeds2->show();
 
}



 void ClusterPage::ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label,Seed3D* seeds)
{
	sWin->SetChannelImage(data);
	sWin->SetLabelImage(label);
	
	seeds2 =  seeds;
	//Raghav for Seed Editor
	//Seeds->LoadSeedEditor(data);

	ftk::Image::Pointer null = NULL;
	if(data == null && label == null)
		hasImage = false;
	else
		hasImage = true;
	emit completeChanged();

	
}


//void Seed3D::GetImage(ftk::Image::Pointer data, vector<Seed> seeds)
void Seed3D::GetImage(ftk::NuclearSegmentation *seg,vector<Seed> seeds) //modified by Yousef
{
	if(this->iRender==1)
		DeleteObjects();
   segPtr = seg; 
   
   vtkSmartPointer<vtkImageData> vtkim = vtkSmartPointer<vtkImageData>::New();   
   //vtkim = data->GetVtkPtr(0,0);
   vtkim = segPtr->getDataImage()->GetVtkPtr(0,0);
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

	
// // Create seeds (spheres) only at coordinates specified by the text file.
// // I also create Red (Markedpoints) and Green Spheres (Markedpoints2add) initially and delete them.
////  So that the first delete and addition by the user would be 
//// "inserts" into the vector !
//// Note: 1) All the seeds form a single glyph.
////       2) The green spheres form a single glyph. These are seeds added but not applied permanently.
////	   3) The red spheres form a single glyph.These seeds are marked for Deletion but can be unmarked

//	//spPoint is the vector of seeds 
	spPoint = GetSeedpts(seeds);
    this->dup_points = spPoint; //dup_points contains the upto date location of seeds 
		
    // Create a float array which represents the points.
     this->delpcoords = vtkFloatArray::New();	
    
	

 	this->pcoords = vtkFloatArray::New();
 /*  Note that by default, an array has 1 component.
   We have to change it to 3 for points
    */
    this->pcoords->SetNumberOfComponents(3);
    this->pcoords->SetNumberOfTuples(spPoint.size());
    this->delpcoords->SetNumberOfComponents(3);
    this->delpcoords->SetNumberOfTuples(1);		
     for (unsigned int j=0; j<spPoint.size(); j++)
    {
    float pts[3] = {spPoint[j].x,spPoint[j].y , spPoint[j].z };	
    this->pcoords->SetTuple(j, pts);
    }
    float pts[3] = {0.0,0.0,0.0};
    
	this->delpcoords->SetTuple(0,pts);	
     
// Create vtkPoints and assign pcoords as the internal data array.
  this->point1 = vtkPoints::New(); //Original Points
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

   
// create the glyphing the sphere with a cone.  Create the mapper
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
    

	// The Pipeline
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
  this->sliderWidget3->AddObserver(vtkCommand::InteractionEvent,callback_seedsize);
  this->sliderWidget3->EnabledOn();
  
  this->iRender = 1;
  this->QVTK->GetRenderWindow()->Render();   
  this->QVTK1->GetRenderWindow()->Render();   
  this->QVTK2->GetRenderWindow()->Render(); 

}

FinalizePage::FinalizePage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Finalization"));
	QVBoxLayout *layout = new QVBoxLayout;
	sWin = new SegmentationWindow();
	layout->addWidget(sWin);
	setLayout(layout);
	hasImage = false;
}

bool FinalizePage::isComplete() const
{
	if(hasImage)
		return true;
	else
		return false;
}

void FinalizePage::ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label)
{
	sWin->SetChannelImage(data);
	sWin->SetLabelImage(label);
	ftk::Image::Pointer null = NULL;
	if(data == null && label == null)
		hasImage = false;
	else
		hasImage = true;
	emit completeChanged();
}

SavePage::SavePage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Save Results"));
	QVBoxLayout *layout = new QVBoxLayout;

	topLabel = new QLabel(tr("The Result will now be saved. Please choose desired option:"));
	layout->addWidget(topLabel);

	imageOnlyRadio = new QRadioButton(tr("Save Label Image Only"));
	layout->addWidget(imageOnlyRadio);

	xmlRadio = new QRadioButton(tr("Compute Features and save XML File and Label Image"));
	registerField("save.xmlRadio",xmlRadio);
	connect(xmlRadio,SIGNAL(toggled(bool)),this, SLOT(radioChanged(bool)));
	layout->addWidget(xmlRadio);

	saveLabel = new QLabel(tr("Please choose XML Filename to save results as:"));
	layout->addWidget(saveLabel);

	xmlFileCombo = new QComboBox();
	xmlFileCombo->addItem(tr(""));
	xmlFileCombo->addItem(tr("Browse..."));
	registerField("save.xmlFile",xmlFileCombo,"currentText", SIGNAL(currentIndexChanged(QString)));
	connect(xmlFileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(SaveAsBrowse(QString)));
	layout->addWidget(xmlFileCombo);

	setLayout(layout);

	xmlRadio->setChecked(true);
	saved = false;
}

void SavePage::radioChanged(bool val)
{
	if(val == true)
	{
		xmlFileCombo->setEnabled(true);
		if( xmlFileCombo->currentText() != "" )
			//wizard()->button(QWizard::CustomButton1)->setVisible(true);
			emit readyToSave(true);
		else
			emit readyToSave(false);
			//wizard()->button(QWizard::CustomButton1)->setVisible(false);
	}
	else
	{
		xmlFileCombo->setEnabled(false);
		emit readyToSave(true);
		//wizard()->button(QWizard::CustomButton1)->setVisible(true);
	}
}

bool SavePage::isComplete() const
{
	if(saved)
		return true;
	else
		return false;
}

void SavePage::ImagesSaved(bool s)
{
	saved = s;
	emit completeChanged();
}

void SavePage::SaveAsBrowse(QString comboSelection)
{
	//First check to see if we have selected browse
	if (comboSelection != tr("Browse..."))
		return;

	QString newfilename  = QFileDialog::getSaveFileName(this,"Choose an XML File",lastPath, 
			tr("XML Files (*.xml)\n"));

	if (newfilename == "")
	{
		xmlFileCombo->setCurrentIndex(0);
		return;
	}

	if( newfilename == xmlFileCombo->currentText() )
		return;

	lastPath = QFileInfo(newfilename).absolutePath();
	xmlFileCombo->setItemText(0,newfilename);
	xmlFileCombo->setCurrentIndex(0);
	emit readyToSave(true);
}


//***********************************************************
//Seed Editor classes ( Seed3D) start here
//***********************************************************

Seed3D::Seed3D(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{

	////Create the GUI
	/*createMenus();
	createStatusBar();*/
	setWindowTitle(tr("Seed Editor"));
	lastPath = ".";

//Variable instantiation
	this->mode = 0;
	this->counter =0;
	this->stateAdd = 0;
	this->stateDelete = 0;
	this->stateSplit = 0;
	this->stateMerge = 0;
	this->flag =1;
	segPtr = NULL;
	
	
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
  this->AddBox = new QCheckBox("&Add Seeds", this->browse);
  this->DeleteBox = new QCheckBox("&Delete Seeds", this->browse);
  this->MergeBox = new QCheckBox("&Merge Nuclei", this->browse);
  this->SplitBox = new QCheckBox("&Split Nuclei", this->browse);
  this->UndoDelBox = new QCheckBox("&Undo Selections for Delete", this->browse);	
  this->PlaceButton = new QPushButton("&Place Seed", this->browse);
  this->ApplyButton = new QPushButton("&Apply", this->browse);  
  
//********************************************************************************
// Layouts 

  QGridLayout *buttonLayout = new QGridLayout(this);
  buttonLayout->addWidget(this->AddBox,10, 2);
  buttonLayout->addWidget(this->DeleteBox,10, 3);
  buttonLayout->addWidget(this->SplitBox,10, 4);
  buttonLayout->addWidget(this->MergeBox,10, 5);
  buttonLayout->addWidget(this->UndoDelBox,10, 6);

  buttonLayout->addWidget(this->PlaceButton,11, 2);
  buttonLayout->addWidget(this->ApplyButton,11, 3);

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
  connect(this->AddBox, SIGNAL(stateChanged(int)), this, SLOT(AddSeed()));
  connect(this->DeleteBox, SIGNAL(stateChanged(int)), this, SLOT(DeleteSeed()));
  connect(this->ApplyButton, SIGNAL(clicked()), this, SLOT(Apply()));		
  connect(this->UndoDelBox, SIGNAL(stateChanged(int)), this, SLOT(UndoDeleteSeeds()));
  connect(this->SplitBox, SIGNAL(stateChanged(int)), this, SLOT(SplitSeeds()));
  connect(this->MergeBox, SIGNAL(stateChanged(int)), this, SLOT(MergeSeeds()));

//Resize the Window 
  this->resize(800,400);
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
    cout<<"Point Selected " << pickPos[0]<<"-"<<pickPos[1]<<"-"<<pickPos[2]<<endl;

    //Useful to check if clicked on a seed !   
//If I am in Delete mode 

if((seed->mode == 2||seed->mode == 4) && seed->flag == 1){

    vtkDataArray* pointIds = seed->Glyph->GetOutput()->GetPointData()->GetArray("InputPointIds"); 
    int pID = (int)pointIds->GetTuple1(seed->PointPicker->GetPointId()); 
    cout<<pID<<" is the iD value"<<endl;		

if((unsigned int)pID<=seed->dup_points.size())    //The ids of non-seed points is much greater than the ids of the seed points 
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
	cout<<index<<" is the inD value"<<endl;
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

if(seed->MarkedPoints.size()==seed->spPoint.size() || seed->dup_points.size()==0) 
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

if((unsigned int)pID<=seed->MarkedPoints.size())    //The ids of non-seed points is much greater than the ids of the seed points 

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
seed->Check();	
}





void Seed3D::PlaceSeed()
{

	if(this->mode==1 || this->mode ==3)
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
if(this->stateDelete)
this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes before performing this action?");
if(this->stateMerge)
this->msgBox->setText("You have marked seeds for merging. Do you want to apply the changes?"); 
this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
this->msgBox->setDefaultButton(QMessageBox::Yes);
int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
    if(this->stateDelete)
    this->mode =2;
    if(this->stateMerge)
    this->mode =4;
	this->Apply();
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
	if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
        this->mode = 1;
				 }
else{
this->mode=0;
}
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
		//Yousef_Code_Add_Seed
		//I will work on it later... (Yousef)
				
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
		//Yousef_Code_Delete_Seed
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
	//Yousef_Code_Split_Seeds
	//Assume that each pair of points (in order) represent two new seeds in a cell that we need to split
	for(int i=0; i<Id; i+=2)
	{
		std::vector< int > newIDs = segPtr->SplitInit(tobeSplit.at(i), tobeSplit.at(i+1));
		//newIDs will be used later for recording the edits
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
		//Yousef_Code_Merge_Seeds
		//I will work on it later... (Yousef)
	}

	vtkDataArray* points2merge = this->point2->GetData();    
	//Remove the glyph		
	for(int i =0 ;i<Id;i++)
	{		
		points2merge->RemoveTuple((vtkIdType)0);
		////cout<<"merge red"<<endl;		
	}
    this->point2->SetData(points2merge);
    this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001);    
    this->QVTK->GetRenderWindow()->Render();
    this->MarkedPoints.erase(this->MarkedPoints.begin(),this->MarkedPoints.end());
    //this->delglyph->SetScaleFactor(this->Glyph->GetScaleFactor()+0.0001); 	
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
if(this->stateAdd) { this->AddBox->setCheckState((Qt::CheckState)0); }////cout<<"hi"<<endl; }
if(this->stateUndoDel){this->UndoDelBox->setCheckState((Qt::CheckState)0); }
this->mode = 2;
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

void Seed3D::DeleteObjects(){
  
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->Volume);
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->AddSphereActor);	
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->DelSphereActor);
  this->QVTK->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this->SphereActor);
  
  this->Volume->Delete();
  this->AddSphereActor->Delete();
  this->DelSphereActor->Delete();	
  this->SphereActor->Delete();

  this->AddSphereMapper->Delete();
  this->DelSphereMapper->Delete();	
  this->SphereMapper->Delete();



  this->polydata1->Delete();
  this->polydata2->Delete();	
  this->polydata3->Delete();


  this->handle->Delete();
  this->handle1->Delete();
  this->handle2->Delete();
  this->widget->Delete();
  this->widget1->Delete();	
  this->widget2->Delete();

  this->sliderRep->Delete();
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
  this->QVTK2->GetRenderWindow()->Render();    
}

vector <ftk::Object::Point> Seed3D::GetSeedpts(vector<Seed> seeds)
{

	point p;
	Seed q;
	std::vector<point> spPoint;
	spPoint.clear();
	for (std::vector<Seed>::iterator i = seeds.begin (); i != seeds.end (); i++)
	{
		q = *i;
		p.x = q.x();
		
		p.y = q.y();
		
		p.z = q.z();
		
		p.t = 0;
		spPoint.push_back(p);
		cout<<spPoint.size()<<"--size"<<endl;
	}
	cout<<spPoint.size()<<"--Total size"<<endl;
	return (spPoint);
}

void Seed3D::SplitSeeds()
{

this->mode =3;
this->stateSplit = this->SplitBox->checkState();
if(MarkedPoints.size()!=0 && this->stateSplit){

	this->msgBox = new QMessageBox(this);
	this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes?");
	if(this->stateDelete)	
	this->msgBox->setText("You have marked seeds for deletion. Do you want to apply the changes?");
	if(this->stateMerge)	
	this->msgBox->setText("You have seeds marked for merging. Do you want to apply the changes?");
	this->msgBox->setStandardButtons(QMessageBox::Yes | QMessageBox::No);
	this->msgBox->setDefaultButton(QMessageBox::Yes);
	int ret = this->msgBox->exec();

switch (ret) {
   case QMessageBox::Yes:{
       // Save was clicked
   
   if(this->stateDelete)
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




