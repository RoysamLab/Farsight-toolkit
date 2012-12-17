#include "ImageViewer.h"



ImageViewer::ImageViewer(QWidget *parent)
: QMainWindow(parent)
{
	
	std::cout<<"in ImageViewer constructor"<<endl;
	setWindowTitle(tr("Image Viewer"));
	
	this->CentralWidget = new QWidget(this);
	this->setCentralWidget(this->CentralWidget);

	
	this->accept = new QAction("Accept", this->CentralWidget);
	this->accept->setObjectName(tr("Accept"));
	connect(this->accept, SIGNAL(triggered()), this, SLOT(saveImage()));
	connect(this->accept, SIGNAL(triggered()), this, SLOT(closeWindow()));
	


	this->reject = new QAction("Reject",this->CentralWidget);
	this->reject->setObjectName(tr("Reject"));
	connect(this->reject, SIGNAL(triggered()), this, SLOT(closeWindow()));
	

	EditsToolBar = addToolBar(tr("Edit Toolbar"));
	this->EditsToolBar->addAction(this->accept);
	this->EditsToolBar->addAction(this->reject);

	this->EditsToolBar->addWidget(new QLabel("Cost Threshold: "));
	this->costThresholdLabel = new QLabel(this);	
	this->EditsToolBar->addWidget(this->costThresholdLabel);

	this->EditsToolBar->addWidget(new QLabel("Intensity Threshold: "));
	this->intensityThresholdLabel = new QLabel(this);	
	this->EditsToolBar->addWidget(this->intensityThresholdLabel);

	this->EditsToolBar->addWidget(new QLabel("Contrast Threshold: "));
	this->contrastThresholdLabel = new QLabel(this);	
	this->EditsToolBar->addWidget(this->contrastThresholdLabel);

	this->EditsToolBar->addWidget(new QLabel("Debris Threshold: "));
	this->debrisThresholdLabel = new QLabel(this);	
	this->EditsToolBar->addWidget(this->debrisThresholdLabel);


	this->EditsToolBar->addWidget(new QLabel("No. Of Trace Lines: "));
	this->noOfTraceLinesLabel= new QLabel(this);	
	this->EditsToolBar->addWidget(this->noOfTraceLinesLabel);

	RaycastBar = new QToolBar("RayCast Tools", this);
	this->RaycastBar->setObjectName("RaycastBar");
	this->RaycastBar->setAllowedAreas(Qt::BottomToolBarArea);
	//this->RaycastBar->setMovable(false);
	this->RaycastBar->setToolTip("Raycaster settings");
	this->addToolBarBreak(Qt::BottomToolBarArea);
	this->addToolBar(Qt::BottomToolBarArea,this->RaycastBar);
	
	
	QStringList ColorProfileList;
	ColorProfileList << "RGB" << "Red" << "Green" << "Blue" << "Gray";
	this->ColorProfileCombo = new QComboBox;
	this->ColorProfileCombo->setObjectName("ColorProfileCombo");
	this->ColorProfileCombo->addItems(ColorProfileList);
	connect(this->ColorProfileCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(RayCastColorValueChanged(int)));

	this->OpacityValueSpin = new QDoubleSpinBox(this);
	this->OpacityValueSpin->setObjectName("OpacityValueSpin");
	this->OpacityValueSpin->setRange(0,1);
	this->OpacityValueSpin->setSingleStep(.01);
	this->OpacityValueSpin->setValue(0.1);
	connect (this->OpacityValueSpin, SIGNAL(valueChanged(double)), this, SLOT(RayCastOpacityValueChanged(double)));

	this->BrightnessSpin = new QSpinBox(this);
	this->BrightnessSpin->setObjectName("BrightnessSpin");
	this->BrightnessSpin->setRange(0,250);
	this->BrightnessSpin->setValue(240);
	connect (this->BrightnessSpin, SIGNAL(valueChanged(int)), this , SLOT(RayCastBrightnessChanged(int)));

	this->OpacitySpin = new QSpinBox(this);
	this->OpacitySpin->setObjectName("OpacitySpin");
	this->OpacitySpin->setRange(0,250);
	

	this->OpacitySlider = new QSlider(Qt::Horizontal);
	this->OpacitySlider->setObjectName("OpacitySlider");
	this->OpacitySlider->setRange(0,250);
	this->OpacitySlider->setSingleStep(1);
	this->OpacitySlider->setTickInterval(5);
	this->OpacitySlider->setTickPosition(QSlider::TicksAbove);

	connect (this->OpacitySlider, SIGNAL(valueChanged(int)), this->OpacitySpin, SLOT(setValue(int)));
	this->OpacitySlider->setValue(240);
	connect (this->OpacitySpin, SIGNAL(valueChanged(int)), this->OpacitySlider, SLOT(setValue(int)));
	connect (this->OpacitySpin, SIGNAL(valueChanged(int)), this, SLOT(RayCastOpacityChanged(int)));


	this->RaycastBar->addWidget(new QLabel("Color Profile"));
	this->RaycastBar->addWidget(this->ColorProfileCombo);
	this->RaycastBar->addWidget(new QLabel("Opacity Threshold"));
	this->RaycastBar->addWidget(this->OpacitySpin);
	this->RaycastBar->addWidget(this->OpacitySlider);
	this->RaycastBar->addSeparator();
	this->RaycastBar->addWidget(new QLabel("Opacity Value"));
	this->RaycastBar->addWidget(this->OpacityValueSpin);
	this->RaycastBar->addSeparator();
	this->RaycastBar->addWidget(new QLabel("Brightness"));
	this->RaycastBar->addWidget(this->BrightnessSpin);
	
	/*accept = new QMenu(tr("Accept"));
	menuBar()->addMenu(accept);
	
	reject = new QMenu(tr("Reject"));
	menuBar()->addMenu(reject);*/

	QVTK = new QVTKWidget(this->CentralWidget);
	
	QGridLayout *viewerLayout = new QGridLayout(this->CentralWidget);
	viewerLayout->addWidget(this->QVTK, 0, 0);
	
	Renderer = vtkSmartPointer<vtkRenderer>::New();
	ImageActors = new ImageRenderActors();
	tobj = new TraceObject;
	LineActor = vtkSmartPointer<vtkActor>::New();
	this->LineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();


	Renderer->SetBackground(0.5,0.5,0.5);
	QVTK->GetRenderWindow()->AddRenderer(Renderer);

	//this->setCentralWidget(QVTK);
	
}

ImageViewer::~ImageViewer()
{
		if(this->QVTK)
	{
		delete this->QVTK;
	}

}






void ImageViewer::RayCastOpacityChanged(int value)
{
	ImageActors->setOpacity(value);
	QVTK->GetRenderWindow()->Render();
}


void ImageViewer::RayCastBrightnessChanged(int value)
{
	this->ImageActors->setBrightness(value);
	this->QVTK->GetRenderWindow()->Render();
}

void ImageViewer::RayCastOpacityValueChanged(double value)
{
	this->ImageActors->setOpacityValue(value);
	this->QVTK->GetRenderWindow()->Render();
}

void ImageViewer::RayCastColorValueChanged(int value)
{
	this->ImageActors->setColorValues(value);
	this->QVTK->GetRenderWindow()->Render();
}




void ImageViewer::showImage(){

	//std::string inputFilename = "C:/Data/mnt_data/del/before_robj_test.tif";
	//std::cout<<"before Create Sphere"<<endl;
	////picutreFlowTest();
	//std::cout<<"done with Create Sphere"<<endl;
	std::string inputFilename = _featureThreshold->getCropImageFilePath();
	ImageActors->loadImage(inputFilename,"Image");
	ImageActors->setOpacityValue(0.15);
	ImageActors->setOpacity(241);
	ImageActors->setRenderStatus(-1, true);
	Renderer->AddVolume(ImageActors->RayCastVolume(-1));
	QVTK->GetRenderWindow()->Render();
	QVTK->show();
	QVTK->resize(1600, 1000);
	this->resize(1600, 1000);
}

void ImageViewer::showTraces(){

	//std::string inputFilename = "C:\\Data\\mnt_data\\del\\MNT1_test_ANT.swc";
	std::string inputFilename = _featureThreshold->getSWCFilePath();
	//std::cout<<_featureThreshold->getContrastThreshold()<<endl;
	tobj->ReadFromSWCFile((char*)inputFilename.c_str());
	
	
	std::vector<TraceLine*> tracelines = tobj->GetTraceLines();
	std::cout<<"trace lines size in show traces----"<<tracelines.size()<<std::endl;
	if(tracelines.size() != 0){
		this->noOfTraceLinesLabel->setText(QString::number(tracelines.size()));
		
		UpdateLineActor();
		LineActor->SetPickable(1);
		Renderer->AddActor(LineActor);
		QVTK->GetRenderWindow()->Render();
		QVTK->show();
	}
}
//void ImageViewer::showTracesXML(){
//
//	this->tobj->ReadFromRPIXMLFile((char*)inputFilename.c_str());
//	std::vector<TraceLine*> tracelines = tobj->GetTraceLines();
//	this->noOfTraceLinesLabel->setText(QString::number(tracelines.size()));
//	UpdateLineActor();
//	LineActor->SetPickable(1);
//	Renderer->AddActor(LineActor);
//	QVTK->GetRenderWindow()->Render();
//	QVTK->show();
//}

void ImageViewer::UpdateLineActor()
{
	this->poly_line_data = this->tobj->GetVTKPolyData();
	this->poly_line_data->Modified();
	this->LineMapper->SetInput(this->poly_line_data);
	this->LineActor->SetMapper(this->LineMapper);
	this->LineActor->GetProperty()->SetColor(0,1,0);
	this->LineActor->GetProperty()->SetPointSize(2);
	this->LineActor->GetProperty()->SetLineWidth(2.5);
}

void ImageViewer::UpdateBranchActor()
{
	this->poly = tobj->generateBranchIllustrator();
	this->polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->polymap->SetInput(this->poly);
	this->BranchActor = vtkSmartPointer<vtkActor>::New();
	this->BranchActor->SetMapper(this->polymap);
	this->BranchActor->SetPickable(0);
}


void ImageViewer::saveImage(){
	emit acceptImage(_featureThreshold);
}
void ImageViewer::closeWindow(){
	delete this->QVTK;	
	close();
}


void ImageViewer::setImageFeatureThreshold(ImageFeatureThreshold* featureThreshold){
	_featureThreshold = featureThreshold;
	this->costThresholdLabel->setText(QString::number(_featureThreshold->getCostThreshold()));
	this->intensityThresholdLabel->setText(QString::number(_featureThreshold->getintensity_threshold()));
	this->contrastThresholdLabel->setText(QString::number(_featureThreshold->getContrastThreshold()));
	this->debrisThresholdLabel->setText(QString::number(_featureThreshold->getDebrisThreshold()));
}
void ImageViewer::saveScreenShot(std::string outputFileName){
	// Code to save the screenshot
	vtkSmartPointer<vtkWindowToImageFilter> ScreenShot = vtkSmartPointer<vtkWindowToImageFilter>::New();
	ScreenShot->SetInput(QVTK->GetRenderWindow());
	ScreenShot->SetMagnification(0);
	ScreenShot->ShouldRerenderOff();
	ScreenShot->Update();
	
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(outputFileName.c_str());
	writer->SetInput(ScreenShot->GetOutput());
	writer->Write();

}
void ImageViewer::saveTraceFeatures(std::string outputFileName){

	this->TreeModel = new TraceModel(this->tobj->GetTraceLines(), this->tobj->FeatureHeaders);
	TreeModel->scaleFactor = 1.0;

	CellTraceModel *CellModel;
	CellModel = new CellTraceModel();
	std::map< int ,CellTrace*> NewCells = this->tobj->CalculateCellFeatures();
	if (NewCells.size() > 0)
	{
		CellModel->setCells(NewCells);
	}
	vtkSmartPointer<vtkTable> table = CellModel->getDataTable();
	//	table->Dump(1);
	
	//Dump out headers
	ofstream myfile;
	myfile.open(outputFileName.c_str());
	for(vtkIdType columnIndex = 0; columnIndex < table->GetNumberOfColumns(); columnIndex++ )
	{	
		myfile << table->GetColumnName(columnIndex) << "\t";
	}
	myfile << std::endl;

	//Dump out data
	for(vtkIdType rowIndex = 0; rowIndex < table->GetNumberOfRows(); rowIndex++ )
	{
		for(vtkIdType columnIndex = 0; columnIndex < table->GetNumberOfColumns(); columnIndex++ )
		{
			vtkVariant value = table->GetValue(rowIndex, columnIndex);
			myfile << value << "\t";
		}
		myfile << endl;
	}
	myfile.close();

}



QStringList ImageViewer::findFiles(const QString& path = QString())
{
  QStringList files;

  QDir dir = QDir::current();
  if(!path.isEmpty())
    dir = QDir(path);

  dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
#if QT_VERSION >= 0x040000
  QFileInfoList list = dir.entryInfoList();
  for (int i = 0; i < list.size(); ++i) 
  {
    QFileInfo fileInfo = list.at(i);
    files.append(dir.absoluteFilePath(fileInfo.fileName()));
  }
#else
  const QFileInfoList* list = dir.entryInfoList();
  if(list) 
  {
    QFileInfoListIterator it( *list );
    QFileInfo * fi;
    while( (fi=it.current()) != 0 ) 
    {
      ++it;
      files.append(dir.absFilePath(fi->fileName()));
    }
  }
#endif

  return files;
}


