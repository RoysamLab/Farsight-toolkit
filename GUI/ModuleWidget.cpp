#include "ModuleWidget.h"

//************************************************************************************
// ModuleWidget
//
// This widget is used for interacting with a module (usually segmentation).  It allows
// for the data image to be selected, the parameters to be set, and each module to be
// executed independently and results shown.
//
// When this widget is closed all associated views are closed as well
//
// The buttons in the widget are populated dynamically based on querying the module
// for the list of modules available, as well as other display options.
//
//************************************************************************************

//Constructor
ModuleWidget::ModuleWidget(vector<string> filenames)
{
	//setAttribute ( Qt::WA_DeleteOnClose );
	segmentation = NULL;
	msgBox = NULL;
	segWin = NULL;

	imageFileCombo = NULL;
	paramFileCombo = NULL;
	possibleFilenames = filenames;
	lastPath = ".";
	
	segmentation = new ftk::NuclearSegmentation( "", "" );

	SetupUi();

	if (imageFileCombo->count() == 1)
	{
		ImageBrowse(tr("Browse..."));
	}
	if (paramFileCombo->count() == 1)
	{
		ParamBrowse(tr("Browse..."));
	}
}

void ModuleWidget::SetupUi(void)
{
	this->setWindowTitle(tr(segmentation->getPackageName().c_str()));

	QGridLayout *gLayout = new QGridLayout();
	int layoutRow = 0;	//Used to keep track of current row in the layout
	QButtonGroup *group = new QButtonGroup();
	gLayout->setSpacing(1);

	//Setup the image filenameCombo
	QLabel *dataLabel = new QLabel(tr("Image File: "));
	gLayout->addWidget(dataLabel,layoutRow,0,1,1);
	imageFileCombo = new QComboBox();
	for (int i=0; i<possibleFilenames.size(); ++i)
	{
		imageFileCombo->addItem(QString::fromStdString(possibleFilenames.at(i)));
	}
	imageFileCombo->addItem(tr("Browse..."));
	connect(imageFileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ImageBrowse(QString)));
	gLayout->addWidget(imageFileCombo,layoutRow++,1,1,3);

	//Setup the Parameter filenameCombo
	QLabel *paramLabel = new QLabel(tr("Parameter File: "));
	gLayout->addWidget(paramLabel,layoutRow,0,1,1);
	paramFileCombo = new QComboBox();
	paramFileCombo->addItem(tr("Browse..."));
	connect(paramFileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));
	gLayout->addWidget(paramFileCombo,layoutRow++,1,1,3);

	//Place Init Button
	QPushButton *initButton = new QPushButton(tr("Set Input and Parameters"));
	connect(initButton,SIGNAL(clicked()),this,SLOT(Initialize()));
	gLayout->addWidget(initButton,layoutRow,0,1,2);

	//Place the Save Button
	QPushButton *saveButton = new QPushButton(tr("Save Results"));
	connect(saveButton,SIGNAL(clicked()),this, SLOT(saveResults()));
	gLayout->addWidget(saveButton,layoutRow++,2,1,2);

	//Setup module buttons
	vector<string> moduleNames = segmentation->getModuleNames();
	int numModules = (int)moduleNames.size();
	for (int i=0; i<numModules; ++i)
	{
		QPushButton *button = new QPushButton;
		button->setText(tr(moduleNames[i].c_str()));
		group->addButton(button,i);
		gLayout->addWidget(button,layoutRow++,0,1,1);
	}

	//Place Features/XML button
	QPushButton *xmlButton = new QPushButton(tr("Generate Features and Export"));
	connect(xmlButton,SIGNAL(clicked()),this, SLOT(featToXML()));
	gLayout->addWidget(xmlButton,layoutRow,0,1,4);

	//Setup message Window
	msgBox = new QTextBrowser();
	gLayout->addWidget(msgBox,layoutRow-numModules,1,numModules,3);

	connect(group,SIGNAL(buttonClicked(int)),\
		    this,SLOT(runModule(int)));

	this->setLayout(gLayout);
}

//******************************************************************************
// Create a new segmentation window.  Shows Images and Segmentation results.
//******************************************************************************
void ModuleWidget::CreateNewSegmentationWindow(void)
{
	segWin = new SegmentationWindow();
	connect(segWin, SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	segWin->setWindowTitle(tr("Segmentation Results Viewer"));
	//segWin->show();
}

//*************************************************************************
// THIS SLOT IS CALLED WHEN A WINDOW IS CLOSED
//*************************************************************************
void ModuleWidget::closeWidget(QWidget *widget)
{
	if(widget == segWin)
	{
		segWin = NULL;
	}
}

//******************************************************************************
//Reimplement closeEvent to also close other module windows
//******************************************************************************
void ModuleWidget::closeEvent(QCloseEvent *event)
{
	//First Close other widgets
	closeChildren();
	//Then close myself
	emit closing(this);
	event->accept();
} 

void ModuleWidget::closeChildren(void)
{
	if(segWin)
	{
		segWin->close();
		segWin = NULL;
	}
}

void ModuleWidget::ImageBrowse(QString comboSelection)
{
	//First check to see if we have selected browse
	if (comboSelection != tr("Browse..."))
		return;

	QString newfilename  = QFileDialog::getOpenFileName(this,"Choose an Image File",lastPath, 
			tr("TIF Files (*.tif *.tiff)\n" 
			   "PIC Files (*.pic)\n" 
			   "LSM Files (*.lsm)\n"));

	if (newfilename == "")
	{
		imageFileCombo->setCurrentIndex(imageFileCombo->count()-2);
		return;
	}

	lastPath = QFileInfo(newfilename).absolutePath();
	int index = imageFileCombo->count()-1;
	imageFileCombo->setCurrentIndex(index-1);
	imageFileCombo->insertItem(index,newfilename);
	imageFileCombo->setCurrentIndex(index);
}

void ModuleWidget::ParamBrowse(QString comboSelection)
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
		paramFileCombo->setCurrentIndex(paramFileCombo->count()-2);
		return;
	}

	lastPath = QFileInfo(newfilename).absolutePath();
	int index = paramFileCombo->count()-1;
	paramFileCombo->setCurrentIndex(index-1);
	paramFileCombo->insertItem(index,newfilename);
	paramFileCombo->setCurrentIndex(index);
}

void ModuleWidget::Initialize(void)
{
	QString imfile = imageFileCombo->currentText();
	QString pafile = paramFileCombo->currentText();

	if ( (imfile == "") || (pafile == "") ) 
	{
		msgBox->insertPlainText(tr("Missing Input\r\n"));
	}
	else
	{
		
		QString path = QFileInfo(imfile).absolutePath();
		QString name = QFileInfo(imfile).baseName();
		QString img = QFileInfo(imfile).fileName();
		QString param = QFileInfo(pafile).fileName();

		msgBox->insertPlainText(tr("Set Image to..."));
		msgBox->insertPlainText(img);
		msgBox->insertPlainText(tr("\r\n"));
		msgBox->insertPlainText(tr("Set Params to..."));
		msgBox->insertPlainText(param);
		msgBox->insertPlainText(tr("\r\n"));
		msgBox->ensureCursorVisible();

		segmentation->SetProjectPath(path.toStdString());
		segmentation->SetProjectName(name.toStdString());
		segmentation->setup(img.toStdString(), param.toStdString());

		closeChildren();
	}
}

void ModuleWidget::saveResults(void)
{
	msgBox->insertPlainText(tr("Saving Results to disk\r\n"));
	msgBox->ensureCursorVisible();
	segmentation->SaveLabel();
}

void ModuleWidget::featToXML(void)
{
	msgBox->insertPlainText(tr("Generating features and writing to XML\r\n"));
	msgBox->ensureCursorVisible();
	if ( !segmentation->LabelsToObjects() )
		cerr << segmentation->GetErrorMessage() << endl;
	if ( !segmentation->WriteToXML() )
		cerr << segmentation->GetErrorMessage() << endl;
	//if ( !segmentation->WriteToMETA() )
	//	cerr << segmentation->GetErrorMessage() << endl;
	if ( !segmentation->WriteToLibSVM() )
		cerr << segmentation->GetErrorMessage() << endl;
}

void ModuleWidget::runModule(int moduleNum)
{
	msgBox->insertPlainText(tr("Running "));
	msgBox->insertPlainText(QString::fromStdString(segmentation->getModuleNames().at(moduleNum)));
	msgBox->insertPlainText(tr("\r\n"));
	msgBox->ensureCursorVisible();

	segmentation->executeModule(moduleNum);

	showResults();
}

void ModuleWidget::showResults(void)
{
	if(!segmentation)
		return;

	//Check for something to show in segmentation window, if so make sure it gets shown
	if ( segmentation->getDataImage() || segmentation->getLabelImage() )
	{
		if(!segWin)
			CreateNewSegmentationWindow();
		segWin->SetChannelImage(segmentation->getDataImage());
		segWin->SetLabelImage(segmentation->getLabelImage());
		segWin->show();
	}
	else
	{
		if(segWin)
		{
			segWin->close();
			segWin = NULL;
		}
	}
}
