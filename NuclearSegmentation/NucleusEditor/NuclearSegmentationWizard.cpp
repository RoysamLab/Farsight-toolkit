#include "NuclearSegmentationWizard.h"

NuclearSegmentationWizard::NuclearSegmentationWizard(QWidget *parent)
	: QWizard(parent)
{
	this->setPage(Page_Input, new InputPage);
	this->setPage(Page_Parameters, new ParametersPage);
	this->setPage(Page_Binarize, new BinarizePage);
	this->setPage(Page_Seeds, new SeedsPage);
	this->setPage(Page_Cluster, new ClusterPage);
	this->setPage(Page_Finalize, new FinalizePage);
	this->setPage(Page_Exit, new ExitPage);

	this->setStartId(Page_Input);
	//this->setModal(true);

	this->setOption(QWizard::HaveCustomButton1);
	this->setButtonText(QWizard::CustomButton1,"Execute");
	connect(this, SIGNAL(customButtonClicked(int)), this, SLOT(executeNextStep(int)));

	//setOption(HaveHelpButton, true);
	//setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	//connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));

	this->setWindowTitle(tr("Nuclear Segmentation Wizard"));

	seg = new ftk::NuclearSegmentation();
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
	case Page_Exit:
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
	case Page_Exit:
		button(QWizard::CustomButton1)->setVisible(false);
		break;
	case Page_Seeds:
	case Page_Cluster:
	case Page_Finalize:
		if( seg )
			button(QWizard::CustomButton1)->setVisible(true);
		else
			button(QWizard::CustomButton1)->setVisible(false);
		break;
	}
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
		((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL );
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
			((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			break;
		case Page_Seeds:
			seg->DetectSeeds();
			((SeedsPage*)page(Page_Seeds))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			break;
		case Page_Cluster:
			seg->RunClustering();
			((ClusterPage*)page(Page_Cluster))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			break;
		case Page_Finalize:
			seg->Finalize();
			((FinalizePage*)page(Page_Finalize))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			break;
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
		return Page_Seeds;
		break;
	case Page_Seeds:
		return Page_Cluster;
		break;
	case Page_Cluster:
		return Page_Finalize;
		break;
	case Page_Finalize:
		return Page_Exit;
		break;
	case Page_Exit:
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
	setLayout(layout);
	hasImage = false;
}

void ClusterPage::ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label)
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
