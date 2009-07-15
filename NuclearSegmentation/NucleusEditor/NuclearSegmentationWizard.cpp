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
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Seeds:
			seg->DetectSeeds();
			((SeedsPage*)page(Page_Seeds))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
			((ClusterPage*)page(Page_Cluster))->ShowImages( NULL, NULL );
			((FinalizePage*)page(Page_Finalize))->ShowImages( NULL, NULL );
			((SavePage*)page(Page_Save))->ImagesSaved(false);
			break;
		case Page_Cluster:
			seg->RunClustering();
			((ClusterPage*)page(Page_Cluster))->ShowImages( seg->getDataImage(), seg->getLabelImage() );
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
	registerField("cluster.jump",jumpBox);
	layout->addWidget(jumpBox);
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
