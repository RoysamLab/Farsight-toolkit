/*********************************************************************************************
// Screenshot dialog to browse for a folder, specify a filename, and save the screenshot.	//
// Additional functions: magnification factor.												//
*********************************************************************************************/

#include "screenshot.h"

ScreenShotDialog::ScreenShotDialog(QWidget* parent, QString fileName, QString imageDir)
: QDialog(parent)
{
	this->magnifynum = 1;

	QLabel *directoryLabel = new QLabel(tr("Directory:"));
	directoryComboBox = createComboBox();	
	directoryComboBox->setEditText(imageDir);
	directoryLabel->setBuddy(directoryComboBox);
	browseButton = createButton(tr("&Browse..."), SLOT(Browse()));

	QLabel *fileNameLabel = new QLabel(tr("Filename:"));
	fileNameLine = new QLineEdit();
	fileNameLabel->setBuddy(fileNameLine);

	QLabel * ZoomLabel = new QLabel(tr("Magnify factor: "));
	ZoomSpinBox = new QSpinBox();
	ZoomSpinBox->setValue(1);
	ZoomLabel->setBuddy(ZoomSpinBox);

	OkButton = new QPushButton(tr("Ok"));
	OkButton->setDefault(true);
	connect(OkButton, SIGNAL(clicked()), this, SLOT(save()));
	
	CancelButton = new QPushButton(tr("Cancel"));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(close()));

	// QGridLayout - the 1st number: initial row index, the 2nd number: initial column index, 
	//					the 3rd number: rowSpan, and the 4th number: columnSpan.
	QGridLayout *mainLayout = new QGridLayout;
	mainLayout->addWidget(directoryLabel,0,0);
	mainLayout->addWidget(directoryComboBox,0,1,1,2);
	mainLayout->addWidget(browseButton,0,3);
	mainLayout->addWidget(fileNameLabel,1,0);
	mainLayout->addWidget(fileNameLine,1,1,1,3);
	mainLayout->addWidget(ZoomLabel,2,0);
	mainLayout->addWidget(ZoomSpinBox,2,1);
	#ifdef USE_QT_TESTING
	this->BaselineBox = new QCheckBox("New testing baseline?", this);
	this->BaselineBox->setChecked(true);
	mainLayout->addWidget(this->BaselineBox,2,3);	
	#endif
	mainLayout->addWidget(OkButton,3,2);
	mainLayout->addWidget(CancelButton,3,3);
	
	setLayout(mainLayout);
	setWindowTitle(tr("Save ScreenShot"));
}
QPushButton *ScreenShotDialog::createButton(const QString &text, const char *member)
{
	QPushButton *button = new QPushButton(text);
	connect(button, SIGNAL(clicked()), this, member);
	return button;
}
QComboBox *ScreenShotDialog::createComboBox(const QString &text)
{
	QComboBox *comboBox = new QComboBox;
	comboBox->setEditable(true);
	comboBox->addItem(text);
	comboBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
	return comboBox;
}
// Specify folder to save file
void ScreenShotDialog::Browse()
{
	curdirectory = 	QFileDialog::getExistingDirectory(this, tr("Directory"), 
															QFileInfo(imageDir).dir().canonicalPath());
	if (!curdirectory.isEmpty())
	{
		if (directoryComboBox->findText(curdirectory) == -1)
		{
			directoryComboBox->addItem(curdirectory);
		}
		directoryComboBox->setCurrentIndex(directoryComboBox->findText(curdirectory));
	}
}
void ScreenShotDialog::save()
{
	//check if directory exist, otherwise make directory
	curdirectory = directoryComboBox->currentText();
	QDir directory(curdirectory);
	if(!directory.exists())
	{
		directory.mkdir(curdirectory);
	}
	this->magnifynum = ZoomSpinBox->value();

	QDialog::accept(); //save and close export cell dialog
}
// Return the directory, filename, and magnification factor
QString ScreenShotDialog::getDir()
{
	return curdirectory;
}
QString ScreenShotDialog::getfileName()
{
	fileName = fileNameLine->text();
	return fileName;
}
int ScreenShotDialog::getMagnification()
{
	return magnifynum;
}
bool ScreenShotDialog::getSave()
{
	return saveclicked;
}

bool ScreenShotDialog::getBaseline()
{
  #ifdef USE_QT_TESTING
  return this->BaselineBox->isChecked();
  #endif
  return false;
}