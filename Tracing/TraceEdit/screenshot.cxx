/* Screenshot dialog for saving the screenshot and allowing for customizable features such as zooming.
*/

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

	QGridLayout *mainLayout = new QGridLayout;
	mainLayout->addWidget(directoryLabel,0,0);
	mainLayout->addWidget(directoryComboBox,0,1,1,2);
	mainLayout->addWidget(browseButton,0,3);
	mainLayout->addWidget(fileNameLabel,1,0);
	mainLayout->addWidget(fileNameLine,1,1,1,3);
	mainLayout->addWidget(ZoomLabel,2,0);
	mainLayout->addWidget(ZoomSpinBox,2,1);
	mainLayout->addWidget(OkButton,3,2);
	mainLayout->addWidget(CancelButton,3,3);
	
	setLayout(mainLayout);
	setWindowTitle(tr("Save ScreenShot"));
}
void ScreenShotDialog::Browse()
{
	//std::cout << curdirectoryswc.toStdString() << std::endl;
	//std::cout << curdirectoryjpg.toStdString() << std::endl;

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