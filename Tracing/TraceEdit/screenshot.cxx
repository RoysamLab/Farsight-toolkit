/* Screenshot dialog for saving the screenshot and allowing for customizable features such as zooming.
*/

#include "screenshot.h"

ScreenShotDialog::ScreenShotDialog(QWidget* parent, QString fileName, QString imageDir)
: QDialog(parent)
{
	//QLabel *directoryLabel = new QLabel(tr("Directory:"));
	directoryComboBox = createComboBox();	
	//directoryLabel->setBuddy(directoryComboBox);
	browseButton = createButton(tr("&Browse..."), SLOT(Browse()));

	//QLabel *fileNameLabel = new QLabel(tr("Filename:"));
	fileNameLine = new QLineEdit();
	//fileNameLabel->setBuddy(fileNameLine);

	//QLabel * ZoomLabel = new QLabel(tr("Zoom: "));
	ZoomSpinBox = new QSpinBox();
	//ZoomLabel->setBuddy(ZoomSpinBox);

	OkButton = new QPushButton(tr("Ok"));
	OkButton->setDefault(true);
	connect(OkButton, SIGNAL(clicked()), this, SLOT(save()));
	
	CancelButton = new QPushButton(tr("Cancel"));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(close()));

	QHBoxLayout * browseLayout = new QHBoxLayout();
	browseLayout->addWidget(directoryComboBox);
	browseLayout->addWidget(browseButton);

	QHBoxLayout * ConcludeLayout = new QHBoxLayout();
	ConcludeLayout->addWidget(OkButton);
	ConcludeLayout->addWidget(CancelButton);

	//QVBoxLayout * verticalLayout = new QVBoxLayout();
	//verticalLayout->addWidget(directoryComboBox);
	//verticalLayout->addWidget(fileNameLine);
	//verticalLayout->addWidget(ZoomSpinBox);
	//verticalLayout->addWidget(OkButton);
	//verticalLayout->addWidget(CancelButton);

	//QGroupBox *displaySettings = new QGroupBox();
	QFormLayout *DisplayLayout = new QFormLayout();
	DisplayLayout->addRow(tr("Browse:"),browseLayout);
	DisplayLayout->addRow(tr("Filename:"),fileNameLine);
	DisplayLayout->addRow(tr("Zoom:"),ZoomSpinBox);

	QVBoxLayout * mainLayout = new QVBoxLayout();
	mainLayout->addLayout(DisplayLayout);
	mainLayout->addLayout(ConcludeLayout);
	
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