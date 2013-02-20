/*****************************************************************************************
// File dialog for autocellexport to save individual traces in swc and jpg files.		//
// You are given the options of saving swc files, jpg files, or both in the chosen		//
// directory you choose.  You can also name the files as you wish or use default names.	//
*****************************************************************************************/

#include "cellexport.h"

SaveCellExportDialog::SaveCellExportDialog(QWidget* parent, QString curdirectoryswc, QString curdirectoryjpg, QString swcfileName, QString jpgfileName, bool changeswcfileName, bool changejpgfileName)
: QDialog(parent)
{
	this->curdirectoryswc = curdirectoryswc;
	this->curdirectoryjpg = curdirectoryjpg;
	saveclicked = false;

// SWC files directory setup
	QLabel *swclabel = new QLabel(tr("Directory for SWC files:"));
	swcdirectoryComboBox = createComboBox(curdirectoryswc);
	swclabel->setBuddy(swcdirectoryComboBox);
	swcbrowseButton = createButton(tr("&Browse..."), SLOT(swcBrowse()));
	swcmoreButton = new QPushButton(tr("&More..."));
	swcmoreButton->setCheckable(true);

// JPG files directory setup
	QLabel *jpglabel = new QLabel(tr("Directory for JPG files:"));
	jpgdirectoryComboBox = createComboBox(curdirectoryjpg);
	jpglabel->setBuddy(jpgdirectoryComboBox);
	jpgbrowseButton = createButton(tr("&Browse..."), SLOT(jpgBrowse()));
	jpgmoreButton = new QPushButton(tr("&More..."));
	jpgmoreButton->setCheckable(true);

// SWC files: customize naming files setup (original filename, "cell_1", or new name with numbers)
	swcextension = new QWidget;
	originalswcfileNameButton = new QRadioButton(tr("Keep original filename"), swcextension);
	originalswcfileNameButton->setChecked(true);
	renumberswcfileNameButton = new QRadioButton(tr("Label files by xyz: (ex. cell_20_32_3)"), swcextension);
	renameswcfileNameButton = new QRadioButton(tr("Rename and auto-assign number to swc files"), swcextension);
	QLabel *swcLineEditLabel = new QLabel(tr("Custom name:"));
	nameswcfileNameLine = new QLineEdit(swcextension);
	connect(swcmoreButton, SIGNAL(toggled(bool)), swcextension, SLOT(setVisible(bool)));
	connect(nameswcfileNameLine, SIGNAL(textChanged(const QString &)), this, SLOT(swcfilenaming()));

// JPG files: customize naming files setup (original filename, "cell_1", or new name with numbers)
	jpgextension = new QWidget;
	originaljpgfileNameButton = new QRadioButton(tr("Keep original filename"), jpgextension);
	originaljpgfileNameButton->setChecked(true);
	renumberjpgfileNameButton = new QRadioButton(tr("Label files by xyz: (ex. cell_20_32_3)"), jpgextension);
	renamejpgfileNameButton = new QRadioButton(tr("Rename and auto-assign number to jpg files"), jpgextension);
	QLabel *jpgLineEditLabel = new QLabel(tr("Custom name:"));
	namejpgfileNameLine = new QLineEdit(jpgextension);
	connect(jpgmoreButton, SIGNAL(toggled(bool)), jpgextension, SLOT(setVisible(bool)));
	connect(namejpgfileNameLine, SIGNAL(textChanged(const QString &)), this, SLOT(jpgfilenaming()));
	
// Create layout
	QHBoxLayout *swcdirLayout = new QHBoxLayout();
	swcdirLayout->addWidget(swclabel);
	swcdirLayout->addWidget(swcdirectoryComboBox);

	QHBoxLayout *jpgdirLayout = new QHBoxLayout();
	jpgdirLayout->addWidget(jpglabel);
	jpgdirLayout->addWidget(jpgdirectoryComboBox);

	QVBoxLayout* swcextensionLayout = new QVBoxLayout();
	swcextensionLayout->addWidget(originalswcfileNameButton);
	swcextensionLayout->addWidget(renumberswcfileNameButton);
	swcextensionLayout->addWidget(renameswcfileNameButton);
	swcextensionLayout->addWidget(nameswcfileNameLine);
	swcextension->setLayout(swcextensionLayout);

	QVBoxLayout* jpgextensionLayout = new QVBoxLayout();
	jpgextensionLayout->addWidget(originaljpgfileNameButton);
	jpgextensionLayout->addWidget(renumberjpgfileNameButton);
	jpgextensionLayout->addWidget(renamejpgfileNameButton);
	jpgextensionLayout->addWidget(namejpgfileNameLine);
	jpgextension->setLayout(jpgextensionLayout);

	OkButton = new QPushButton(tr("Ok"));
	OkButton->setDefault(true);
	connect(OkButton, SIGNAL(clicked()), this, SLOT(save()));
	
	CancelButton = new QPushButton(tr("Cancel"));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(close()));

// Different layout
	//QVBoxLayout *leftswcLayout = new QVBoxLayout();
	//leftswcLayout->addLayout(swcdirLayout);
	//leftswcLayout->addWidget(swcextension);
	//leftswcLayout->addStretch();

	//QVBoxLayout *swcButtonLayout = new QVBoxLayout();
	//swcButtonLayout->addWidget(swcbrowseButton);
	//swcButtonLayout->addWidget(swcmoreButton);
	//swcButtonLayout->addStretch();

	//QHBoxLayout *fullswcLayout = new QHBoxLayout();
	//fullswcLayout->addLayout(leftswcLayout);
	//fullswcLayout->addLayout(swcButtonLayout);


	//QVBoxLayout *leftjpgLayout = new QVBoxLayout();
	//leftjpgLayout->addLayout(jpgdirLayout);
	//leftjpgLayout->addWidget(jpgextension);
	//leftjpgLayout->addStretch();

	//QVBoxLayout *jpgButtonLayout = new QVBoxLayout();
	//jpgButtonLayout->addWidget(jpgbrowseButton);
	//jpgButtonLayout->addWidget(jpgmoreButton);
	//jpgButtonLayout->addStretch();

	//QHBoxLayout *fulljpgLayout = new QHBoxLayout();
	//fulljpgLayout->addLayout(leftjpgLayout);
	//fulljpgLayout->addLayout(jpgButtonLayout);

	QHBoxLayout *bottomLayout = new QHBoxLayout();
	bottomLayout->addWidget(OkButton);
	bottomLayout->addWidget(CancelButton);

	//QVBoxLayout *MainLayout = new QVBoxLayout();
	//MainLayout->addLayout(fullswcLayout);
	//MainLayout->addLayout(fulljpgLayout);
	//MainLayout->addLayout(bottomLayout);
	//MainLayout->setSizeConstraint(QLayout::SetFixedSize);

	// QGridLayout - the 1st number: initial row index, the 2nd number: initial column index, 
	//					the 3rd number: rowSpan, and the 4th number: columnSpan.
	saveSWCGroupBox = new QGroupBox(tr("Save SWC files"));
	saveSWCGroupBox->setCheckable(true);
	QGridLayout *swcLayout = new QGridLayout(saveSWCGroupBox);
	swcLayout->addLayout(swcdirLayout,0,0);
	swcLayout->addWidget(swcbrowseButton,0,2);
	swcLayout->addWidget(swcmoreButton,1,2,Qt::AlignTop);
	swcLayout->addWidget(swcextension,1,0);

	saveJPGGroupBox = new QGroupBox(tr("Save JPG files"));
	saveJPGGroupBox->setCheckable(true);
	QGridLayout *jpgLayout = new QGridLayout(saveJPGGroupBox);
	jpgLayout->addLayout(jpgdirLayout,0,0);
	jpgLayout->addWidget(jpgbrowseButton,0,2);
	jpgLayout->addWidget(jpgmoreButton,1,2,Qt::AlignTop);
	jpgLayout->addWidget(jpgextension,1,0);

	QGridLayout *cellexportLayout = new QGridLayout;
	cellexportLayout->addWidget(saveSWCGroupBox,0,0);
	cellexportLayout->addWidget(saveJPGGroupBox,1,0);
	cellexportLayout->addLayout(bottomLayout,5,0,Qt::AlignRight);
	//cellexportLayout->addWidget(OkButton,5,0);
	//cellexportLayout->addWidget(CancelButton,5,1);
	//cellexportLayout->setSizeConstraint(QLayout::SetFixedSize);
	
	setLayout(cellexportLayout);

	setWindowTitle(tr("Export Cells"));
	resize(500,250);

	swcextension->hide();
	jpgextension->hide();
}
// Specify folder to save swc file
void SaveCellExportDialog::swcBrowse()
{
	curdirectoryswc = QFileDialog::getExistingDirectory(this, tr("Choose Directory for SWC files"), 
															QFileInfo(curdirectoryswc).dir().canonicalPath());
	if (!curdirectoryswc.isEmpty())
	{
		if (swcdirectoryComboBox->findText(curdirectoryswc) == -1)
		{
			swcdirectoryComboBox->addItem(curdirectoryswc);
		}
		swcdirectoryComboBox->setCurrentIndex(swcdirectoryComboBox->findText(curdirectoryswc));	
	}
}
// Specify folder to save jpg file
void SaveCellExportDialog::jpgBrowse()
{
	curdirectoryjpg = QFileDialog::getExistingDirectory(this, tr("Choose Directory for JPG files"), 
															QFileInfo(curdirectoryjpg).dir().canonicalPath());
	if (!curdirectoryjpg.isEmpty()) 
	{
		if (jpgdirectoryComboBox->findText(curdirectoryjpg) == -1)
		{
			jpgdirectoryComboBox->addItem(curdirectoryjpg);
		}
		jpgdirectoryComboBox->setCurrentIndex(jpgdirectoryComboBox->findText(curdirectoryjpg));
	}
}
//checkmark rename button if rename text is changed
void SaveCellExportDialog::swcfilenaming()
{
	renameswcfileNameButton->setChecked(true);
}
void SaveCellExportDialog::jpgfilenaming()
{
	renamejpgfileNameButton->setChecked(true);
}

void SaveCellExportDialog::save()
{
	saveclicked = true;
	//check if directory exist, otherwise make directory
	curdirectoryswc = swcdirectoryComboBox->currentText();
	QDir SWCdirectory(curdirectoryswc);
	if(!SWCdirectory.exists())
	{
		SWCdirectory.mkdir(curdirectoryswc);
	}
	curdirectoryjpg = jpgdirectoryComboBox->currentText();
	QDir JPGdirectory(curdirectoryjpg);
	if(!JPGdirectory.exists())
	{
		JPGdirectory.mkdir(curdirectoryjpg);
	}

	QDialog::accept(); //save and close export cell dialog
}
QPushButton *SaveCellExportDialog::createButton(const QString &text, const char *member)
{
	QPushButton *button = new QPushButton(text);
	connect(button, SIGNAL(clicked()), this, member);
	return button;
}
QComboBox *SaveCellExportDialog::createComboBox(const QString &text)
{
	QComboBox *comboBox = new QComboBox;
	comboBox->setEditable(true);
	comboBox->addItem(text);
	comboBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
	return comboBox;
}
// return directory to save files
QString SaveCellExportDialog::getSWCDir()
{
	if (!saveSWCGroupBox->isChecked())
	{
		curdirectoryswc.clear();
	}
	return curdirectoryswc;
}
QString SaveCellExportDialog::getJPGDir()
{
	if (!saveJPGGroupBox->isChecked())
	{
		curdirectoryjpg.clear();
	}
	return curdirectoryjpg;
}
//customize the swc filename
QString SaveCellExportDialog::getSWCfileName()
{
	if (renameswcfileNameButton->isChecked())
	{
		swcfileName = nameswcfileNameLine->text();
	}
	else
	{
		swcfileName.clear();
	}
	return swcfileName;
}
//customize the jpg filename
QString SaveCellExportDialog::getJPGfileName()
{
	if (renamejpgfileNameButton->isChecked())
	{
		jpgfileName = namejpgfileNameLine->text();
	}
	else
	{
		jpgfileName.clear();
	}
	return jpgfileName;
}
//decide whether to keep original filename or change it
bool SaveCellExportDialog::differentSWCfileName()
{
	if (originalswcfileNameButton->isChecked())
	{
		changeswcfileName = false;
	}
	else
	{
		changeswcfileName = true;
	}
	return changeswcfileName;
}
bool SaveCellExportDialog::differentJPGfileName()
{
	if (originaljpgfileNameButton->isChecked())
	{
		changejpgfileName = false;
	}
	else
	{
		changejpgfileName = true;
	}
	return changejpgfileName;
}
bool SaveCellExportDialog::getSave()
{
	return saveclicked;
}