#include "TrackingDialog.h"

TrackingDialog::TrackingDialog(QWidget *parent)
: QDialog(parent)
{
	 tabWidget = new QTabWidget;
	 paramtertab = new ParametersTab();
	 foldertab = new FoldersTab();


     tabWidget->addTab(paramtertab, tr("Parameters"));
     tabWidget->addTab(foldertab, tr("Folders"));

     trackdlgButton = new QDialogButtonBox(QDialogButtonBox::Ok
                                      | QDialogButtonBox::Cancel);
     connect(trackdlgButton, SIGNAL(accepted()), this, SLOT(accept()));
     connect(trackdlgButton, SIGNAL(rejected()), this, SLOT(reject()));

     QVBoxLayout *mainLayout = new QVBoxLayout;
     mainLayout->addWidget(tabWidget);
     mainLayout->addWidget(trackdlgButton);
     setLayout(mainLayout);
	 this->setWindowTitle(tr("Tracking"));
	 this->setFixedSize(400,450);
}
ParametersTab::ParametersTab(QWidget *parent)
{
	QGroupBox * volumespacingBox = new QGroupBox(tr("Spacing/Volume"));
	QGridLayout * volumespacingLayout = new QGridLayout;

	xspacingLabel = new QLabel(tr("X Spacing"));
	xspacingSpinBox = new QDoubleSpinBox;
	xspacingSpinBox->setValue(0.64);
	xspacingLabel->setBuddy(xspacingSpinBox);

	yspacingLabel = new QLabel(tr("Y Spacing"));
	yspacingSpinBox = new QDoubleSpinBox;
	yspacingSpinBox->setValue(0.64);
	yspacingLabel->setBuddy(yspacingSpinBox);


	zspacingLabel = new QLabel(tr("Z Spacing"));
	zspacingSpinBox = new QDoubleSpinBox;
	zspacingSpinBox->setValue(2.0);
	zspacingLabel->setBuddy(zspacingSpinBox);

	volumevarLabel = new QLabel(tr("Volume Variance"));
	volumevarSpinBox = new QDoubleSpinBox;
	volumevarSpinBox->setRange(0.0,100000.0);
	volumevarSpinBox->setValue(90000.0);
	volumevarLabel->setBuddy(volumevarSpinBox);

	volumespacingLayout->addWidget(xspacingLabel, 0, 0);
	volumespacingLayout->addWidget(xspacingSpinBox, 0, 1);
	volumespacingLayout->addWidget(yspacingLabel, 0, 2 );
	volumespacingLayout->addWidget(yspacingSpinBox, 0, 3);
	volumespacingLayout->addWidget(zspacingLabel, 1, 0);
	volumespacingLayout->addWidget(zspacingSpinBox, 1, 1);
	volumespacingLayout->addWidget(volumevarLabel, 1, 2);
	volumespacingLayout->addWidget(volumevarSpinBox, 1, 3);
	volumespacingBox->setLayout(volumespacingLayout);

	QGroupBox * meanBox = new QGroupBox(tr("Mean"));
	QGridLayout * meanLayout = new QGridLayout;

	distmeanLabel = new QLabel(tr("Distance"));
	distmeanSpinBox = new QDoubleSpinBox;
	distmeanSpinBox->setValue(2.0);
	distmeanLabel->setBuddy(distmeanSpinBox);

	timemeanLabel = new QLabel(tr("Time"));
	timemeanSpinBox = new QDoubleSpinBox;
	timemeanSpinBox->setValue(0.0);
	timemeanLabel->setBuddy(timemeanSpinBox);

	overlapmeanLabel = new QLabel(tr("Overlap"));
	overlapmeanSpinBox = new QDoubleSpinBox;
	overlapmeanSpinBox->setValue(0.0);
	overlapmeanLabel->setBuddy(overlapmeanSpinBox);

	boundarydistmeanLabel = new QLabel(tr("Bound Distance"));
	boundarydistmeanSpinBox = new QDoubleSpinBox;
	boundarydistmeanSpinBox->setValue(2.0);
	boundarydistmeanLabel->setBuddy(boundarydistmeanSpinBox);

	meanLayout->addWidget(distmeanLabel, 0, 0);
	meanLayout->addWidget(distmeanSpinBox, 0, 1);
	meanLayout->addWidget(timemeanLabel, 1, 0 );
	meanLayout->addWidget(timemeanSpinBox, 1, 1);
	meanLayout->addWidget(overlapmeanLabel, 2, 0);
	meanLayout->addWidget(overlapmeanSpinBox, 2, 1);
	meanLayout->addWidget(boundarydistmeanLabel, 3, 0);
	meanLayout->addWidget(boundarydistmeanSpinBox, 3, 1);
	meanBox->setLayout(meanLayout);

	QGroupBox * varBox = new QGroupBox(tr("Variance"));
	QGridLayout * varLayout = new QGridLayout;

	distvarLabel = new QLabel(tr("Distance"));
	distvarSpinBox = new QDoubleSpinBox;
	distvarSpinBox->setValue(25.0);
	distvarLabel->setBuddy(distvarSpinBox);

	timevarLabel = new QLabel(tr("Time"));
	timevarSpinBox = new QDoubleSpinBox;
	timevarSpinBox->setValue(0.1);
	timevarLabel->setBuddy(timevarSpinBox);

	overlapvarLabel = new QLabel(tr("Overlap"));
	overlapvarSpinBox = new QDoubleSpinBox;
	overlapvarSpinBox->setValue(1.0);
	overlapvarLabel->setBuddy(overlapvarSpinBox);

	boundarydistvarLabel = new QLabel(tr("Bound Distance"));
	boundarydistvarSpinBox = new QDoubleSpinBox;
	boundarydistvarSpinBox->setValue(12.0);
	boundarydistvarLabel->setBuddy(boundarydistvarSpinBox);

	varLayout->addWidget(distvarLabel, 0, 0);
	varLayout->addWidget(distvarSpinBox, 0, 1);
	varLayout->addWidget(timevarLabel, 1, 0 );
	varLayout->addWidget(timevarSpinBox, 1, 1);
	varLayout->addWidget(overlapvarLabel, 2, 0);
	varLayout->addWidget(overlapvarSpinBox, 2, 1);
	varLayout->addWidget(boundarydistvarLabel, 3, 0);
	varLayout->addWidget(boundarydistvarSpinBox, 3, 1);
	varBox->setLayout(varLayout);

	QGroupBox * probBox = new QGroupBox(tr("Prior Probabilities"));
	QGridLayout * probLayout = new QGridLayout;

	mergesplitpriorLabel = new QLabel(tr("Merge/Split"));
	mergesplitpriorSpinBox = new QDoubleSpinBox;
	mergesplitpriorSpinBox->setValue(1.0);
	mergesplitpriorSpinBox->setRange(0.0,1.0);
	mergesplitpriorLabel->setBuddy(mergesplitpriorSpinBox);

	appeardisappearpriorLabel = new QLabel(tr("Appear/Disappear"));
	appeardisappearpriorSpinBox = new QDoubleSpinBox;
	appeardisappearpriorSpinBox->setValue(1.0);
	appeardisappearpriorSpinBox->setRange(0.0,1.0);
	appeardisappearpriorLabel->setBuddy(appeardisappearpriorSpinBox);


	translationpriorLabel = new QLabel(tr("Translation"));
	translationpriorSpinBox = new QDoubleSpinBox;
	translationpriorSpinBox->setValue(1.0);
	translationpriorSpinBox->setRange(0.0,1.0);
	translationpriorLabel->setBuddy(translationpriorSpinBox);

	probLayout->addWidget(mergesplitpriorLabel, 0, 0);
	probLayout->addWidget(mergesplitpriorSpinBox, 0, 1);
	probLayout->addWidget(appeardisappearpriorLabel, 0, 2 );
	probLayout->addWidget(appeardisappearpriorSpinBox, 0, 3);
	probLayout->addWidget(translationpriorLabel, 0, 4);
	probLayout->addWidget(translationpriorSpinBox, 0, 5);
	probBox->setLayout(probLayout);


    QGridLayout * mainLayout = new QGridLayout;
	mainLayout->addWidget(volumespacingBox, 0, 0, 1, 2);
	mainLayout->addWidget(probBox, 1, 0, 1, 2);
	mainLayout->addWidget(meanBox, 2, 0);
	mainLayout->addWidget(varBox, 2, 1);
	setLayout(mainLayout);

}
FoldersTab::FoldersTab(QWidget *parent)
{
	lastdirpath = QDir::homePath();
	lastdirpath = "C:\\Lab\\ArunFiles\\Data\\Tracking";

	browseNumbersDirectoryButton = createButton(tr("&Browse..."), SLOT(browseNumbersDirectory()));
	browseEntropyDirectoryButton = createButton(tr("&Browse..."), SLOT(browseEntropyDirectory()));
	browseDebugDirectoryButton = createButton(tr("&Browse..."), SLOT(browseDebugDirectory()));
	browseResultDirectoryButton	= createButton(tr("&Browse..."), SLOT(browseResultDirectory()));

	numbersdirectoryLabel = new QLabel(tr("Numbers File:"));
	entropyfileLabel = new QLabel(tr("Entropy Filename:"));
	entropydirectoryLabel = new QLabel(tr("Entropy Folder:"));
	debugdprefixLabel = new QLabel(tr("Debug Prefix:"));
	debugfileLabel = new QLabel(tr("Debug Filename:"));
	debugdirectoryLabel = new QLabel(tr("Debug Folder:"));
	loadLabel = new QLabel(tr("Load Files:"));
	saveLabel = new QLabel(tr("Save Files:"));
	resultprefixLabel = new QLabel(tr("Result Prefix:"));
	resultdirectoryLabel = new QLabel(tr("Result Folder:"));


	numbersdirectoryLineEdit = createLineEdit("C:\\Lab\\ArunFiles\\Data\\Tracking\\numbers.bmp");
	entropyfileLineEdit = createLineEdit("Entropy.txt");
	entropydirectoryLineEdit = createLineEdit(lastdirpath);
	debugfileLineEdit = createLineEdit("debug.txt");
	//debugdirectoryLineEdit = createLineEdit(lastdirpath);
	debugdirectoryLineEdit = createLineEdit("C:\\Lab\\ArunFiles\\Data\\Tracking\\debug\\");
	
	debugdprefixLineEdit = createLineEdit("1color");
	resultprefixLineEdit= createLineEdit("labeled_tracks");
	resultdirectoryLineEdit = createLineEdit(lastdirpath);

	QGroupBox * loadgroupBox = new QGroupBox(tr("Load Files"));
	QGridLayout * loadLayout = new QGridLayout;
	loadLayout->addWidget(numbersdirectoryLabel, 0, 0);
	loadLayout->addWidget(numbersdirectoryLineEdit, 0, 1);
	loadLayout->addWidget(browseNumbersDirectoryButton, 0, 2);
	loadgroupBox->setLayout(loadLayout);

	QGroupBox * savegroupBox = new QGroupBox(tr("Save Files"));
	QGridLayout * saveLayout = new QGridLayout;

	QGroupBox * reslutgroupBox = new QGroupBox;
	QGridLayout * reslutLayout = new QGridLayout;
	reslutLayout->addWidget(resultprefixLabel, 2, 0);
	reslutLayout->addWidget(resultprefixLineEdit, 2, 1, 1, 2);
	reslutLayout->addWidget(resultdirectoryLabel, 3, 0);
	reslutLayout->addWidget(resultdirectoryLineEdit, 3, 1);
	reslutLayout->addWidget(browseResultDirectoryButton, 3, 2);
	reslutgroupBox->setLayout(reslutLayout);

	QGroupBox * entropygroupBox = new QGroupBox;
	QGridLayout * entropyLayout = new QGridLayout;
	entropyLayout->addWidget(entropyfileLabel, 2, 0);
	entropyLayout->addWidget(entropyfileLineEdit, 2, 1, 1, 2);
	entropyLayout->addWidget(entropydirectoryLabel, 3, 0);
	entropyLayout->addWidget(entropydirectoryLineEdit, 3, 1);
	entropyLayout->addWidget(browseEntropyDirectoryButton, 3, 2);
	entropygroupBox->setLayout(entropyLayout);

	QGroupBox * debuggroupBox = new QGroupBox;
	QGridLayout * debugLayout = new QGridLayout;
	debugLayout->addWidget(debugdprefixLabel, 0, 0);
	debugLayout->addWidget(debugdprefixLineEdit, 0, 1, 1, 2);
	debugLayout->addWidget(debugfileLabel, 1, 0);
	debugLayout->addWidget(debugfileLineEdit, 1, 1, 1, 2);
	debugLayout->addWidget(debugdirectoryLabel, 2, 0);
	debugLayout->addWidget(debugdirectoryLineEdit, 2, 1);
	debugLayout->addWidget(browseDebugDirectoryButton, 2, 2);
	debuggroupBox->setLayout(debugLayout);

	saveLayout->addWidget(reslutgroupBox,0,0);
	saveLayout->addWidget(entropygroupBox,1,0);
	saveLayout->addWidget(debuggroupBox,2,0);
	savegroupBox->setLayout(saveLayout);

	QVBoxLayout * mainLayout = new QVBoxLayout;
	mainLayout->addWidget(loadgroupBox,0,0);
	mainLayout->addWidget(savegroupBox,1,0);
	setLayout(mainLayout);

}
QPushButton * FoldersTab::createButton(const QString &text, const char *member)
{
	 QPushButton * button = new QPushButton(text);
	 connect(button, SIGNAL(clicked()), this, member);
	 return button;
}
QLineEdit * FoldersTab::createLineEdit(const QString &text)
 {
     QLineEdit * lineEdit = new QLineEdit;
     lineEdit->setText(text);
     lineEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
     return lineEdit;
 }
void FoldersTab::browseNumbersDirectory()
{
	QString numbersfileName = QFileDialog::getOpenFileName(this,
									tr("Load File"),lastdirpath,
									tr("Images (*.bmp)"));
	if (numbersfileName.isEmpty()) return;
	numbersdirectoryLineEdit->setText(numbersfileName);

}
void FoldersTab::browseEntropyDirectory()
{
	QString entropyDirectoryName = QFileDialog::getExistingDirectory(this,
									tr("Select Entropy Folder"),lastdirpath);
	if (entropyDirectoryName.isEmpty()) return;
	entropydirectoryLineEdit->setText(entropyDirectoryName);
	lastdirpath.clear();
	lastdirpath.push_back(entropydirectoryLineEdit->text());
}
void FoldersTab::browseDebugDirectory()
{
	QString debugDirectoryName = QFileDialog::getExistingDirectory(this,
									tr("Select Debug Folder"), lastdirpath);
	if (debugDirectoryName.isEmpty()) return;
	debugdirectoryLineEdit->setText(debugDirectoryName);
	lastdirpath.clear();
	lastdirpath.push_back(debugdirectoryLineEdit->text());
}

void FoldersTab::browseResultDirectory()
{
	QString resultDirectoryName = QFileDialog::getExistingDirectory(this,
									tr("Select Reslut Folder"), lastdirpath);
	if (resultDirectoryName.isEmpty()) return;
	resultdirectoryLineEdit->setText(resultDirectoryName);
	lastdirpath.clear();
	lastdirpath.push_back(resultdirectoryLineEdit->text());
}

std::vector<std::pair<std::string,std::string> > TrackingDialog::getFolders(void)
{
	std::vector<std::pair<std::string,std::string> > filefolders;
	std::pair<std::string,std::string> tmpfilefolders;

	tmpfilefolders.first = "numbers file:";
	tmpfilefolders.second = foldertab->numbersdirectoryLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "entropy filename:";
	tmpfilefolders.second = foldertab->entropyfileLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "entropy file directory:";
	tmpfilefolders.second = foldertab->entropydirectoryLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "debug prefix:";
	tmpfilefolders.second = foldertab->debugdprefixLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "debug filename:";
	tmpfilefolders.second = foldertab->debugfileLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "debug file directory:";
	tmpfilefolders.second = foldertab->debugdirectoryLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "result filename:";
	tmpfilefolders.second = foldertab->resultprefixLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);

	tmpfilefolders.first = "result file directory:";
	tmpfilefolders.second = foldertab->resultdirectoryLineEdit->text().toStdString();
	filefolders.push_back(tmpfilefolders);	

	return filefolders;
}



std::vector<std::pair<std::string,float> > TrackingDialog::getParameters()
{

	 std::vector<std::pair<std::string,float> > trackparameters;
	 std::pair<std::string,float> tmppair;
	 

	tmppair.first =  "xspacing";
	tmppair.second = static_cast<float> (paramtertab->xspacingSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "yspacing";
	tmppair.second = static_cast<float> (paramtertab->yspacingSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "zspacing";
	tmppair.second = static_cast<float> (paramtertab->zspacingSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "distvar";
	tmppair.second = static_cast<float> (paramtertab->distvarSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "distmean";
	tmppair.second = static_cast<float> (paramtertab->distmeanSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "timevar";
	tmppair.second = static_cast<float> (paramtertab->timevarSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first =  "timemean";
	tmppair.second = static_cast<float> (paramtertab->timemeanSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first = "overlapvar";
	tmppair.second = static_cast<float> (paramtertab->overlapvarSpinBox->value());
	trackparameters.push_back(tmppair);

	tmppair.first = "overlapmean";
	tmppair.second = static_cast<float> (paramtertab->overlapmeanSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "volumevar";
	tmppair.second = static_cast<float> (paramtertab->volumevarSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "mergesplitprior";
	tmppair.second = static_cast<float> (paramtertab->mergesplitpriorSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "appeardisappearprior";
	tmppair.second = static_cast<float> (paramtertab->appeardisappearpriorSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "translationprior";
	tmppair.second = static_cast<float> (paramtertab->translationpriorSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "boundarydistmean";
	tmppair.second = static_cast<float> (paramtertab->boundarydistmeanSpinBox->value());
	trackparameters.push_back(tmppair);
	
	tmppair.first = "boundarydistvar";
	tmppair.second = static_cast<float> (paramtertab->boundarydistvarSpinBox->value());
	trackparameters.push_back(tmppair);

	return trackparameters;

}

