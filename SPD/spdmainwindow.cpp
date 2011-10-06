#include "spdmainwindow.h"
#include <QGridLayout>
#include <QLCDNumber>
#include <QFileDialog>
#include <QDir>
#include <fstream>
#include <QMessageBox>
//define NDEBUG
#include <assert.h>

using std::ifstream;
using std::endl;

SPDMainWindow::SPDMainWindow(QWidget *parent) :
    QWidget(parent)
{
	SPDModel = NULL;

    dataFileLabel = new QLabel(tr("Choose file:"));

    int frameStyle = QFrame::Sunken | QFrame::Panel;
    dataFileName = new QLabel;
    dataFileName->setFrameStyle(frameStyle);

    browseButton = new QPushButton(tr("Browse"));
    loadButton = new QPushButton(tr("Load"));

    featureNumLabel = new QLabel(tr("Feature size:"));
    featureNum = new QLabel;
    featureNum->setFrameStyle(frameStyle);
    sampleNumLabel = new QLabel(tr("Sample size:"));
    sampleNum = new QLabel;
    sampleNum->setFrameStyle(frameStyle);

    clusterCoherenceLabel = new QLabel(tr("Coherence:"));
    clusterCoherenceBox = new QLineEdit;
    clusterMergeLabel = new QLabel(tr("Merge Coherence:"));
    clusterMergeBox = new QLineEdit;
    clusterButton = new QPushButton(tr("Agglomerate"));
    clusterResultButton = new QPushButton(tr("Show Results"));

	generateMSTButton = new QPushButton(tr("MST"));
	mstState = new QLabel();
	showMSTButton = new QPushButton(tr("Show MST"));

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
    connect(clusterButton, SIGNAL(clicked()), this, SLOT(clusterFunction()));
    connect(clusterResultButton, SIGNAL(clicked()), this, SLOT(showResult()));
	connect(generateMSTButton, SIGNAL(clicked()), this, SLOT(generateMST()));
	connect(showMSTButton, SIGNAL(clicked()), this, SLOT(showMST()));

    QGridLayout *mainLayout = new QGridLayout;

    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }

    for ( int row = 1; row <= 7; row++)
    {
        mainLayout->setRowMinimumHeight(row,20);
        mainLayout->setRowStretch(row, 1);
    }

    mainLayout->addWidget(dataFileLabel, 0, 0);

    mainLayout->addWidget(dataFileName, 1, 0, 1, 2);
    mainLayout->addWidget(browseButton, 1, 2);
    mainLayout->addWidget(loadButton, 2, 2);

    mainLayout->addWidget(featureNumLabel, 2, 0);
    mainLayout->addWidget(featureNum, 2, 1);
    mainLayout->addWidget(clusterCoherenceLabel, 4, 0);
    mainLayout->addWidget(clusterCoherenceBox, 4, 1);

    mainLayout->addWidget(sampleNumLabel, 3, 0);
    mainLayout->addWidget(sampleNum, 3, 1);
    mainLayout->addWidget(clusterMergeLabel, 5, 0);
    mainLayout->addWidget(clusterMergeBox, 5, 1);

    mainLayout->addWidget(clusterButton, 6, 0);
    mainLayout->addWidget(clusterResultButton, 6, 1);

    mainLayout->addWidget(generateMSTButton, 7, 0);
    mainLayout->addWidget(mstState, 7, 2);
	mainLayout->addWidget(showMSTButton, 7, 1);

    setLayout(mainLayout);

	SPDModel = SPDAnalysisModel::InitInstance();

	assert(SPDModel!=NULL);
}

SPDMainWindow::~SPDMainWindow()
{
	SPDAnalysisModel::DeInstance();
}

void SPDMainWindow::browse()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileName->setText(fileName);
        this->FileName = fileName;
    }
}

void SPDMainWindow::load()
{
	std::string file = this->FileName.toStdString();

	if ( true == this->SPDModel->ReadCellTraceFile(file.c_str()))
	{
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		this->SPDModel->NormalizeData();
	}
}

void SPDMainWindow::clusterFunction()
{
	if ( this->SPDModel->GetFeatureNum() <= 0 && this->SPDModel->GetSampleNum() <= 0)
	{
		QMessageBox mes;
		mes.setText("You haven't loaded the data file!");
		mes.exec();
	}

	std::string clusterCor = this->clusterCoherenceBox->text().toStdString();
	std::string clusterMer = this->clusterMergeBox->text().toStdString();
	if ( clusterCor.length() > 0 && clusterMer.length() > 0)
	{
		if ( atof(clusterCor.c_str()) >= 0 && atof(clusterCor.c_str()) <= 1
			&& atof(clusterMer.c_str()) >= 0 && atof(clusterMer.c_str()) <= 1)
		{
			this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()));
			this->SPDModel->ClusterMerge( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
		}
		else
		{
			QMessageBox mes;
			mes.setText("The coherence must be bewtween 0 and 1!");
			mes.exec();
		}
	}
	else
	{
		QMessageBox mes;
		mes.setText("Set coherence!");
		mes.exec();
	}
}

void SPDMainWindow::showResult()
{

}

void SPDMainWindow::generateMST()
{
	this->SPDModel->GenerateMST();
}

void SPDMainWindow::showMST()
{

}