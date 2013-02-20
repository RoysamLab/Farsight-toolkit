#include "clusclusMainwindow.h"
#include <QGridLayout>
#include <QLCDNumber>
#include <QFileDialog>
#include <QDir>
#include <fstream>
#include <QMessageBox>

ClusClusMainWindow::ClusClusMainWindow(QWidget *parent) :
    QWidget(parent)
{
	QGridLayout *mainLayout = new QGridLayout;

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

	rclusterNumLabel = new QLabel(tr("Samcluster number:"));
    rclusterNum = new QLabel;
    rclusterNum->setFrameStyle(frameStyle);
	cclusterNumLabel = new QLabel(tr("Feacluster number:"));
    cclusterNum = new QLabel;
    cclusterNum->setFrameStyle(frameStyle);

    rtrialsNumLable = new QLabel(tr(" Samtrials number:"));
    rtrialsNumBox = new QLineEdit;
    rgapsNumLable = new QLabel(tr("Samgap number:"));
    rgapsNumBox = new QLineEdit;

	ctrialsNumLable = new QLabel(tr(" Featrials number:"));
    ctrialsNumBox = new QLineEdit;
    cgapsNumLable = new QLabel(tr("Feagap number:"));
    cgapsNumBox = new QLineEdit;


    clusterButton = new QPushButton(tr("Runcluster"));
    biclusterButton = new QPushButton(tr("Runbicluster"));
	rcomputeGapButton = new QPushButton(tr("ComputeSamGaps"));
	ccomputeGapButton = new QPushButton(tr("ComputeFeaGaps"));
	rgenerateDendroButton = new QPushButton(tr("GenerateDendrogram1"));
	cgenerateDendroButton = new QPushButton(tr("GenerateDendrogram2"));

	singleradio = new QRadioButton("single link");
    averageradio = new QRadioButton("average link");
	completeradio = new QRadioButton("complete link");
    radiogrp = new QButtonGroup(mainLayout);

	radiogrp->addButton(singleradio);
	radiogrp->addButton(averageradio);
	radiogrp->addButton(completeradio);
	radiogrp->setId(singleradio, 1);
	radiogrp->setId(averageradio, 2);
	radiogrp->setId(completeradio, 3);

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
    connect(clusterButton, SIGNAL(clicked()), this, SLOT(runcluster()));
    connect(biclusterButton, SIGNAL(clicked()), this, SLOT(runbicluster()));
	connect(rcomputeGapButton, SIGNAL(clicked()), this, SLOT(samplecomputegap()));
	connect(ccomputeGapButton, SIGNAL(clicked()), this, SLOT(featurecomputegap()));
	connect(radiogrp, SIGNAL(buttonClicked (int)), this, SLOT(determinlinkmode(int)));
	connect(rgenerateDendroButton, SIGNAL(clicked()), this, SLOT(generatedendrogram1()));
	connect(cgenerateDendroButton, SIGNAL(clicked()), this, SLOT(generatedendrogram2()));


    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }
    for ( int row = 1; row <= 10; row++)
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
    mainLayout->addWidget(sampleNumLabel, 3, 0);
    mainLayout->addWidget(sampleNum, 3, 1);
	mainLayout->addWidget(rclusterNumLabel, 4, 0);
    mainLayout->addWidget(rclusterNum, 4, 1);
	mainLayout->addWidget(singleradio, 4, 2);

	mainLayout->addWidget(cclusterNumLabel, 5, 0);
    mainLayout->addWidget(cclusterNum, 5, 1);

    mainLayout->addWidget(rtrialsNumLable, 6, 0);
    mainLayout->addWidget(rtrialsNumBox, 6, 1);
	mainLayout->addWidget(averageradio, 5, 2);
    mainLayout->addWidget(rgapsNumLable, 7, 0);
    mainLayout->addWidget(rgapsNumBox, 7, 1);
	mainLayout->addWidget(completeradio, 6, 2);

	mainLayout->addWidget(ctrialsNumLable, 8, 0);
    mainLayout->addWidget(ctrialsNumBox, 8, 1);
    mainLayout->addWidget(cgapsNumLable, 9, 0);
    mainLayout->addWidget(cgapsNumBox, 9, 1);
    mainLayout->addWidget(rcomputeGapButton, 7, 2);

    mainLayout->addWidget(clusterButton, 10, 0);
    mainLayout->addWidget(biclusterButton, 10, 1);
    mainLayout->addWidget(ccomputeGapButton, 8, 2);

	mainLayout->addWidget(rgenerateDendroButton, 9, 2);
	mainLayout->addWidget(cgenerateDendroButton, 10, 2);

    setLayout(mainLayout);

	cc1 = NULL;
	cc2 = NULL;
	cg1 = NULL;
	cg2 = NULL;
}

ClusClusMainWindow::~ClusClusMainWindow()
{
	if(cc1)delete cc1;
	if(cc2)delete cc2;
	if(cg1)delete cg1;
	if(cg2)delete cg2;
}

void ClusClusMainWindow::browse()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileName->setText(fileName);
        this->FileName = fileName;
    }
}

void ClusClusMainWindow::load()
{
	std::string file = this->FileName.toStdString();

	this->cc1 = new clusclus();
	this->cc1->ReadFile(file.c_str());
	this->featureNum->setText( QString::number(this->cc1->num_features));
	this->sampleNum->setText( QString::number(this->cc1->num_samples));
}

void ClusClusMainWindow::runcluster()
{
	std::string rnumgaps = this->rgapsNumBox->text().toStdString();

	if(rnumgaps.length()>0)
		cc1->num_gaps = atoi(rnumgaps.c_str());

	int num_cluster;

	cc1->RunClusClus();

	cc1->MergersToProgress();

	cc1->PrepareTreeData();

	cc1->GetOptimalLeafOrderD();

	num_cluster = cc1->ComputeGapStatistics();

	this->rclusterNum->setText( QString::number(num_cluster));

	cc1->GetMembers(num_cluster);

	cc1->WriteClusteringOutputToFile("mergers.txt","features.txt","progress.txt", "members.txt",
		"gap.txt", "treedata.txt", "Optimalleaforder.txt");

}

void ClusClusMainWindow::runbicluster()
{
	cc1->Transpose();

	cc2 = new clusclus(cc1->transposefeatures,cc1->num_features, cc1->num_samples);

	std::string cnumgaps = this->cgapsNumBox->text().toStdString();

	if(cnumgaps.length()>0)
		cc2->num_gaps = atoi(cnumgaps.c_str());

	int num_cluster;

	cc2->RunClusClus();

	cc2->MergersToProgress();

	cc2->PrepareTreeData();

	cc2->GetOptimalLeafOrderD();

	num_cluster = cc2->ComputeGapStatistics();

	this->cclusterNum->setText( QString::number(num_cluster));

	cc2->GetMembers(num_cluster);

	cc2->WriteClusteringOutputToFile("mergers2.txt","features2.txt","progress2.txt", "members2.txt",
		"gap2.txt", "treedata2.txt", "Optimalleaforder2.txt");
}

void ClusClusMainWindow::samplecomputegap()
{
	std::string rnumgaps = this->rgapsNumBox->text().toStdString();
	std::string rnumtrials = this->rtrialsNumBox->text().toStdString();

	cg1 = new clusgap(cc1);

	if(rnumgaps.length()>0)
		cg1->num_gaps = atoi(rnumgaps.c_str());

	if(rnumtrials.length()>0)
		cg1->num_trials = atoi(rnumtrials.c_str());

	int num_cluster;

	num_cluster = cg1->ComputeGap();

	this->rclusterNum->setText( QString::number(num_cluster));

	this->cc1->WriteClusteringOutputToFile("Gap1.txt");
}

void ClusClusMainWindow::featurecomputegap()
{
	std::string cnumgaps = this->cgapsNumBox->text().toStdString();
	std::string cnumtrials = this->ctrialsNumBox->text().toStdString();

	cg2 = new clusgap(cc2);

	if(cnumgaps.length()>0)
		cg2->num_gaps = atoi(cnumgaps.c_str());

	if(cnumtrials.length()>0)
		cg2->num_trials = atoi(cnumtrials.c_str());

	int num_cluster;

	num_cluster = cg2->ComputeGap();

	this->cclusterNum->setText( QString::number(num_cluster));

	this->cc2->WriteClusteringOutputToFile("Gap2..txt");
}

void ClusClusMainWindow ::determinlinkmode(int id)
{
	this->cc1->linkmode = id;
	cout<<".............."<<id<<endl;
}

void ClusClusMainWindow::generatedendrogram1()
{
	Dendrogram Dendro1(0);
	Dendro1.setTreeData(this->cc1->num_samples, this->cc1->treedata, this->cc1->optimalleaforder);
	Dendro1.setModels();
}

void ClusClusMainWindow::generatedendrogram2()
{
	Dendrogram Dendro2(0);
	Dendro2.setTreeData(this->cc2->num_samples, this->cc2->treedata, this->cc2->optimalleaforder);
	Dendro2.setModels();

}