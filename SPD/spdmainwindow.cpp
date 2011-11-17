#include "spdmainwindow.h"
#include <QGridLayout>
#include <QLCDNumber>
#include <QFileDialog>
#include <QDir>
#include <fstream>
#include <QMessageBox>
//define NDEBUG
#include <assert.h>
#include "ClusClus/clusclus.h"

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
	loadTestButton = new QPushButton(tr("Test Load"));

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

	listWidget = new QListWidget( this);
	generateMSTButton = new QPushButton(tr("MST"));
	autoProcButton = new QPushButton(tr("Auto-processing"));
	emdButton = new QPushButton(tr("EMD"));
	
	emdThresBox = new QLineEdit;
	emdPercentageBox = new QLineEdit;
	psmLable = new QLabel(tr("PSM Threshold:"));
	psmPerLable = new QLabel(tr("PSM Selected Blocks' Percentage:"));
    psmButton = new QPushButton(tr("Show PSM"));
	psmHisButton = new QPushButton(tr("PSM Histogram"));

	psdtLable = new QLabel(tr("Input hand-picked modules(seperate by comma):"));
	psdModuleSelectBox = new QLineEdit;
    psdtButton = new QPushButton(tr("View Progression"));

	//saveFeatureButton = new QPushButton(tr("Save Selected Features"));

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
	connect(loadTestButton, SIGNAL(clicked()), this, SLOT(loadTestData()));
    connect(clusterButton, SIGNAL(clicked()), this, SLOT(clusterFunction()));
	connect(generateMSTButton, SIGNAL(clicked()), this, SLOT(generateMST()));
	connect(autoProcButton, SIGNAL(clicked()), this, SLOT(autoProcess()));
	connect(emdButton, SIGNAL(clicked()), this, SLOT(emdFunction()));
	connect(psmButton, SIGNAL(clicked()), this, SLOT(showPSM()));
	connect(psmHisButton, SIGNAL(clicked()), this, SLOT(showPSMHist()));
	connect(psdtButton, SIGNAL(clicked()), this, SLOT(viewProgression()));
	//connect(saveFeatureButton, SIGNAL(clicked()), this, SLOT(saveSelectedFeatures()));
	connect(emdThresBox, SIGNAL(editingFinished()), this, SLOT(editThreshold()));
	connect(emdPercentageBox, SIGNAL(editingFinished()), this, SLOT(editPercentage()));
	
    QGridLayout *mainLayout = new QGridLayout;

    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }

    for ( int row = 1; row <= 12; row++)
    {
        mainLayout->setRowMinimumHeight(row,20);
        mainLayout->setRowStretch(row, 1);
    }

    mainLayout->addWidget(dataFileLabel, 0, 0);

    mainLayout->addWidget(dataFileName, 1, 0, 1, 2);
    mainLayout->addWidget(browseButton, 1, 2);
    mainLayout->addWidget(loadButton, 2, 2);
	mainLayout->addWidget(loadTestButton, 3, 2);

    mainLayout->addWidget(featureNumLabel, 2, 0);
    mainLayout->addWidget(featureNum, 2, 1);


    mainLayout->addWidget(sampleNumLabel, 3, 0);
    mainLayout->addWidget(sampleNum, 3, 1);

	mainLayout->addWidget(clusterCoherenceLabel, 4, 0);
    mainLayout->addWidget(clusterCoherenceBox, 4, 1);
	mainLayout->addWidget(clusterButton, 4, 2);

    mainLayout->addWidget(clusterMergeLabel, 5, 0);
    mainLayout->addWidget(clusterMergeBox, 5, 1);

	mainLayout->addWidget(listWidget, 6, 0, 3, 2);

    mainLayout->addWidget(generateMSTButton, 6, 2);
	mainLayout->addWidget(autoProcButton, 7, 2);
	mainLayout->addWidget(emdButton, 8, 2);

	mainLayout->addWidget(psmLable, 9, 0);
	mainLayout->addWidget(emdThresBox, 9, 1);
	mainLayout->addWidget(psmHisButton, 9, 2);

	mainLayout->addWidget(psmPerLable, 10, 0);
	mainLayout->addWidget(emdPercentageBox, 10, 1);
	mainLayout->addWidget(psmButton, 10, 2);

	mainLayout->addWidget(psdtLable, 11, 0);
	mainLayout->addWidget(psdModuleSelectBox, 12, 0, 1, 2);
	mainLayout->addWidget(psdtButton, 12, 2);
	//mainLayout->addWidget(saveFeatureButton, 13, 2);

    setLayout(mainLayout);

	SPDModel = SPDAnalysisModel::InitInstance();

	assert(SPDModel!=NULL);

	graph =  new GraphWindow(this);
	heatmap = new Heatmap(this);
	histo = new HistoWindow(this);
	
	connect( heatmap, SIGNAL( SelChanged()), this, SLOT( updateSelMod()));
}

SPDMainWindow::~SPDMainWindow()
{
	SPDAnalysisModel::DeInstance();
}

void SPDMainWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
{
	data = table;
	selection = sels;
	selection2 = sels2;

	if( data == NULL)
	{
		browseButton->setEnabled(TRUE);
		loadButton->setEnabled(TRUE);
		loadTestButton->setEnabled(TRUE);
	}
	else
	{
		SPDModel->ParseTraceFile( this->data);
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		this->SPDModel->NormalizeData();
		
		browseButton->setEnabled(FALSE);
		loadButton->setEnabled(FALSE);
		loadTestButton->setEnabled(FALSE);
	}
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

	if ( true == this->SPDModel->ReadCellTraceFile(file, false))
	{
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		this->SPDModel->NormalizeData();
	}
}

void SPDMainWindow::loadTestData()
{
	std::string file = this->FileName.toStdString();

	if ( true == this->SPDModel->ReadCellTraceFile(file, true))
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
			this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
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

void SPDMainWindow::generateMST()
{
	this->SPDModel->GenerateMST();
}

void SPDMainWindow::autoProcess()
{
	this->SPDModel->GenerateMST();
	this->SPDModel->RunEMDAnalysis();
	std::string str = "0";
	vtkSmartPointer<vtkTable> table = this->SPDModel->GenerateProgressionTree(str);
	if( table != NULL)
	{
		std::vector<std::string> headers;
		SPDModel->GetTableHeaders( headers);
		this->graph->setModels(data, selection, selection2);
		QString str = SPDModel->GetFileName();
		std::set<long int> featureSelectedIDs;
		SPDModel->GetSelectedFeatures(featureSelectedIDs);
		this->graph->SetTreeTable( table, headers[0], headers[1], headers[2], featureSelectedIDs, str);
		this->graph->ShowGraphWindow();
	}
}

void SPDMainWindow::emdFunction()
{
	this->SPDModel->RunEMDAnalysis();
}

void SPDMainWindow::editThreshold()
{
	std::string emdThres = this->emdThresBox->text().toStdString();
	double thres = atof(emdThres.c_str());
	double per = 0;
	if( thres >= 0 && thres <= 1)
	{
		per = this->SPDModel->GetEMDSelectedPercentage( thres);
	}
	emdPercentageBox->setText(QString::number(per));
}

void SPDMainWindow::editPercentage()
{
	std::string emdPer = this->emdPercentageBox->text().toStdString();
	double per = atof(emdPer.c_str());
	double thres = 0;
	if( per >= 0 && per <=1)
	{
		thres = this->SPDModel->GetEMDSelectedThreshold( per);
	}
	emdThresBox->setText(QString::number(thres));
}

void SPDMainWindow::showPSMHist()
{
	vtkSmartPointer<vtkTable> emdTable = vtkSmartPointer<vtkTable>::New();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName("earth mover distance");
	emdTable->AddColumn(column);
	vnl_matrix<double> emdMatrix;
	this->SPDModel->GetEMDMatrixDivByMax(emdMatrix);
	for( int i = 0; i < emdMatrix.rows(); i++)
	{
		for( int j = 0; j < emdMatrix.cols(); j++)
		{
			if( i != j)
			{
				vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
				row->InsertNextValue( vtkVariant( emdMatrix(i,j)));
				emdTable->InsertNextRow(row);
			}
		}
	}
	this->histo->setModels(emdTable);
	this->histo->show();
}

void SPDMainWindow::showPSM()
{
	std::string emdThres = this->emdThresBox->text().toStdString();
	if ( emdThres.length() > 0)
	{
		if ( atof(emdThres.c_str()) >= 0 && atof(emdThres.c_str()) <= 1)
		{
			clusclus clus1, clus2;
			this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()));
			optimalleaforder.set_size(clus1.num_samples);
			for( int i = 0; i < clus1.num_samples; i++)
			{
				optimalleaforder[i] = clus1.optimalleaforder[i];
			}
			this->heatmap->setModels();
			this->heatmap->setDataForSimilarMatrixHeatmap(clus1.features, clus1.optimalleaforder, clus2.optimalleaforder, clus1.num_samples, clus2.num_samples);	
			this->heatmap->creatDataForSimilarMatrixHeatmap();
			this->heatmap->showSimilarMatrixGraph();
		}
		else
		{
			QMessageBox mes;
			mes.setText("The threshold must be bewtween 0 and 1!");
			mes.exec();
		}
	}
	else
	{
		QMessageBox mes;
		mes.setText("Set threshold!");
		mes.exec();
	}
}

void SPDMainWindow::viewProgression()
{
	std::string selectModulesID = this->psdModuleSelectBox->text().toStdString();
	vtkSmartPointer<vtkTable> table = this->SPDModel->GenerateProgressionTree(selectModulesID);
	if( table != NULL)
	{
		
		std::vector<std::string> headers;
		SPDModel->GetTableHeaders( headers);
		this->graph->setModels(data, selection, selection2);
		QString str = SPDModel->GetFileName();
		std::set<long int> featureSelectedIDs;
		SPDModel->GetSelectedFeatures(featureSelectedIDs);
		SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", featureSelectedIDs);
		std::cout<< "Features saved in SelFeatures.txt"<<endl;
		this->graph->SetTreeTable( table, headers[0], headers[1], headers[2], featureSelectedIDs, str);
		//this->graph->SetGraphTable( table, headers[0], headers[1]);
		this->graph->ShowGraphWindow();
	}
}

void SPDMainWindow::updateSelMod()
{
	int r1 = 0;
	int r2 = 0;
	int c1 = 0;
	int c2 = 0;
	int size = optimalleaforder.size();

	this->heatmap->GetSelRowCol(r1, c1, r2, c2);
	//this->heatmap->SetSelRowCol(r1, size - 1 - r1, r2, size - 1 - r2);   // make the selection block symetric

	int num = abs(r1 - r2) + 1;
	int max = r1 > r2 ? r1 : r2;

	selMod.set_size(num);
	QString str;

	for( int i = 0; i < num; i++)
	{
		selMod[i] = optimalleaforder[size - 1 - max + i];
		if( i != num - 1)
		{
			str += QString::number(selMod[i])+",";
		}
	}
	str += QString::number(selMod[num - 1]);
	psdModuleSelectBox->setText(str);
}