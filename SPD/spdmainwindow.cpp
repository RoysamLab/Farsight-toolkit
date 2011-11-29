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
    clusterButton = new QPushButton(tr("Feature Cluster"));
	cellClusterButton = new QPushButton(tr("Cell Cluster"));
	
	listWidget = new QListWidget( this);
	generateMSTButton = new QPushButton(tr("MST"));
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
	connect(cellClusterButton , SIGNAL(clicked()), this, SLOT(clusterCells()));
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

    for ( int row = 1; row <= 11; row++)
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
	mainLayout->addWidget(clusterButton, 5, 2);

    mainLayout->addWidget(clusterMergeLabel, 5, 0);
    mainLayout->addWidget(clusterMergeBox, 5, 1);
	mainLayout->addWidget(cellClusterButton, 4, 2);
	
	mainLayout->addWidget(listWidget, 6, 0, 2, 2);
    mainLayout->addWidget(generateMSTButton, 6, 2);
	mainLayout->addWidget(emdButton, 7, 2);

	mainLayout->addWidget(psmLable, 8, 0);
	mainLayout->addWidget(emdThresBox, 8, 1);
	mainLayout->addWidget(psmHisButton, 8, 2);

	mainLayout->addWidget(psmPerLable, 9, 0);
	mainLayout->addWidget(emdPercentageBox, 9, 1);
	mainLayout->addWidget(psmButton, 9, 2);

	mainLayout->addWidget(psdtLable, 10, 0);
	mainLayout->addWidget(psdModuleSelectBox, 11, 0, 1, 2);
	mainLayout->addWidget(psdtButton, 11, 2);
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
		//this->SPDModel->NormalizeData();
		
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
		//this->SPDModel->NormalizeData();
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
			try
			{
				this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
				this->SPDModel->ClusterMerge( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
			}
			catch(...)
			{
				std::cout<< "Clustering exception, please try again!"<<endl;
			}
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
	try
	{
		this->SPDModel->GenerateMST();
	}
	catch(...)
	{
		std::cout<< "MST construction failure, please try again!"<<endl;
	}
}

void SPDMainWindow::clusterCells()
{
	std::string clusterCor = this->clusterCoherenceBox->text().toStdString();
	if ( atof(clusterCor.c_str()) >= 0 && atof(clusterCor.c_str()) <= 1)
	{
		this->SPDModel->ClusterCells(atof(clusterCor.c_str()));
	}
}

void SPDMainWindow::emdFunction()
{
	try
	{
		this->SPDModel->RunEMDAnalysis();
	}
	catch(...)
	{
		std::cout<< "EMD construction failure, please try again!"<<endl;
	}
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
		std::vector<int> index;

		SPDModel->GetTableHeaders( headers);
		SPDModel->GetClusterMapping(index);
		if( index.size() > 0)
		{
			this->graph->setModels(data, selection, selection2, &index);
		}
		else
		{
			this->graph->setModels(data, selection, selection2);
		}

		QString str = SPDModel->GetFileName();
		std::set<long int> featureSelectedIDs;
		SPDModel->GetSelectedFeatures(featureSelectedIDs);
		SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", featureSelectedIDs);
		std::cout<< "Features saved in SelFeatures.txt"<<endl;
		this->graph->SetTreeTable( table, headers[0], headers[1], headers[2], featureSelectedIDs, str);
		//this->graph->SetGraphTable( table, headers[0], headers[1]);
		try
		{
			this->graph->ShowGraphWindow();
		}
		catch(...)
		{
			std::cout<< "Graph window error!"<<endl;
		}
	}
}

void SPDMainWindow::updateSelMod()   // possible bugs
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

	if( num < size)
	{
		selMod.set_size(num);
		QString str;
		
		for( int i = 0; i < num; i++)
		{
			int index = size - 1 - max + i;
			if( index >= 0 && index < size)
			{
				selMod[i] = optimalleaforder[index];
				if( i != num - 1)
				{
					str += QString::number(selMod[i])+",";
				}
			}
		}
		str += QString::number(selMod[num - 1]);
		psdModuleSelectBox->setText(str);
	}
}

void SPDMainWindow::GetProgressionTreeOrder(std::vector<long int> &order)
{
	this->graph->GetProgressionTreeOrder(order);
}