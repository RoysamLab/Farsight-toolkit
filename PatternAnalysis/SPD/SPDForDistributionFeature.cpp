#include "SPDForDistributionFeature.h"
#include <QGridLayout>
#include <QLCDNumber>
#include <QFileDialog>
#include <QDir>
#include <fstream>
#include <QMessageBox>
//define NDEBUG
#include <assert.h>
#include "ClusClus/clusclus.h"

#define MSTSPD 0

using std::ifstream;
using std::endl;

SPDForDistributionWindow::SPDForDistributionWindow(QWidget *parent) :
    QWidget(parent)
{
	SPDModel = NULL;
	selection = NULL;
	selection2 = NULL;
	graph = NULL;
	simHeatmap = NULL;
	histo = NULL;
	progressionHeatmap = NULL;
	HeatmapWin = NULL;
	plot = NULL;

    dataFileLabel = new QLabel(tr("Choose file:"));

    int frameStyle = QFrame::Sunken | QFrame::Panel;
    dataFileName = new QLabel;
    dataFileName->setFrameStyle(frameStyle);

    browseButton = new QPushButton(tr("Browse"));
    loadButton = new QPushButton(tr("Load"));
	loadTestButton = new QPushButton(tr("Load Contrast Data"));

    featureNumLabel = new QLabel(tr("Feature size:"));
    featureNum = new QLabel;
    featureNum->setFrameStyle(frameStyle);
    sampleNumLabel = new QLabel(tr("Sample size:"));
    sampleNum = new QLabel;
    sampleNum->setFrameStyle(frameStyle);
	sampleNumLabel = new QLabel(tr("Bin size:"));
    sampleNum = new QLabel;
    sampleNum->setFrameStyle(frameStyle);

	emdThresBox = new QDoubleSpinBox;
	emdThresBox->setRange(0,1);
	emdThresBox->setSingleStep(0.1);

	emdPercentageBox = new QLineEdit;
	psmLable = new QLabel(tr("PSM Threshold(0.0 ~ 1.0):"));
	psmPerLable = new QLabel(tr("PSM Selected Blocks' Percentage:"));
    psmButton = new QPushButton(tr("Show PSM"));

	psdtLable = new QLabel(tr("Input hand-picked modules(seperate by comma):"));
	psdModuleSelectBox = new QLineEdit;
	
    psdtButton = new QPushButton(tr("View Progression"));
	heatmapLabel = new QLabel(tr("View Progression Heatmap:"));
	
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);

	//saveFeatureButton = new QPushButton(tr("Save Selected Features"));

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
	
	connect(psmButton, SIGNAL(clicked()), this, SLOT(showPSM()));
	connect(psdtButton, SIGNAL(clicked()), this, SLOT(viewProgression()));
	//connect(saveFeatureButton, SIGNAL(clicked()), this, SLOT(saveSelectedFeatures()));
	connect(emdThresBox, SIGNAL(editingFinished()), this, SLOT(editThreshold()));
	connect(emdPercentageBox, SIGNAL(editingFinished()), this, SLOT(editPercentage()));
	connect(heatmapButton, SIGNAL(clicked()), this, SLOT(showProgressionHeatmap()));
	
    QGridLayout *mainLayout = new QGridLayout;

    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }

     for ( int row = 1; row <= 9; row++)
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

	mainLayout->addWidget(binNumLabel, 4, 0);
    mainLayout->addWidget(binNum, 4, 1);

	mainLayout->addWidget(psmLable, 5, 0);
	mainLayout->addWidget(emdThresBox, 5, 1);

	mainLayout->addWidget(psmPerLable, 6, 0);
	mainLayout->addWidget(emdPercentageBox, 6, 1);
	mainLayout->addWidget(psmButton, 6, 2);

	mainLayout->addWidget(psdtLable, 7, 0);
	mainLayout->addWidget(psdModuleSelectBox, 8, 0, 1, 2);
	mainLayout->addWidget(psdtButton, 8, 2);

	mainLayout->addWidget(heatmapLabel, 9, 0);
	mainLayout->addWidget(heatmapButton, 9, 2);

    setLayout(mainLayout);

	SPDModel = SPDAnalysisModel::InitInstance();
}

SPDForDistributionWindow::~SPDForDistributionWindow()
{
	SPDAnalysisModel::DeInstance();
}

void SPDForDistributionWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
{
	data = table;

	if( sels == NULL)
	{
		selection = new ObjectSelection();
	}
	else
	{
		selection = sels;
	}

	if( sels2 == NULL)
	{
		selection2 = new ObjectSelection();
	}
	else
	{
		selection2 = sels2;
	}

	if( data == NULL)
	{
		browseButton->setEnabled(TRUE);
		loadButton->setEnabled(TRUE);
	}
	else
	{
		SPDModel->ParseTraceFile( data);
		
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		//this->SPDModel->NormalizeData();
		
		browseButton->setEnabled(FALSE);
		loadButton->setEnabled(FALSE);
	}

	assert(SPDModel!=NULL);

	if(this->simHeatmap)
	{
		delete this->simHeatmap;
	}
	this->simHeatmap = new ProgressionHeatmap( this);
	connect( simHeatmap, SIGNAL( SelChanged()), this, SLOT( updateSelMod()));

	if(this->graph)
	{
		delete this->graph;
	}
	this->graph = new GraphWindow( this);


}

void SPDForDistributionWindow::browse()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileName->setText(fileName);
        this->FileName = fileName;
    }
}

void SPDForDistributionWindow::load()
{
	std::string file = this->FileName.toStdString();

	if ( true == this->SPDModel->ReadCellTraceFile(file, false))
	{
		data = this->SPDModel->GetDataTable();
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		//this->SPDModel->NormalizeData();
	}
}

void SPDForDistributionWindow::loadContrastData()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        std::string file = fileName.toStdString();
		if ( true == this->SPDModel->ReadCellTraceFile(file, true))  // add the contrast data to the previous data
		{
			this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
			this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum() + this->SPDModel->GetContrastDataSampleNum()));
			//this->SPDModel->NormalizeData();
		}
    }
}


void SPDForDistributionWindow::editThreshold()
{
	QString emdThres = this->emdThresBox->text();
	double thres = this->emdThresBox->valueFromText(emdThres);
	double per = 0;
	if( thres >= 0 && thres <= 1)
	{
#if MSTSPD
		per = this->SPDModel->GetEMDSelectedPercentage( thres);
#else
		per = this->SPDModel->GetCorMatSelectedPercentage( thres);
#endif
	}
	emdPercentageBox->setText(QString::number(per));
}

void SPDForDistributionWindow::editPercentage()
{
	//std::string emdPer = this->emdPercentageBox->text().toStdString();
	//double per = atof(emdPer.c_str());
	//double thres = 0;
	//if( per >= 0 && per <=1)
	//{
	//	thres = this->SPDModel->GetEMDSelectedThreshold( per);
	//}
	//emdThresBox->setText(QString::number(thres));
}

void SPDForDistributionWindow::showPSM()
{
	std::string emdThres = this->emdThresBox->text().toStdString();

	clusclus *clus1 = new clusclus();
	clusclus *clus2 = new clusclus();
	std::vector< unsigned int> moduleIDs;
	if( SPDModel->GetProgressionType())
	{
#if MSTSPD
		this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
#else

		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
#endif

		QString str;
		int i = 0;
		if( moduleIDs.size() > 0)
		{
			for( i = 0; i < moduleIDs.size() - 1; i++)
			{
				str += QString::number(moduleIDs[i])+",";
			}
			str += QString::number(moduleIDs[i]);
			psdModuleSelectBox->setText(str);
			psdtButton->setEnabled(TRUE);
		}
	}
	else
	{
#if MSTSPD
		this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
#else

		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
#endif
		optimalleaforder.set_size(clus1->num_samples);
		for( int i = 0; i < clus1->num_samples; i++)
		{
			optimalleaforder[i] = clus1->optimalleaforder[i];
		}
		this->simHeatmap->setModels();
		this->simHeatmap->setDataForSimilarMatrixHeatmap(clus1->features, clus1->optimalleaforder, clus2->optimalleaforder, clus1->num_samples, clus2->num_samples);	
		this->simHeatmap->creatDataForSimilarMatrixHeatmap();
		this->simHeatmap->showSimilarMatrixGraph();
	}
	delete clus1;
	delete clus2;
}

void SPDForDistributionWindow::viewProgression()
{
	/* needs to be changed:
	   get orders of the features
	   setmodels heatmap: table after dimension reduced in sample space
	   selection: threshold selection, node selection
	   show default heatmap and tree here
	*/
	
	/** heatmap set models */

	if( this->HeatmapWin)
	{
		delete this->HeatmapWin;
	}
	this->HeatmapWin = new Heatmap(this);


	std::string selectModulesID = this->psdModuleSelectBox->text().toStdString();
	std::vector< unsigned int> selModuleID;

	std::vector< int> clusterSize;
	selFeatureID.clear();
	selOrder.clear();
	unselOrder.clear();

	split( selectModulesID, ',', selModuleID);
	SPDModel->GetFeatureIdbyModId(selModuleID, selFeatureID);
	if( SPDModel->GetProgressionType())
	{
		SPDModel->SaveSelectedFeatureNames("DistanceProgressionSelFeatures.txt", selFeatureID);
	}
	else
	{
		SPDModel->SaveSelectedFeatureNames("SPDtestWindow_ProgressionSelFeatures.txt", selFeatureID);
	}

	GetFeatureOrder( selFeatureID, selOrder, unselOrder);

	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();
	connect(selection, SIGNAL( thresChanged()), this, SLOT( regenerateProgressionTree()));
	connect(selection, SIGNAL( ItemDeleted()), this, SLOT( ReRunSPDAnlysis()));
	connect(HeatmapWin, SIGNAL(columnToColorChanged(int)), this, SLOT( ReColorProgressionTree(int)));

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);
	this->HeatmapWin->setModelsforSPD( tableAfterCellCluster, selection, selOrder, unselOrder, &indexMap);
	this->HeatmapWin->showGraphforSPD( selOrder.size(), unselOrder.size());


	//vtkSmartPointer<vtkTable> table = this->SPDModel->GenerateProgressionTree(selectModulesID);
	//std::vector<std::string> headers;
	//SPDModel->GetTableHeaders( headers);
	//QString str = SPDModel->GetFileName();
	//std::set<long int> featureSelectedIDs;
	//SPDModel->GetSelectedFeatures(featureSelectedIDs);
	//SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", featureSelectedIDs);


	//this->graph->setModels(data, selection);
	//this->graph->SetTreeTable( table, headers[0], headers[1], headers[2]);
	//try
	//{
	//	this->graph->ShowGraphWindow();
	//	//this->HeatmapWin->showGraph();
	//}
	//catch(...)
	//{
	//	std::cout<< "Graph window error!"<<endl;
	//}
}

void SPDForDistributionWindow::split(std::string& s, char delim, std::vector< unsigned int>& indexVec)
{
	size_t last = 0;
	size_t index = s.find_first_of(delim,last);
	std::vector< std::string > stringVec;
	while( index!=std::string::npos)
	{
		stringVec.push_back(s.substr(last,index-last));
		last = index+1;
		index=s.find_first_of(delim,last);
	}
	if( index-last>0)
	{
		stringVec.push_back(s.substr(last,index-last));
	}

	for( int i = 0; i < stringVec.size(); i++)
	{
		unsigned int index = atoi( stringVec[i].c_str());
		indexVec.push_back( index);
	}
} 

void SPDForDistributionWindow::updateSelMod()   // possible bugs
{
	int r1 = 0;
	int r2 = 0;
	int c1 = 0;
	int c2 = 0;
	int size = optimalleaforder.size();

	this->simHeatmap->GetSelRowCol(r1, c1, r2, c2);
	//this->simHeatmap->SetSelRowCol(r1, size - 1 - r1, r2, size - 1 - r2);   // make the selection block symetric

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
		psdtButton->setEnabled(TRUE);
		heatmapButton->setEnabled(FALSE);
	}
}

void SPDForDistributionWindow::GetProgressionTreeOrder(std::vector<long int> &order)
{
	this->graph->GetProgressionTreeOrder(order);
}

void SPDForDistributionWindow::showProgressionHeatmap()
{
	//ofstream ofs("SPDHeatmapOptimalOrder.txt");
	if( this->progressionHeatmap)
	{
		delete this->progressionHeatmap;
	}
	this->progressionHeatmap = new Heatmap(this);
	
	std::vector<long int> TreeOrder;
	this->graph->GetProgressionTreeOrder(TreeOrder);   // order of the cluster 
	if( TreeOrder.size() <=0)
	{          
		std::cout<< "progression tree hasn't been built yet"<<endl;
		return;
	}

	std::vector< std::vector< long int> > sampleIndex;
	selection->GetSampleIndex( sampleIndex);
	std::vector< std::vector< long int> > clusIndex;
	selection->GetClusterIndex( clusIndex);
	std::vector< int> clusterOrder;
	SPDModel->GetClusterOrder(clusIndex, TreeOrder, clusterOrder);

	// module feature and percentage plot
	std::vector< double> percentageOfSamples;
	std::vector< double> percentageOfNearDeviceSamples;
	std::string distanceThres = this->distanceThres->text().toStdString();
	SPDModel->GetPercentage(sampleIndex, percentageOfSamples);
	SPDModel->GetCloseToDevicePercentage(sampleIndex, percentageOfNearDeviceSamples, atof(distanceThres.c_str()));

	vtkSmartPointer<vtkTable> tableForAverModulePlot = SPDModel->GetAverModuleTable(clusIndex, TreeOrder, percentageOfSamples, percentageOfNearDeviceSamples, selOrder, unselOrder);
	if( plot)
	{
		delete plot;
	}
	plot = new PlotWindow(this);
	plot->setModels(tableForAverModulePlot, selection);
	plot->show();

	// progression heatmap
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);
	this->progressionHeatmap->setModelsforSPD( tableAfterCellCluster, selection, clusterOrder, selOrder, unselOrder, &indexMap);
	this->progressionHeatmap->showGraphforSPD( selOrder.size(), unselOrder.size(), true);
}

void SPDForDistributionWindow::closeEvent(QCloseEvent *event)
{
	closeSubWindows();
	event->accept();
}

void SPDForDistributionWindow::closeSubWindows()
{
	if(graph)
	{
		graph->close();
	}
	if(simHeatmap)
	{
		simHeatmap->close();
	}
	if(progressionHeatmap)
	{
		progressionHeatmap->close();
	}
	if(HeatmapWin)
	{
		HeatmapWin->close();
	}
	if(plot)
	{
		plot->close();
	}
}