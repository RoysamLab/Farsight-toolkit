#include "spdtestwindow.h"
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

#ifndef NEW_DELETE
#define NEW_DELETE
void operator delete(void *memP)
{
   if(memP)
   {
	   free(memP);
	   memP = NULL;
   }
};
#endif

using std::ifstream;
using std::endl;

SPDtestWindow::SPDtestWindow(QWidget *parent) :
    QWidget(parent)
{
	SPDModel = NULL;
	selection = NULL;
	selection2 = NULL;
	graph = NULL;
	simHeatmap = NULL;
	histo = NULL;
	originalHeatmap = NULL;
	progressionHeatmap = NULL;
	HeatmapWin = NULL;
	plot = NULL;
	connectedNum = 0;
	bconnected = true;

    dataFileLabel = new QLabel(tr("Choose file:"), this);

    int frameStyle = QFrame::Sunken | QFrame::Panel;
    dataFileName = new QLabel(this);
    dataFileName->setFrameStyle(frameStyle);

    browseButton = new QPushButton(tr("Browse"), this);
    loadButton = new QPushButton(tr("Load"), this);
	//loadTestButton = new QPushButton(tr("Raw Data Heatmap"), this);

    featureNumLabel = new QLabel(tr("Feature size:"), this);
    featureNum = new QLabel(this);
    featureNum->setFrameStyle(frameStyle);
    sampleNumLabel = new QLabel(tr("Sample size:"), this);
    sampleNum = new QLabel(this);
    sampleNum->setFrameStyle(frameStyle);

    clusterCoherenceLabel = new QLabel(tr("Feature Coherence(0.0 ~ 1.0):"), this);
    clusterCoherenceBox = new QDoubleSpinBox(this);
	clusterCoherenceBox->setValue(0.95);
	clusterCoherenceBox->setRange(0,1); 
	clusterCoherenceBox->setSingleStep(0.1);

    kNearestNeighborLabel = new QLabel(tr("Nearest Neighbor Number:"), this);
    kNearestNeighborBox = new QSpinBox(this);
	kNearestNeighborBox->setValue(5);
	kNearestNeighborBox->setMinimum (0);
	kNearestNeighborBox->setSingleStep(1);

    //clusterButton = new QPushButton(tr("Feature Cluster"), this);
	
	emdLabel = new QLabel(tr("KNNG based EMD module matching:"), this);
	//progressionOverDistance = new QLabel(tr("Progression over distance to device:"), this);
	//bcheckBox = new QCheckBox(this);
	emdButton = new QPushButton(tr("Match"), this);

	emdThresBox = new QDoubleSpinBox(this);
	emdThresBox->setRange(0,1);
	emdThresBox->setValue(0.3);

	emdPercentageBox = new QLineEdit(this);
	emdPercentageBox->setReadOnly( true);
	emdPercentageBox->setEnabled( false);
	psmLable = new QLabel(tr("PSM Threshold(0.0 ~ 1.0):"), this);
	psmPerLable = new QLabel(tr("PSM Selected Blocks' Percentage:"), this);
    psmButton = new QPushButton(tr("Show PSM"), this);

	psdtLable = new QLabel(tr("Input hand-picked modules(seperate by comma):"), this);
	psdModuleSelectBox = new QLineEdit(this);
	maxVetexIdLabel = new QLabel(tr("Id to seperate:"), this);
	maxVetexIdEdit = new QSpinBox(this);
	maxVetexIdEdit->setRange(0,1000000);
	maxVetexIdEdit->setSingleStep(100);
	maxVetexIdEdit->setValue(4500);

	//searchSubsetsButton = new QPushButton(tr("Search subsets"), this);

	connectedGraphLabel = new QLabel(tr("Connected Graphs:"), this);
	connectedGraphEdit = new QLineEdit(this);
	connectedGraphEdit->setText(QString::number(0));
	connectedGraphEdit->setReadOnly( true);
	connectedGraphEdit->setEnabled( false);
	updateConnectedNumButton = new QPushButton(tr("Update"), this);
	
    psdtButton = new QPushButton(tr("View Progression"), this);
	heatmapLabel = new QLabel(tr("View Progression Heatmap:"), this);
	heatmapButton = new QPushButton(tr("Heatmap"), this);

	distanceLabel = new QLabel(tr("Distance Threshold:"), this);
	distanceThres = new QDoubleSpinBox(this);
	distanceThres->setRange(0,2000);
	distanceThres->setSingleStep(0.1);
	distanceThres->setValue(700.0);
	
	emdButton->setEnabled(TRUE);
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	updateConnectedNumButton->setEnabled(FALSE);
	//searchSubsetsButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);
	//bcheckBox->setChecked(FALSE);
	//saveFeatureButton = new QPushButton(tr("Save Selected Features"));

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
	//connect(loadTestButton, SIGNAL(clicked()), this, SLOT(showOriginalHeatmap()));
    //connect(clusterButton, SIGNAL(clicked()), this, SLOT(clusterFunction()));
	//connect( bcheckBox, SIGNAL(clicked()), this, SLOT(updateProgressionType()));
	connect( kNearestNeighborBox, SIGNAL(editingFinished()), this, SLOT(editNearestNeighbor()));
	connect( updateConnectedNumButton, SIGNAL(clicked()), this, SLOT(UpdateConnectedNum()));
	//connect( searchSubsetsButton, SIGNAL(clicked()), this, SLOT(searchSubsetsOfFeatures()));
	connect(heatmapButton, SIGNAL(clicked()), this, SLOT(showProgressionHeatmap()));

	connect(emdButton, SIGNAL(clicked()), this, SLOT(emdFunction()));
	connect(psmButton, SIGNAL(clicked()), this, SLOT(showPSM()));
	connect(psdtButton, SIGNAL(clicked()), this, SLOT(viewProgression()));
	connect(emdThresBox, SIGNAL(editingFinished()), this, SLOT(editThreshold()));
	//connect(emdPercentageBox, SIGNAL(editingFinished()), this, SLOT(editPercentage()));
	
    QGridLayout *mainLayout = new QGridLayout(this);

    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }

     for ( int row = 1; row <= 14; row++)
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
	//mainLayout->addWidget(loadTestButton, 3, 2);

	mainLayout->addWidget(clusterCoherenceLabel, 4, 0);
	mainLayout->addWidget(clusterCoherenceBox, 4, 1);
	//mainLayout->addWidget(clusterButton, 4, 2);

    mainLayout->addWidget(kNearestNeighborLabel, 5, 0);
    mainLayout->addWidget(kNearestNeighborBox, 5, 1);
	
	//mainLayout->addWidget(progressionOverDistance, 6, 0);
	//mainLayout->addWidget(bcheckBox, 6, 1);

	mainLayout->addWidget(emdLabel, 6, 0);
	mainLayout->addWidget(emdButton, 6, 2);

	mainLayout->addWidget(psmLable, 7, 0);
	mainLayout->addWidget(emdThresBox, 7, 1);

	mainLayout->addWidget(psmPerLable, 8, 0);
	mainLayout->addWidget(emdPercentageBox, 8, 1);
	mainLayout->addWidget(psmButton, 8, 2);

	mainLayout->addWidget(psdtLable, 10, 0);
	mainLayout->addWidget(distanceLabel, 9, 0);
	mainLayout->addWidget(distanceThres, 9, 1);
	mainLayout->addWidget(psdModuleSelectBox, 11, 0, 1, 2);
	//mainLayout->addWidget(searchSubsetsButton, 11, 2);

	mainLayout->addWidget(connectedGraphLabel, 12, 0);
	mainLayout->addWidget(connectedGraphEdit, 12, 1);
	mainLayout->addWidget(updateConnectedNumButton, 12, 2);
	mainLayout->addWidget(maxVetexIdLabel, 13, 0);
	mainLayout->addWidget(maxVetexIdEdit, 13, 1);
	mainLayout->addWidget(psdtButton, 13, 2);

	mainLayout->addWidget(heatmapLabel, 14, 0);
	mainLayout->addWidget(heatmapButton, 14, 2);


    setLayout(mainLayout);

	SPDModel = new SPDAnalysisModel();
}

SPDtestWindow::~SPDtestWindow()
{
	delete(SPDModel);
	//delete(selection);
	delete(selection2);
}

void SPDtestWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
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
		//loadTestButton->setEnabled(TRUE);
	}
	else
	{
		SPDModel->ParseTraceFile( data);
		
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		
		browseButton->setEnabled(FALSE);
		loadButton->setEnabled(FALSE);
		//loadTestButton->setEnabled(TRUE);
	}

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

void SPDtestWindow::browse()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileName->setText(fileName);
        this->FileName = fileName;
    }
}

void SPDtestWindow::load()
{
	std::string file = this->FileName.toStdString();

	if ( true == this->SPDModel->ReadRawData(file))
	{
		data = this->SPDModel->GetDataTable();
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
	}
}

void SPDtestWindow::loadContrastData()
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

void SPDtestWindow::showOriginalHeatmap()
{
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();
	if( this->originalHeatmap)
	{
		delete this->originalHeatmap;
		this->originalHeatmap = NULL;
	}

	this->originalHeatmap = new Heatmap(this);
	std::vector< int> sampleOrder;
	for( int i = 0; i < tableAfterCellCluster->GetNumberOfRows(); i++)
	{
		sampleOrder.push_back(i);
	}
	std::vector< int> selOrder;
	for( int i = 0; i < tableAfterCellCluster->GetNumberOfColumns(); i++)
	{
		selOrder.push_back(i);
	}
	std::vector< int> unselOrder;

	this->originalHeatmap->setModelsforSPD( tableAfterCellCluster, selection, sampleOrder, selOrder, unselOrder);
	this->originalHeatmap->showGraphforSPD( selOrder.size(), unselOrder.size(), true);
}

void SPDtestWindow::showHeatmapAfterFeatureClustering()
{
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();
	if( this->originalHeatmap)
	{
		delete this->originalHeatmap;
		this->originalHeatmap = NULL;
	}

	this->originalHeatmap = new Heatmap(this);
	std::vector< int> sampleOrder;
	for( int i = 0; i < tableAfterCellCluster->GetNumberOfRows(); i++)
	{
		sampleOrder.push_back(i);
	}
	std::vector< int> selOrder;
	std::vector< unsigned int> inputSelOrder;
	for( unsigned int i = 0; i < tableAfterCellCluster->GetNumberOfColumns(); i++)
	{
		inputSelOrder.push_back(i);
	}
	GetFeatureOrder(inputSelOrder, selOrder, unselOrder);

	this->originalHeatmap->setModelsforSPD( tableAfterCellCluster, selection, sampleOrder, selOrder, unselOrder);
	this->originalHeatmap->showGraphforSPD( selOrder.size(), unselOrder.size(), true);
}

void SPDtestWindow::clusterFunction()
{
	if ( this->SPDModel->GetFeatureNum() <= 0 && this->SPDModel->GetSampleNum() <= 0)
	{
		QMessageBox mes;
		mes.setText("You haven't loaded the data file!");
		mes.exec();
	}

	std::string clusterCor = this->clusterCoherenceBox->text().toStdString();
	
	try
	{
		this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()), atof(clusterCor.c_str()));
		//showHeatmapAfterFeatureClustering();
	}
	catch(...)
	{
		std::cout<< "Clustering exception, please try again!"<<endl;
	}
}

//void SPDtestWindow::updateProgressionType()
//{  
//	if( bcheckBox->isChecked())
//	{
//		bool rtn = this->SPDModel->SetProgressionType( true);   // progression over distance to device
//		if( rtn == false)
//		{
//			QMessageBox mes;
//			mes.setText("Distance to device is not available!");
//			mes.exec();
//		}
//	}
//	else
//	{
//		this->SPDModel->SetProgressionType( false);  // overall progression
//	}
//}

void SPDtestWindow::emdFunction()
{
	clusterFunction();
	try
	{
#if MSTSPD
		this->SPDModel->GenerateMST();
		this->SPDModel->GenerateDistanceMST();
		this->SPDModel->RunEMDAnalysis();
#else
		//this->SPDModel->ModuleCoherenceMatchAnalysis();
		std::string kNeighbor = this->kNearestNeighborBox->text().toStdString();
		this->SPDModel->ModuleCorrelationMatrixMatch(atoi(kNeighbor.c_str()));
#endif
		psmButton->setEnabled(TRUE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "EMD construction failure, please try again!"<<endl;
	}
}

void SPDtestWindow::editThreshold()
{
	QString emdThres = this->emdThresBox->text();
	double thres = emdThres.toDouble();
	double per = 0;
	if( thres >= 0 && thres <= 1)
	{
//#if MSTSPD
		per = this->SPDModel->GetEMDSelectedPercentage( thres);
//#else
//		per = this->SPDModel->GetCorMatSelectedPercentage( thres);
//#endif
	}
	emdPercentageBox->setText(QString::number(per));
}

//void SPDtestWindow::editPercentage()
//{
//	//std::string emdPer = this->emdPercentageBox->text().toStdString();
//	//double per = atof(emdPer.c_str());
//	//double thres = 0;
//	//if( per >= 0 && per <=1)
//	//{
//	//	thres = this->SPDModel->GetEMDSelectedThreshold( per);
//	//}
//	//emdThresBox->setText(QString::number(thres));
//}

void SPDtestWindow::editNearestNeighbor()
{
	QString str = this->kNearestNeighborBox->text();
	int numNeighbor = str.toInt();
	int kNum = this->SPDModel->GetKNeighborNum();
	if( kNum <= 0 || kNum != numNeighbor)
	{
		psmButton->setEnabled(FALSE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
}

void SPDtestWindow::showPSM()
{
	editThreshold();

	std::string emdThres = this->emdThresBox->text().toStdString();

	clusclus *clus1 = new clusclus();
	clusclus *clus2 = new clusclus();
	std::vector< unsigned int> moduleIDs;
	if( SPDModel->GetProgressionType())
	{
//#if MSTSPD
		this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
//#else
//		
//		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
//#endif

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
//#if MSTSPD
		this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
//#else
//
//		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);
//#endif
		optimalleaforder.set_size(clus1->num_samples);


		for( int i = 0; i < clus1->num_samples; i++)
		{
			optimalleaforder[i] = clus1->optimalleaforder[i];
		}

		//for( int i = 0; i < clus1->num_samples; i++)
		//{
		//	clus1->optimalleaforder[i] = i;
		//}
		//for( int i = 0; i < clus2->num_samples; i++)
		//{
		//	clus2->optimalleaforder[i] = i;
		//}

		this->simHeatmap->setModels();
		this->simHeatmap->setDataForSimilarMatrixHeatmap(clus1->features, clus1->optimalleaforder, clus2->optimalleaforder, clus1->num_samples, clus2->num_samples);	
		this->simHeatmap->creatDataForSimilarMatrixHeatmap();
		this->simHeatmap->showSimilarMatrixGraph();
	}
	delete clus1;
	delete clus2;
}

void SPDtestWindow::viewProgression()
{	
	UpdateConnectedNum();  // update connected component
	if(connectedNum <= 1e-9)
	{
		return;
	}

	if( this->originalHeatmap)
	{
		delete this->originalHeatmap;
		this->originalHeatmap = NULL;
	}

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

	GetFeatureOrder( selFeatureID, selOrder, unselOrder);
	//SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", selOrder);

	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();

	connect(selection, SIGNAL( thresChanged()), this, SLOT( regenerateProgressionTree()));
	connect(selection, SIGNAL( ItemDeleted()), this, SLOT( ReRunSPDAnlysis()));
	connect(HeatmapWin, SIGNAL(columnToColorChanged(int)), this, SLOT( ReColorProgressionTree(int)));

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);

	//SPDModel->BuildMSTForConnectedComponent(selFeatureID, connectedComponent, connectedNum);

	//std::vector< Tree> sampleTree = SPDModel->mstTreeList;
	//double **streedata = new double*[sampleTree.size()];
	//for(int i = 0; i < sampleTree.size(); i++)
	//{
	//	streedata[i] = new double[4];
	//	streedata[i][0] = sampleTree[i].first;
	//	streedata[i][1] = sampleTree[i].second;
	//	streedata[i][2] = sampleTree[i].cor;
	//	streedata[i][3] = sampleTree[i].parent;
	//}

	vnl_matrix<double> subTreeDistance;
	SPDModel->GetComponentMinDistance(selFeatureID, connectedComponent, connectedNum, subTreeDistance);
	this->HeatmapWin->setModelsforSPD( tableAfterCellCluster, selection, selOrder, unselOrder, &indexMap,\
										&connectedComponent, connectedNum, &subTreeDistance);
	this->HeatmapWin->showGraphforSPD( selOrder.size(), unselOrder.size());

	//for(int i = 0; i < sampleTree.size(); i++)
	//{
	//	delete streedata[i];
	//}
	//delete streedata;
}


void SPDtestWindow::searchSubsetsOfFeatures()
{
	//std::string selectModulesID = this->psdModuleSelectBox->text().toStdString();
	std::vector< unsigned int> selModuleID;
	//split( selectModulesID, ',', selModuleID);
	bool bstate = SPDModel->SearchSubsetsOfFeatures(selModuleID);
	if( false == bstate)
	{
		QMessageBox mes;
		mes.setText("Distance not available, target unclear.");
		mes.exec();
	}
}

void SPDtestWindow::split(std::string& s, char delim, std::vector< unsigned int>& indexVec)
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

void SPDtestWindow::GetFeatureOrder(std::vector< unsigned int> &selID, std::vector<int> &selIdOrder, std::vector<int> &unselIdOrder)
{
	SPDModel->HierachicalClustering();
	std::vector< Tree> FeatureTreeData = SPDModel->PublicTreeData;
	double **ftreedata = new double*[FeatureTreeData.size()];

	for(int i = 0; i < FeatureTreeData.size(); i++)
	{
		ftreedata[i] = new double[4];
		ftreedata[i][0] = FeatureTreeData[i].first;
		ftreedata[i][1] = FeatureTreeData[i].second;
		ftreedata[i][2] = (1 - FeatureTreeData[i].cor + 0.01) * 100;
		ftreedata[i][3] = FeatureTreeData[i].parent;
	}

	clusclus *cc2 = new clusclus();
	cc2->Initialize(ftreedata, FeatureTreeData.size() + 1);
	cc2->GetOptimalLeafOrderD();

	for( int i = 0; i < cc2->num_samples; i++)
	{
		if( IsExist(selID, (unsigned int)cc2->optimalleaforder[i]))
		{
			selIdOrder.push_back( cc2->optimalleaforder[i]);
		}
		else
		{
			unselIdOrder.push_back( cc2->optimalleaforder[i]);
		}
	}

	//ofstream ofs("FeatureOrder.txt");
	//ofs<< "feature optimal order:"<<endl;
	//for( int i = 0; i < cc2->num_samples; i++)
	//{
	//	ofs<< cc2->optimalleaforder[i]<<"\t";
	//}
	//ofs<<endl;
	//ofs<< "Selected features optimal order:"<<endl;
	//for( int i = 0; i < selIdOrder.size(); i++)
	//{
	//	ofs<< selIdOrder[i]<<"\t";
	//}
	//ofs<<endl;
	//ofs<< "UnSelected features optimal order:"<<endl;
	//for( int i = 0; i < unselIdOrder.size(); i++)
	//{
	//	ofs<< unselIdOrder[i]<<"\t";
	//}
	//ofs<<endl;
	//ofs.close();

	for( int i = 0; i < FeatureTreeData.size(); i++)
	{
		delete ftreedata[i];
	}
	delete ftreedata;
	delete cc2;
}

bool SPDtestWindow::IsExist(std::vector< unsigned int> vec, unsigned int value)
{
	for( int i = 0; i < vec.size(); i++)
	{
		if( value == vec[i])
		{
			return true;
		}
	}
	return false;
}

void SPDtestWindow::regenerateProgressionTree()
{
	heatmapButton->setEnabled(TRUE);
	if( selection && this->HeatmapWin)
	{
		std::cout<< "rerender progression view"<<endl;
		selection->clear();
		std::vector< std::vector< long int> > sampleIndex;
		selection->GetSampleIndex( sampleIndex);
		
		vnl_matrix<double> clusAverageMat;
		std::vector< double> colorVec;
		std::vector< double> percentVec;
		SPDModel->GetSingleLinkageClusterAverage(sampleIndex, clusAverageMat);

		int maxId = this->maxVetexIdEdit->value();
		SPDModel->SetMaxVertexID(maxId);
		SPDModel->GetPercentage(sampleIndex, colorVec);
		std::string distanceThres = this->distanceThres->text().toStdString();
		SPDModel->GetCloseToDevicePercentage(sampleIndex, percentVec, atof(distanceThres.c_str()));

		std::vector<int> clusterNum;
		this->HeatmapWin->GetSubTreeClusterNum(clusterNum);
		//SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", selFeatureID);
		vtkSmartPointer<vtkTable> newtable = SPDModel->GenerateMST( clusAverageMat, selFeatureID, clusterNum);

		/** graph window set models */
		std::vector<int> index;
		SPDModel->GetSingleLinkageClusterMapping(sampleIndex, index);
		vtkSmartPointer<vtkTable> dataTable = vtkSmartPointer<vtkTable>::New();
		SPDModel->GetCombinedDataTable(dataTable);
		this->graph->setModels(dataTable, selection, &index);

		std::vector<std::string> headers;
		SPDModel->GetTableHeaders( headers);
		this->graph->SetTreeTable( newtable, headers[0], headers[1], headers[2], &colorVec, &percentVec);
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

void SPDtestWindow::ReColorProgressionTree(int nfeature)
{
	if( this->graph)
	{
		std::vector< std::vector< long int> > clusIndex;
		vnl_vector<double> featureValue;
		std::string featureName;
		selection->GetClusterIndex( clusIndex);
		SPDModel->GetClusterFeatureValue(clusIndex, nfeature, featureValue, featureName);
		this->graph->ColorTreeAccordingToFeatures(featureValue, featureName.c_str());
	}
}

void SPDtestWindow::updateSelMod()   
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

	std::vector<unsigned int> selMod;
	if( num <= size)   
	{
		selMod.resize(num);
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
		//searchSubsetsButton->setEnabled(TRUE);
		updateConnectedNumButton->setEnabled(TRUE);
		heatmapButton->setEnabled(FALSE);
	}
}

void SPDtestWindow::UpdateConnectedNum()
{
	std::string selectModulesID = this->psdModuleSelectBox->text().toStdString();
	std::vector< unsigned int> selModuleID;
	split( selectModulesID, ',', selModuleID);
	bool bchanged = false;
	for(int i = 0; i < selModuleID.size(); i++)
	{
		if(m_selModuleID.find(selModuleID[i]) == m_selModuleID.end())
		{
			bchanged = true;
		}
	}

	if( bchanged)
	{
		std::cout<< "Update Connected Component..."<<std::endl;
		m_selModuleID.clear();
		for(int i = 0; i < selModuleID.size(); i++)
		{
			m_selModuleID.insert(selModuleID[i]);
		}
		selFeatureID.clear();
		SPDModel->GetFeatureIdbyModId(selModuleID, selFeatureID);

		connectedNum = this->SPDModel->GetConnectedComponent(selFeatureID, connectedComponent);
		QString str = QString::number(connectedNum);
		connectedGraphEdit->setText(str);
	}
}

void SPDtestWindow::GetProgressionTreeOrder(std::vector<long int> &order)
{
	this->graph->GetProgressionTreeOrder(order);
}

void SPDtestWindow::showProgressionHeatmap()
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

	std::vector< std::vector< long int> > clusIndexByTreeOrder;
	for( int i = 0; i < TreeOrder.size(); i++)
	{
		clusIndexByTreeOrder.push_back(clusIndex[ TreeOrder[i] ]);
	}

	//vtkSmartPointer<vtkTable> tableForHistPlot = SPDModel->GetTableForHist(selOrder, unselOrder);
	//if( histo)
	//{
	//	delete histo;
	//}
	//histo = new HistoWindow(this);
	//histo->setModels(tableForHistPlot, selection, &clusIndexByTreeOrder);
	//histo->show();

	// progression heatmap
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);
	this->progressionHeatmap->setModelsforSPD( tableAfterCellCluster, selection, clusterOrder, selOrder, unselOrder, &indexMap);
	this->progressionHeatmap->showGraphforSPD( selOrder.size(), unselOrder.size(), true);
}

void SPDtestWindow::closeEvent(QCloseEvent *event)
{
	closeSubWindows();
	event->accept();
}

void SPDtestWindow::closeSubWindows()
{
	if(originalHeatmap)
	{
		originalHeatmap->close();
		delete originalHeatmap;
		originalHeatmap = NULL;
	}
	if(graph)
	{
		graph->close();
		delete graph;
		graph = NULL;
	}
	if(simHeatmap)
	{
		simHeatmap->close();
		delete simHeatmap;
		simHeatmap = NULL;
	}
	if(progressionHeatmap)
	{
		progressionHeatmap->close();
		delete progressionHeatmap;
		progressionHeatmap = NULL;
	}
	if(HeatmapWin)
	{
		HeatmapWin->close();
		delete HeatmapWin;
		HeatmapWin = NULL;
	}
	if(plot)
	{
		plot->close();
		delete plot;
		plot = NULL;
	}
	if(histo)
	{
		histo->close();
		delete histo;
		histo = NULL;
	}
}

void SPDtestWindow::ReRunSPDAnlysis()
{
	closeSubWindows();
	std::set<long int> selItems = this->selection->getSelections();
	disconnect(selection, SIGNAL( thresChanged()), this, SLOT( regenerateProgressionTree()));
	disconnect(selection, SIGNAL( ItemDeleted()), this, SLOT( ReRunSPDAnlysis()));
	this->selection->clear();
	vtkSmartPointer<vtkTable> table = GetSubTableExcludeItems(data, selItems);
	setModels( table, this->selection);
}

vtkSmartPointer<vtkTable> SPDtestWindow::GetSubTableExcludeItems(vtkSmartPointer<vtkTable> table, std::set<long int> &IDs)
{
	excludedIds = IDs;
	vtkSmartPointer<vtkTable> newTable = vtkSmartPointer<vtkTable>::New();
	for( vtkIdType i = 0; i < table->GetNumberOfColumns(); i++)
	{
		vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( table->GetColumnName(i));
		newTable->AddColumn(column);
    }

	for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
	{
		long int id = table->GetValue(i, 0).ToLong();
		if( IDs.find(id) == IDs.end())
		{
			newTable->InsertNextRow( table->GetRow(i));
		}
	}
	return newTable;
}

vtkSmartPointer<vtkTable> SPDtestWindow::NormalizeTable(vtkSmartPointer<vtkTable> table)
{
	SPDModel->ParseTraceFile( table, false);
	vtkSmartPointer<vtkTable> normalTable = SPDModel->GetDataTableAfterCellCluster();
	return normalTable;
}
