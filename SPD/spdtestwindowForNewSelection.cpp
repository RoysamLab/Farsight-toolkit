#include "spdtestwindowForNewSelection.h"
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

SPDWindowForNewSelection::SPDWindowForNewSelection(QWidget *parent) :
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

	sampleCoherenceLabel = new QLabel(tr("Sample Coherence(0.0 ~ 1.0):"));
	sampleCoherenceBox = new QDoubleSpinBox;
	sampleCoherenceBox->setValue( 0.9);	
	sampleCoherenceBox->setRange(0,1); 
	sampleCoherenceBox->setSingleStep(0.1);

    clusterCoherenceLabel = new QLabel(tr("Feature Coherence(0.0 ~ 1.0):"));
    clusterCoherenceBox = new QDoubleSpinBox;
	clusterCoherenceBox->setValue(0.9);
	clusterCoherenceBox->setRange(0,1); 
	clusterCoherenceBox->setSingleStep(0.1);

    kNearestNeighborLabel = new QLabel(tr("Feature Merge Coherence(0.0 ~ 1.0):"));
    kNearestNeighborBox = new QDoubleSpinBox;
	kNearestNeighborBox->setValue(0.9);
	kNearestNeighborBox->setRange(0,1);
	kNearestNeighborBox->setSingleStep(0.1);

    clusterButton = new QPushButton(tr("Feature Cluster"));
	cellClusterButton = new QPushButton(tr("Cell Cluster"));
	
	emdLabel = new QLabel(tr("Matching modules based on coherence:"));
	progressionOverDistance = new QLabel(tr("Progression over distance to device:"));
	bcheckBox = new QCheckBox();
	emdButton = new QPushButton(tr("Match"));

	emdThresBox = new QDoubleSpinBox;
	emdThresBox->setRange(0,1);
	emdThresBox->setSingleStep(0.1);

	emdPercentageBox = new QLineEdit;
	psmLable = new QLabel(tr("PSM Threshold(0.0 ~ 1.0):"));
	psmPerLable = new QLabel(tr("PSM Selected Blocks' Percentage:"));
    psmButton = new QPushButton(tr("Show PSM"));
	psmHisButton = new QPushButton(tr("PSM Histogram"));

	psdtLable = new QLabel(tr("Input hand-picked modules(seperate by comma):"));
	psdModuleSelectBox = new QLineEdit;
	maxVetexIdLabel = new QLabel(tr("Id to seperate:"));
	maxVetexIdEdit = new QSpinBox();
	maxVetexIdEdit->setSingleStep(100);
	maxVetexIdEdit->setRange(0,10000);
	maxVetexIdEdit->setValue(4500);

    psdtButton = new QPushButton(tr("View Progression"));
	heatmapLabel = new QLabel(tr("View Progression Heatmap:"));
	heatmapButton = new QPushButton(tr("Heatmap"));
	distanceThres = new QDoubleSpinBox;
	distanceThres->setRange(0,2000);
	distanceThres->setSingleStep(0.1);
	distanceThres->setValue(700.0);
	
	clusterButton->setEnabled(FALSE);
	cellClusterButton->setEnabled(FALSE);
	emdButton->setEnabled(FALSE);
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);
	bcheckBox->setChecked(FALSE);
	//saveFeatureButton = new QPushButton(tr("Save Selected Features"));

    connect(browseButton, SIGNAL(clicked()), this, SLOT(browse()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
	connect(loadTestButton, SIGNAL(clicked()), this, SLOT(loadContrastData()));
    connect(clusterButton, SIGNAL(clicked()), this, SLOT(clusterFunction()));
	connect(cellClusterButton , SIGNAL(clicked()), this, SLOT(clusterCells()));
	connect( bcheckBox, SIGNAL(clicked()), this, SLOT(updateProgressionType()));
	connect(emdButton, SIGNAL(clicked()), this, SLOT(emdFunction()));
	connect(psmButton, SIGNAL(clicked()), this, SLOT(showPSM()));
	connect(psmHisButton, SIGNAL(clicked()), this, SLOT(showPSMHist()));
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

    for ( int row = 1; row <= 14; row++)
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

	mainLayout->addWidget(sampleCoherenceLabel, 4, 0);
    mainLayout->addWidget(sampleCoherenceBox, 4, 1);
	//mainLayout->addWidget(cellClusterButton, 4, 2);

	mainLayout->addWidget(clusterCoherenceLabel, 5, 0);
    mainLayout->addWidget(clusterCoherenceBox, 5, 1);
	mainLayout->addWidget(clusterButton, 6, 2);

    mainLayout->addWidget(kNearestNeighborLabel, 6, 0);
    mainLayout->addWidget(kNearestNeighborBox, 6, 1);
	
	mainLayout->addWidget(progressionOverDistance, 7, 0);
	mainLayout->addWidget(bcheckBox, 7, 1);

	mainLayout->addWidget(emdLabel, 8, 0);
	mainLayout->addWidget(emdButton, 8, 2);

	mainLayout->addWidget(psmLable, 9, 0);
	mainLayout->addWidget(emdThresBox, 9, 1);

	mainLayout->addWidget(psmPerLable, 10, 0);
	mainLayout->addWidget(emdPercentageBox, 10, 1);
	mainLayout->addWidget(psmButton, 10, 2);

	mainLayout->addWidget(psdtLable, 11, 0);
	mainLayout->addWidget(distanceThres, 11, 1);
	mainLayout->addWidget(psdModuleSelectBox, 12, 0, 1, 2);
	mainLayout->addWidget(maxVetexIdLabel, 13, 0);
	mainLayout->addWidget(maxVetexIdEdit, 13, 1);
	mainLayout->addWidget(psdtButton, 13, 2);

	mainLayout->addWidget(heatmapLabel, 14, 0);
	mainLayout->addWidget(heatmapButton, 14, 2);

    setLayout(mainLayout);

	SPDModel = new SPDAnalysisModel();
}

SPDWindowForNewSelection::~SPDWindowForNewSelection()
{

}

void SPDWindowForNewSelection::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, SelectiveClustering * clusterSelection, ObjectSelection * sels2)
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

	if( clusterSelection == NULL)
	{
		clusterSelection = new SelectiveClustering();
	}
	else
	{
		ClusterSelections = clusterSelection;
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
		loadTestButton->setEnabled(TRUE);
	}
	else
	{
		SPDModel->ParseTraceFile( data);
		
		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		//this->SPDModel->NormalizeData();
		
		browseButton->setEnabled(FALSE);
		loadButton->setEnabled(FALSE);
		loadTestButton->setEnabled(TRUE);

		clusterButton->setEnabled(TRUE);
		cellClusterButton->setEnabled(TRUE);
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
	this->graph = new GraphWindowForNewSelection( this);


}

void SPDWindowForNewSelection::browse()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileName->setText(fileName);
        this->FileName = fileName;
    }
}

void SPDWindowForNewSelection::load()
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

void SPDWindowForNewSelection::loadContrastData()
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

void SPDWindowForNewSelection::clusterFunction()
{
	if ( this->SPDModel->GetFeatureNum() <= 0 && this->SPDModel->GetSampleNum() <= 0)
	{
		QMessageBox mes;
		mes.setText("You haven't loaded the data file!");
		mes.exec();
	}

	std::string clusterCor = this->clusterCoherenceBox->text().toStdString();
	std::string clusterMer = this->kNearestNeighborBox->text().toStdString();

	try
	{
		this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
		this->SPDModel->ClusterMerge( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
		emdButton->setEnabled(TRUE);
		psmButton->setEnabled(FALSE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "Clustering exception, please try again!"<<endl;
	}
}

void SPDWindowForNewSelection::clusterCells()
{
	std::string clusterCor = this->sampleCoherenceBox->text().toStdString();
	this->SPDModel->ClusterSamples(atof(clusterCor.c_str()));
	emdButton->setEnabled(TRUE);
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);
}

void SPDWindowForNewSelection::updateProgressionType()
{
	if( bcheckBox->isChecked())
	{
		bool rtn = this->SPDModel->SetProgressionType( true);   // progression over distance to device
		if( rtn == false)
		{
			QMessageBox mes;
			mes.setText("Distance to device is not available!");
			mes.exec();
		}
	}
	else
	{
		this->SPDModel->SetProgressionType( false);  // overall progression
	}
}

void SPDWindowForNewSelection::emdFunction()
{
	try
	{
		this->SPDModel->ModuleCoherenceMatchAnalysis();
		psmButton->setEnabled(TRUE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "EMD construction failure, please try again!"<<endl;
	}
}

void SPDWindowForNewSelection::editThreshold()
{
	QString emdThres = this->emdThresBox->text();
	double thres = this->emdThresBox->valueFromText(emdThres);
	double per = 0;
	if( thres >= 0 && thres <= 1)
	{
		per = this->SPDModel->GetCorMatSelectedPercentage( thres);
	}
	emdPercentageBox->setText(QString::number(per));
}

void SPDWindowForNewSelection::editPercentage()
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

void SPDWindowForNewSelection::showPSMHist()
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

void SPDWindowForNewSelection::showPSM()
{
	std::string emdThres = this->emdThresBox->text().toStdString();

	clusclus *clus1 = new clusclus();
	clusclus *clus2 = new clusclus();

	if( SPDModel->GetProgressionType())
	{
		std::vector< unsigned int> moduleIDs;
		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()), &moduleIDs);

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
		this->SPDModel->GetClusClusDataForCorMatrix(clus1, clus2, atof(emdThres.c_str()));
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

void SPDWindowForNewSelection::viewProgression()
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
	this->HeatmapWin = new HeatmapForNewSelection(this);

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
	connect(this->ClusterSelections, SIGNAL(selectionFinished()), this, SLOT(regenerateProgressionTree()));
	connect(selection, SIGNAL( ItemDeleted()), this, SLOT( ReRunSPDAnlysis()));
	connect(HeatmapWin, SIGNAL(columnToColorChanged(int)), this, SLOT( ReColorProgressionTree(int)));

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);
	this->HeatmapWin->setModelsforSPD( tableAfterCellCluster, selOrder, unselOrder, &indexMap, ClusterSelections);
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

void SPDWindowForNewSelection::split(std::string& s, char delim, std::vector< unsigned int>& indexVec)
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

void SPDWindowForNewSelection::GetFeatureOrder(std::vector< unsigned int> &selID, std::vector<int> &selIdOrder, std::vector<int> &unselIdOrder)
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

bool SPDWindowForNewSelection::IsExist(std::vector< unsigned int> vec, unsigned int value)
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

//void SPDWindowForNewSelection::regenerateProgressionTree()
//{
//	if( selection)
//	{
//		heatmapButton->setEnabled(TRUE);
//
//		std::string distanceThres = this->distanceThres->text().toStdString();
//		std::cout<< "rerender progression view"<<endl;
//		selection->clear();
//		std::vector< std::vector< long int> > sampleIndex;
//		selection->GetSampleIndex( sampleIndex);
//		
//		vnl_matrix<double> clusAverageMat;
//		std::vector<int> modSize;
//		std::vector< double> sampleVec;
//		std::vector< double> percentVec;
//		SPDModel->GetSingleLinkageClusterAverage(sampleIndex, clusAverageMat);
//		SPDModel->GetPercentage(sampleIndex, sampleVec);
//		SPDModel->GetCloseToDevicePercentage(sampleIndex, percentVec, atof(distanceThres.c_str()));
//
//		SPDModel->SaveSelectedFeatureNames("ReGenProgressionSelFeatures.txt", selFeatureID);
//		vtkSmartPointer<vtkTable> newtable = SPDModel->GenerateMST( clusAverageMat, selFeatureID);
//
//		/** graph window set models */
//		std::vector<int> index;
//		SPDModel->GetSingleLinkageClusterMapping(sampleIndex, index);
//		vtkSmartPointer<vtkTable> dataTable = vtkSmartPointer<vtkTable>::New();
//		SPDModel->GetCombinedDataTable(dataTable);
//		this->graph->setModels(dataTable, selection, &index);
//
//		std::vector<std::string> headers;
//		SPDModel->GetTableHeaders( headers);
//		this->graph->SetTreeTable( newtable, headers[0], headers[1], headers[2], &sampleVec, &percentVec);
//		try
//		{
//			this->graph->ShowGraphWindow();
//		}
//		catch(...)
//		{
//			std::cout<< "Graph window error!"<<endl;
//		}
//		
//	}
//}

void SPDWindowForNewSelection::regenerateProgressionTree()
{
	heatmapButton->setEnabled(TRUE);
	vtkSmartPointer<vtkTable> averageClusterTable = ClusterSelections->ClusterFeatureTable();
	vnl_matrix<double> clusAverageMat;
	std::vector<int> index; 
	vnl_vector<double> distance;
	SPDModel->ConvertTableToMatrix(averageClusterTable, clusAverageMat, index, distance);

	std::vector< std::vector< long int> > sampleIndex;
	std::vector< int> sizeVec;

	std::set< vtkIdType > clusIds = ClusterSelections->GetClusterIDs();
	std::set< vtkIdType >::iterator iter;

	for( iter = clusIds.begin(); iter != clusIds.end(); iter++)
	{
		std::set< vtkIdType > ids = ClusterSelections->SelectionFromCluster( *iter);
		std::set< vtkIdType >::iterator idIter;
		std::vector<long int> vecIds;
		for( idIter = ids.begin(); idIter != ids.end(); idIter++)
		{
			vecIds.push_back((long int)*idIter);
		}
		sampleIndex.push_back(vecIds);
		sizeVec.push_back(vecIds.size());
	}

	int maxId = this->maxVetexIdEdit->value();
	SPDModel->SetMaxVertexID(maxId);
	SPDModel->GetPercentage(sampleIndex, sampleVec);
	SPDModel->GetCloseToDevicePercentage(sampleIndex, percentVec, distanceThres->value());

	//std::vector<int> clusterNum;
	//this->HeatmapWin->GetSubTreeClusterNum(clusterNum);   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//SPDModel->SaveSelectedFeatureNames("SelFeatures.txt", selFeatureID);
	//vtkSmartPointer<vtkTable> newtable = SPDModel->GenerateMST( clusAverageMat, selFeatureID, clusterNum);

	//this->graph->setModels(averageClusterTable, ClusterSelections);

	//std::vector<std::string> headers;
	//SPDModel->GetTableHeaders( headers);
	//this->graph->SetTreeTable( newtable, headers[0], headers[1], headers[2], &sizeVec, &sampleVec, &percentVec);

	//try
	//{
	//	this->graph->ShowGraphWindow();
	//}
	//catch(...)
	//{
	//	std::cout<< "Graph window error!"<<endl;
	//}
}

void SPDWindowForNewSelection::ReColorProgressionTree(int nfeature)
{
	//if( this->graph)
	//{
	//	std::vector< std::vector< long int> > clusIndex;
	//	vnl_vector<double> featureValue;
	//	std::string featureName;
	//	selection->GetClusterIndex( clusIndex);
	//	SPDModel->GetClusterFeatureValue(clusIndex, nfeature, featureValue, featureName);
	//	this->graph->ColorTreeAccordingToFeatures(featureValue, featureName.c_str());
	//}
}

void SPDWindowForNewSelection::updateSelMod()   // possible bugs
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

void SPDWindowForNewSelection::GetProgressionTreeOrder(std::vector<long int> &order)
{
	this->graph->GetProgressionTreeOrder(order);
}

void SPDWindowForNewSelection::showProgressionHeatmap()
{
	//ofstream ofs("SPDHeatmapOptimalOrder.txt");
	if( this->progressionHeatmap)
	{
		delete this->progressionHeatmap;
	}
	this->progressionHeatmap = new HeatmapForNewSelection(this);
	
	std::vector<long int> TreeOrder;
	this->graph->GetProgressionTreeOrder(TreeOrder);   // order of the cluster 
	if( TreeOrder.size() <=0)
	{          
		std::cout<< "progression tree hasn't been built yet"<<endl;
		return;
	}

	std::vector< std::vector< long int> > clusIndex;
	std::vector< int> clusterOrder;
	std::set< vtkIdType > clusIds = ClusterSelections->GetClusterIDs();
	std::set< vtkIdType >::iterator iter = clusIds.begin();
	for(; iter != clusIds.end(); iter++)
	{
		std::vector< long int> ids = ClusterSelections->SelectionIDsFromCluster(*iter);
		clusIndex.push_back(ids);
	}

	SPDModel->GetClusterOrder(clusIndex, TreeOrder, clusterOrder);
	int maxId = this->maxVetexIdEdit->value();

	// module feature and percentage plot
	vtkSmartPointer<vtkTable> tableForAverModulePlot = SPDModel->GetAverModuleTable(clusIndex, TreeOrder, sampleVec, percentVec, selOrder, unselOrder, maxId);
	if( plot)
	{
		delete plot;
	}
	plot = new PlotWindow(this);
	plot->setModels(tableForAverModulePlot);
	plot->show();

	// progression heatmap
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();

	std::map< int, int> indexMap;
	SPDModel->GetClusterMapping(indexMap);
	this->progressionHeatmap->setModelsforSPD( tableAfterCellCluster, clusterOrder, selOrder, unselOrder, &indexMap, ClusterSelections);
	this->progressionHeatmap->showGraphforSPD( selOrder.size(), unselOrder.size(), true);
}

void SPDWindowForNewSelection::closeEvent(QCloseEvent *event)
{
	closeSubWindows();
	event->accept();
}

void SPDWindowForNewSelection::closeSubWindows()
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

void SPDWindowForNewSelection::ReRunSPDAnlysis()
{
	closeSubWindows();
	std::set<long int> selItems = this->selection->getSelections();
	disconnect(selection, SIGNAL( thresChanged()), this, SLOT( regenerateProgressionTree()));
	disconnect(selection, SIGNAL( ItemDeleted()), this, SLOT( ReRunSPDAnlysis()));
	this->selection->clear();
	vtkSmartPointer<vtkTable> table = GetSubTableExcludeItems(data, selItems);
	setModels( table, this->selection);
}

vtkSmartPointer<vtkTable> SPDWindowForNewSelection::GetSubTableExcludeItems(vtkSmartPointer<vtkTable> table, std::set<long int> &IDs)
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