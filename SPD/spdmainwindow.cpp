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
	selection = NULL;
	selection2 = NULL;

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

	sampleCoherenceLabel = new QLabel(tr("Sample Coherence(0.0 ~ 1.0):"));;
	sampleCoherenceBox = new QDoubleSpinBox;
	sampleCoherenceBox->setValue( 0.9);
	sampleCoherenceBox->setRange(0,1); 
	sampleCoherenceBox->setSingleStep(0.1);

    clusterCoherenceLabel = new QLabel(tr("Feature Coherence(0.0 ~ 1.0):"));
    clusterCoherenceBox = new QDoubleSpinBox;
	clusterCoherenceBox->setValue(0.9);
	clusterCoherenceBox->setRange(0,1); 
	clusterCoherenceBox->setSingleStep(0.1);

    clusterMergeLabel = new QLabel(tr("Feature Merge Coherence(0.0 ~ 1.0):"));
    clusterMergeBox = new QDoubleSpinBox;
	clusterMergeBox->setValue(0.9);
	clusterMergeBox->setRange(0,1);
	clusterMergeBox->setSingleStep(0.1);

    clusterButton = new QPushButton(tr("Feature Cluster"));
	cellClusterButton = new QPushButton(tr("Cell Cluster"));
	

	//listWidget = new QListWidget( this);
	mstLabel = new QLabel(tr("Generate MST for each module:"));
	generateMSTButton = new QPushButton(tr("MST"));
	emdLabel = new QLabel(tr("Matching MSTs with modules:"));
	emdButton = new QPushButton(tr("EMD"));

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
    psdtButton = new QPushButton(tr("View Progression"));
	heatmapLabel = new QLabel(tr("View Progression Heatmap:"));
	heatmapButton = new QPushButton(tr("Heatmap"));

	clusterButton->setEnabled(FALSE);
	cellClusterButton->setEnabled(FALSE);
	generateMSTButton->setEnabled(FALSE);
	emdButton->setEnabled(FALSE);
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);

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
	connect(heatmapButton, SIGNAL(clicked()), this, SLOT(showProgressionHeatmap()));
	
    QGridLayout *mainLayout = new QGridLayout;

    for ( int col = 0; col<= 2; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }

    for ( int row = 1; row <= 13; row++)
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
	mainLayout->addWidget(cellClusterButton, 4, 2);

	mainLayout->addWidget(clusterCoherenceLabel, 5, 0);
    mainLayout->addWidget(clusterCoherenceBox, 5, 1);
	mainLayout->addWidget(clusterButton, 6, 2);

    mainLayout->addWidget(clusterMergeLabel, 6, 0);
    mainLayout->addWidget(clusterMergeBox, 6, 1);
	
	mainLayout->addWidget(mstLabel, 7, 0);
    mainLayout->addWidget(generateMSTButton, 7, 2);
	mainLayout->addWidget(emdLabel, 8, 0);
	mainLayout->addWidget(emdButton, 8, 2);

	mainLayout->addWidget(psmLable, 9, 0);
	mainLayout->addWidget(emdThresBox, 9, 1);
	//mainLayout->addWidget(psmHisButton, 9, 2);

	mainLayout->addWidget(psmPerLable, 10, 0);
	mainLayout->addWidget(emdPercentageBox, 10, 1);
	mainLayout->addWidget(psmButton, 10, 2);

	mainLayout->addWidget(psdtLable, 11, 0);
	mainLayout->addWidget(psdModuleSelectBox, 12, 0, 1, 2);
	mainLayout->addWidget(psdtButton, 12, 2);

	mainLayout->addWidget(heatmapLabel, 13, 0);
	mainLayout->addWidget(heatmapButton, 13, 2);

    setLayout(mainLayout);

	SPDModel = SPDAnalysisModel::InitInstance();

	assert(SPDModel!=NULL);
	this->HeatmapWin = new Heatmap();
	graph =  new GraphWindow(this);
	simHeatmap = new ProgressionHeatmap(this);
	histo = new HistoWindow(this);
	this->progressionHeatmap = NULL;
	connect( simHeatmap, SIGNAL( SelChanged()), this, SLOT( updateSelMod()));


}

SPDMainWindow::~SPDMainWindow()
{
	SPDAnalysisModel::DeInstance();
	if( this->progressionHeatmap)
	{
		delete this->progressionHeatmap;
	}
}

void SPDMainWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
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
	connect(selection, SIGNAL( thresChanged()), this, SLOT( regenerateProgressionTree()));

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
		data->RemoveColumnByName("Soma_X");
		data->RemoveColumnByName("Soma_Y");
		data->RemoveColumnByName("Soma_Z");
		data->RemoveColumnByName("Trace_File");

		data->RemoveColumnByName("Soma_X_Pos");
		data->RemoveColumnByName("Soma_Y_Pos");
		data->RemoveColumnByName("Soma_Z_Pos");
		
		SPDModel->ParseTraceFile( data);
		std::cout<<"table size after parse "<<data->GetNumberOfRows()<<"\t"<<data->GetNumberOfColumns()<<endl;
		

		this->featureNum->setText( QString::number(this->SPDModel->GetFeatureNum()));
		this->sampleNum->setText( QString::number(this->SPDModel->GetSampleNum()));
		//this->SPDModel->NormalizeData();
		
		browseButton->setEnabled(FALSE);
		loadButton->setEnabled(FALSE);
		loadTestButton->setEnabled(FALSE);

		clusterButton->setEnabled(TRUE);
		cellClusterButton->setEnabled(TRUE);
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
		data = this->SPDModel->GetDataTable();
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
		//this->SPDModel->NormalizeData();
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

	try
	{
		this->SPDModel->ClusterAgglomerate( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
		this->SPDModel->ClusterMerge( atof(clusterCor.c_str()), atof(clusterMer.c_str()));
		generateMSTButton->setEnabled(TRUE);
		emdButton->setEnabled(FALSE);
		psmButton->setEnabled(FALSE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "Clustering exception, please try again!"<<endl;
	}
}

void SPDMainWindow::generateMST()
{
	try
	{
		this->SPDModel->GenerateMST();
		emdButton->setEnabled(TRUE);
		psmButton->setEnabled(FALSE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "MST construction failure, please try again!"<<endl;
	}
}

void SPDMainWindow::clusterCells()
{
	std::string clusterCor = this->sampleCoherenceBox->text().toStdString();
	this->SPDModel->ClusterCells(atof(clusterCor.c_str()));
	generateMSTButton->setEnabled(TRUE);
	emdButton->setEnabled(FALSE);
	psmButton->setEnabled(FALSE);
	psdtButton->setEnabled(FALSE);
	heatmapButton->setEnabled(FALSE);
}

void SPDMainWindow::emdFunction()
{
	try
	{
		this->SPDModel->RunEMDAnalysis();
		psmButton->setEnabled(TRUE);
		psdtButton->setEnabled(FALSE);
		heatmapButton->setEnabled(FALSE);
	}
	catch(...)
	{
		std::cout<< "EMD construction failure, please try again!"<<endl;
	}
}

void SPDMainWindow::editThreshold()
{
	QString emdThres = this->emdThresBox->text();
	double thres = this->emdThresBox->valueFromText(emdThres);
	double per = 0;
	if( thres >= 0 && thres <= 1)
	{
		per = this->SPDModel->GetEMDSelectedPercentage( thres);
	}
	emdPercentageBox->setText(QString::number(per));
}

void SPDMainWindow::editPercentage()
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

	clusclus clus1, clus2;
	this->SPDModel->GetClusClusData(clus1, clus2, atof(emdThres.c_str()));
	optimalleaforder.set_size(clus1.num_samples);
	for( int i = 0; i < clus1.num_samples; i++)
	{
		optimalleaforder[i] = clus1.optimalleaforder[i];
	}
	this->simHeatmap->setModels();
	this->simHeatmap->setDataForSimilarMatrixHeatmap(clus1.features, clus1.optimalleaforder, clus2.optimalleaforder, clus1.num_samples, clus2.num_samples);	
	this->simHeatmap->creatDataForSimilarMatrixHeatmap();
	this->simHeatmap->showSimilarMatrixGraph();
}

void SPDMainWindow::viewProgression()
{
	/* needs to be changed:
	   get orders of the features
	   setmodels heatmap: table after dimension reduced in sample space
	   selection: threshold selection, node selection
	   show default heatmap and tree here
	*/
	
	/** heatmap set models */
	std::string selectModulesID = this->psdModuleSelectBox->text().toStdString();
	std::vector< unsigned int> selModuleID;
	std::vector< int> selOrder;
	std::vector< int> unselOrder;
	std::vector< int> clusterSize;
	selFeatureID.clear();

	split( selectModulesID, ',', selModuleID);
	SPDModel->GetFeatureIdbyModId(selModuleID, selFeatureID);
	GetFeatureOrder( selFeatureID, selOrder, unselOrder);
	
	vtkSmartPointer<vtkTable> tableAfterCellCluster = SPDModel->GetDataTableAfterCellCluster();

	this->HeatmapWin->setModelsforSPD( tableAfterCellCluster, selection, selOrder, unselOrder);
	this->HeatmapWin->showGraphforSPD( selOrder.size());

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

void SPDMainWindow::split(std::string& s, char delim, std::vector< unsigned int>& indexVec)
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

void SPDMainWindow::GetFeatureOrder(std::vector< unsigned int> &selID, std::vector<int> &selIdOrder, std::vector<int> &unselIdOrder)
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

bool SPDMainWindow::IsExist(std::vector< unsigned int> vec, unsigned int value)
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

void SPDMainWindow::regenerateProgressionTree()
{
	if( selection)
	{
		std::cout<< "rerender progression view"<<endl;
		std::vector< std::vector< long int> > clusIndex;
		selection->GetClusterIndex( clusIndex);

		vnl_matrix<double> clusAverageMat;
		std::vector<int> modSize;

		SPDModel->GetSingleLinkageClusterAverage(clusIndex, clusAverageMat);
		SPDModel->GetSingleLinkageClusterModuleSize(clusIndex, modSize);

		SPDModel->SaveSelectedFeatureNames("ReGenProgressionSelFeatures.txt", selFeatureID);
		vtkSmartPointer<vtkTable> newtable = SPDModel->GenerateMST( clusAverageMat, selFeatureID);

		/** graph window set models */
		std::vector<int> index;
		SPDModel->GetSingleLinkageClusterMapping(clusIndex, index);

		this->graph->setModels(data, selection, &index);

		std::vector<std::string> headers;
		SPDModel->GetTableHeaders( headers);
		this->graph->SetTreeTable( newtable, headers[0], headers[1], headers[2]);
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

void SPDMainWindow::GetProgressionTreeOrder(std::vector<long int> &order)
{
	this->graph->GetProgressionTreeOrder(order);
}

void SPDMainWindow::showProgressionHeatmap()
{
	ofstream ofs("SPDHeatmapOptimalOrder.txt");
	if( this->progressionHeatmap)
	{
		delete this->progressionHeatmap;
	}

	this->progressionHeatmap = new ProgressionHeatmap(this);

	this->progressionHeatmap->setModels(data,selection,selection2);
	
	std::vector<long int> TreeOrder;
	this->graph->GetProgressionTreeOrder(TreeOrder);
	//SPDModel->GetDistanceOrder(TreeOrder);

	if( TreeOrder.size() <=0)
	{
		std::cout<< "progression tree hasn't been built yet"<<endl;
		return;
	}

	int *order = new int[TreeOrder.size()];
	ofs<< "Sample Tree Order:"<<endl;
	for( long int i = 0; i < TreeOrder.size(); i++)
	{
		order[i] = TreeOrder[i];
		ofs<< TreeOrder[i]<<"\t";
	}
	ofs<<endl;

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
	
	ofs<< "feature optimal order:"<<endl;
	for( int i = 0; i < cc2->num_samples; i++)
	{
		ofs<< cc2->optimalleaforder[i]<<"\t";
	}
	ofs<<endl;
	ofs.close();
	
	vnl_matrix<double> mat;
	SPDModel->GetMatrixData(mat);
	for(int i = 0; i < mat.cols(); i++)
	{
		vnl_vector<double> coln = mat.get_column(i);
		//double max = abs( coln.max_value());
		double interval = coln.max_value() - coln.min_value();
		if( interval != 0)
		{
			coln = ( coln - coln.min_value()) / interval;
			//coln = coln + 1;
			//coln = coln / max;
			for( int j = 0; j < coln.size(); j++)
			{
				if(coln[j] < 0.5)
				{
					coln[j] = 2 * pow( coln[j], 2);
				}
				else
				{
					coln[j] = - 2 * pow( coln[j] - 1, 2) + 1;
				}
			}		
		}
		mat.set_column(i, coln);
	}
	
	// optimal order is the mst tree order
	this->progressionHeatmap->setDataForHeatmap(mat.data_array(), order, cc2->optimalleaforder, TreeOrder.size(), cc2->num_samples);
	this->progressionHeatmap->setDataForDendrograms(NULL, cc2->treedata);
	this->progressionHeatmap->creatDataForProgressionHeatmap(1);	
	this->progressionHeatmap->showSPDGraph();

	for( int i = 0; i < FeatureTreeData.size(); i++)
	{
		delete ftreedata[i];
	}
	delete ftreedata;

	delete cc2;
	delete order;
	
	cout<<"finish simHeatmap..."<<endl;
}

void SPDMainWindow::closeEvent(QCloseEvent *event)
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
	event->accept();
}

