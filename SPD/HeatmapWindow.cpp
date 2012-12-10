#include "HeatmapWindow.h"

#define pi 3.1415926
#define POWER_PARAM 0.2

Heatmap::Heatmap(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;
	this->view = NULL;
	this->theme = NULL;
	this->graph_Layout = NULL;
	this->aPlane = NULL;
	this->cellData = NULL;
	this->celllut = NULL;
	this->mapper = NULL;
	this->actor = NULL;

	this->v = NULL;
	this->points = NULL;
	this->vertexColors = NULL;
	this->vetexlut = NULL;

	this->myCellPicker = NULL;

	this->ids1 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->ids2 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->ids1->SetNumberOfComponents(1);
	this->ids2->SetNumberOfComponents(1);	
	
	this->dencolors1 = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->denpoints1 = vtkSmartPointer<vtkPoints>::New();
	this->denlines1 = vtkSmartPointer<vtkCellArray>::New();
	this->denlinesPolyData1 =vtkSmartPointer<vtkPolyData>::New();
	this->denmapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->denactor1 = vtkSmartPointer<vtkActor>::New();

	this->dencolors2 = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->denpoints2 = vtkSmartPointer<vtkPoints>::New();
	this->denlines2 = vtkSmartPointer<vtkCellArray>::New();
	this->denlinesPolyData2 =vtkSmartPointer<vtkPolyData>::New();
	this->denmapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->denactor2 = vtkSmartPointer<vtkActor>::New();	

	this->removeActorflag = 0;
	this->denResetflag1 = 0;
	this->denResetflag2 = 0;
	this->continueselectnum = 0;
	this->continueselect = false;
	this->intersectionselect = false;
	this->clusflag = false;

	this->mapdata = NULL;
	this->Optimal_Leaf_Order1 = NULL;
	this->Optimal_Leaf_Order2 = NULL;
	this->connect_Data_Tree1 = NULL;
	this->connect_Data_Tree2 = NULL;
	this->ftreedata = NULL;
}

Heatmap::~Heatmap()
{
	if(this->mapdata)
	{
		for(int i=0; i<num_samples; i++)
			delete this->mapdata[i];
		delete this->mapdata;
	}

	if(this->Optimal_Leaf_Order1)
		delete this->Optimal_Leaf_Order1;
	if(this->Optimal_Leaf_Order2)
		delete this->Optimal_Leaf_Order2;

	if(this->connect_Data_Tree1)
		for(int i = 0; i<this->num_samples - 1; i++)
			delete this->connect_Data_Tree1[i];
	delete this->connect_Data_Tree1;

	if(this->connect_Data_Tree2)
		for(int i = 0; i<this->num_features - 1; i++)
			delete this->connect_Data_Tree2[i];
	delete this->connect_Data_Tree2;

	if(ftreedata)
	{
		for( int i = 0; i < this->table->GetNumberOfRows() - 1; i++)
		{
			delete ftreedata[i];
		}
		delete ftreedata;
	}
}

void Heatmap::setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features)
{
	this->num_samples = num_samples;
	this->num_features = num_features;

	this->rowMapFromOriginalToReorder.clear();
	this->columnMapFromOriginalToReorder.clear();

	this->mapdata = new double*[num_samples];
	for(int i=0; i<num_samples; i++)
	{
		this->mapdata[i] = new double[num_features];
		for(int j = 0 ; j<num_features; j++)
			this->mapdata[i][j] = features[i][j];
	}

	this->Optimal_Leaf_Order1 = new int[num_samples] ;
	for(int i=0; i<num_samples; i++)
	{
		this->Optimal_Leaf_Order1[i] = optimalleaforder1[i];
		this->rowMapFromOriginalToReorder.insert( std::pair< int, int>(optimalleaforder1[i], i));
	}

	this->Optimal_Leaf_Order2 = new int[num_features];
	for(int i=0; i<num_features; i++)
	{
		this->Optimal_Leaf_Order2[i] = optimalleaforder2[i];
		this->columnMapFromOriginalToReorder.insert( std::pair< int, int>(optimalleaforder2[i], i));
	}
}

void Heatmap::creatDataForHeatmap(double powCof)
{
	//double** mustd = new double*[2];
	//mustd[0] = new double[107];
	//mustd[1] = new double[107];

	//this->readmustd(mustd);
	//this->scaleData(mustd);
	this->scaleData();

	std::vector< double > temp;
	temp.resize(num_features);

	double** tempdata;
	tempdata = new double*[this->num_samples];
	for(int i = 0; i < this->num_samples; i++)
		tempdata[i] = new double[this->num_features];

	for(int i = 0; i < this->num_samples; i++)
	{
		double mean = 0.0; 
		double std = 0.0;
		double sum = 0.0;
		for(int j = 0; j < this->num_features; j++)
		{
			temp[j] = mapdata[i][Optimal_Leaf_Order2[j]];
		}

		for(int j = 0; j < this->num_features; j++)
			tempdata[i][j] = temp[j];

	}
	for(int i = 0; i < this->num_samples; i++)
		mapdata[this->num_samples - i - 1] = tempdata[Optimal_Leaf_Order1[i]]; 

	//const char* filename = "mapdata.txt";
	//FILE *fp1 = fopen(filename,"w");
	//for(int i=0; i<num_samples; i++)
	//{
	//	for(int j=0; j<num_features; j++)
	//		fprintf(fp1,"%f\t",mapdata[i][j]);
	//	fprintf(fp1,"\n");
	//}
	//fclose(fp1);

	if( this->connect_Data_Tree1 != NULL)
	{
		this->createDataForDendogram1(powCof);
	}
	if( this->connect_Data_Tree2 != NULL)
	{
		this->createDataForDendogram2(powCof);
	}
	else
	{
		this->createDataForDendogram2();
	}
}

void Heatmap::scaleData()
{
	for(int i = 0; i<this->num_features; i++)
	{
		double mean = 0.0;
		double std = 0.0;
		double sum = 0.0;

		for(int j = 0; j<this->num_samples; j++)
			mean += mapdata[j][i];

		mean /= this->num_samples;

		for(int j = 0; j<this->num_samples; j++)
			sum += (mapdata[j][i] - mean) * (mapdata[j][i] - mean);

		std = sqrt(sum/this->num_samples);

		if(std)
			for(int j = 0; j<this->num_samples; j++)
				mapdata[j][i] = (mapdata[j][i] - mean)/std;
		else
			for(int j = 0; j<this->num_samples; j++)
				mapdata[j][i] = 0;

	}
}

void Heatmap::scaleData(double** mustd)
{
	for(int i = 0; i<this->num_features; i++)
	{
		if(mustd[1][i+1])
			for(int j = 0; j<this->num_samples; j++)
				mapdata[j][Optimal_Leaf_Order2[i]] = (mapdata[j][Optimal_Leaf_Order2[i]] - mustd[0][i+1])/mustd[1][i+1];
		else
			for(int j = 0; j<this->num_samples; j++)
				mapdata[j][i] = 0;

	}
}

void Heatmap::readmustd(double** mustd)
{
	double temp[107];
	const int MAXLINESIZE = 10024;
	char line[MAXLINESIZE];

	ifstream infile;
	infile.open("mustd.txt");
	
	if(infile.is_open())
	{
		int i = 0;
		for(int t = 0 ;t<2; t++)
		{
			infile.getline(line, MAXLINESIZE);
			char* pch = strtok(line, "\t");

			int k = 0;
			while(pch)
			{
				temp[k++] = atof(pch);
				pch = strtok(NULL, "\t");
			}
			for(int j = 0; j<107; j++)
				mustd[i][j] = temp[j];
			i++;
		}
	}
	infile.close();

	for(int  i=0; i<2 ;i++)
	{
		for(int j = 0; j<107; j++)
			cout<<mustd[i][j]<<"\t";
		cout<<endl;
	}
}

void Heatmap::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
{
	this->table = table;
	this->indMapFromVertexToInd.clear();
	this->indMapFromIndToVertex.clear();
	if( this->table)
	{
		for( int i = 0; i < this->table->GetNumberOfRows(); i++)
		{
			int var = this->table->GetValue( i, 0).ToInt();
			this->indMapFromVertexToInd.insert( std::pair< int, int>(var, i));
			this->indMapFromIndToVertex.push_back( var);
		}
	}

	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;

	if(!sels2)
		this->Selection2 = new ObjectSelection();
	else
		this->Selection2 = sels2;

	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));
}

void Heatmap::runClusclus()
{
	double** datas;
	vtkVariant temp; 

	datas = new double*[this->table->GetNumberOfRows()];

	std::cout<<"number of rows"<<this->table->GetNumberOfRows()<<endl;
	std::cout<<"number of columns"<<this->table->GetNumberOfColumns()<<endl;

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		datas[i] = new double[this->table->GetNumberOfColumns() - 1 + 2 ];
	}

	for(int i = 0; i < this->table->GetNumberOfRows(); i++)
	{		
		for(int j = 1; j < this->table->GetNumberOfColumns(); j++)
		{
			temp = this->table->GetValue(i, j);
			datas[i][j-1] = temp.ToDouble();
		}
	}

	cc1 = new clusclus(datas, (int)this->table->GetNumberOfRows(), (int)this->table->GetNumberOfColumns() - 1);
	cc1->Transpose();
	cc2 = new clusclus(cc1->transposefeatures,cc1->num_features, cc1->num_samples);

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			cc1->RunClusClus();
			/*cc1->WriteClusteringOutputToFile("mergers.txt","features.txt","progress.txt", "members.txt",
				"gap.txt", "treedata.txt", "Optimalleaforder.txt");*/
		}
		#pragma omp section
		{
			cc2->RunClusClus();
			/*cc2->WriteClusteringOutputToFile("mergers2.txt","features2.txt","progress2.txt", "members2.txt",
				"gap2.txt", "treedata2.txt", "Optimalleaforder2.txt");*/
		}
	}

	cout<<"finish clusclus....."<<endl;
	this->setDataForHeatmap(cc1->features, cc1->optimalleaforder, cc2->optimalleaforder,cc1->num_samples, cc2->num_samples);
	this->setDataForDendrograms(cc1->treedata, cc2->treedata);
	this->creatDataForHeatmap(POWER_PARAM);	

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		delete datas[i];
	}
	delete datas;

	delete cc1;
	delete cc2;
}

void Heatmap::runClus()
{
	this->clusflag = true;
	double** datas;
	vtkVariant temp; 

	datas = new double*[this->table->GetNumberOfRows()];

	std::cout<<this->table->GetNumberOfRows()<<endl;
	std::cout<<this->table->GetNumberOfColumns()<<endl;

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		datas[i] = new double[this->table->GetNumberOfColumns() - 1 + 2 ];
	}

	for(int i = 0; i < this->table->GetNumberOfRows(); i++)
	{		
		for(int j = 1; j < this->table->GetNumberOfColumns(); j++)
		{
			temp = this->table->GetValue(i, j);
			datas[i][j-1] = temp.ToDouble();
		}
	}

	cc1 = new clusclus(datas, (int)this->table->GetNumberOfRows(), (int)this->table->GetNumberOfColumns() - 1);
	cc1->RunClusClus();
	cout<<"finish clusclus....."<<endl;
	cout<<this->table->GetNumberOfRows();
	cout<<this->table->GetNumberOfColumns();
	cc1->WriteClusteringOutputToFile("mergers.txt","features.txt","progress.txt", "members.txt", "gap.txt", "treedata.txt", "Optimalleaforder.txt");

	int* optimalleaforder2 = new int[cc1->num_features];

	for(int i = 0;i<cc1->num_features; i++)
		optimalleaforder2[i]=i;
	this->setDataForHeatmap(cc1->features, cc1->optimalleaforder, optimalleaforder2,cc1->num_samples, cc1->num_features);
	this->setDataForDendrograms(cc1->treedata);
	this->creatDataForHeatmap(POWER_PARAM);
	
	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		delete datas[i];
	}
	delete datas;

	delete cc1;
}

void Heatmap::showGraph()
{	
	if(this->clusflag == true)
		this->drawPoints3();
	else
		this->drawPoints1();


	this->aPlane = vtkSmartPointer<vtkPlaneSource>::New();
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	this->cellData = vtkSmartPointer<vtkFloatArray>::New();

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	this->celllut = vtkSmartPointer<vtkLookupTable>::New();
	this->celllut->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->celllut->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->celllut->Build();

	//int k = 0;
	//boost::math::normal N;
	//try
	//{
	//	for(int i = 0; i < this->num_samples; i++)
	//	{
	//		for(int j = 0; j < this->num_features; j++)
	//		{
	//			if(mapdata[num_samples - i - 1][j] < 0)
	//				celllut->SetTableValue(k++, 0, 1 - cdf(N,mapdata[num_samples - i - 1][j]), 0);
	//			else if(mapdata[num_samples - i - 1][j] > 0)
	//				celllut->SetTableValue(k++, cdf(N,mapdata[num_samples - i - 1][j]), 0, 0);
	//			else
	//				celllut->SetTableValue(k++, 0, 0, 0);
	//		}
	//	}
	//}
	//catch(...)
	//{
	//	cout<<"Boost call failed! please try again!"<<endl;
	//}

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			rgb rgb = GetRGBValue( mapdata[num_samples - i - 1][j]);
			celllut->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->aPlane->Update();
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(celllut);

	this->actor = vtkSmartPointer<vtkActor>::New();
	this->actor->SetMapper(mapper);

	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetTableRange (-1, 1);
	scalarbarLut->SetNumberOfTableValues(COLOR_MAP_SIZE);
	for(int index = 0; index<COLOR_MAP_SIZE;index++)
	{
		rgb rgbscalar = COLORMAP[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(10);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->SetMaximumHeightInPixels(1000);
	scalarBar->SetMaximumWidthInPixels(100);

	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->AddActor2D(scalarBar);
	this->SetInteractStyle();
	this->view->GetRenderer()->GradientBackgroundOff();
	this->view->GetRenderer()->SetBackground(1,1,1);

	try
	{
		if(this->clusflag == true)
			this->showDendrogram1();
		else
		{
			this->showDendrogram1();
			this->showDendrogram2();
		}
	}
	catch(...)
	{
		cout<<"Draw dendrogram failed ! Please try again!"<<endl;
	}
	this->view->Render();
	this->view->GetInteractor()->Start();
}

void Heatmap::SetInteractStyle()
{
	this->theme->SetCellValueRange(0, this->num_samples*this->num_features - 1);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
	this->view->ApplyViewTheme(theme);

	this->myCellPicker = vtkSmartPointer<vtkCellPicker>::New();
	this->view->GetInteractor()->SetPicker(this->myCellPicker);
	this->myCellPicker->SetTolerance(0.004);
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2 =vtkSmartPointer<vtkCallbackCommand>::New();
	selectionCallback2->SetClientData(this);
	selectionCallback2->SetCallback(SelectionCallbackFunction2 );

	vtkSmartPointer<vtkCallbackCommand> selectionCallback3 =vtkSmartPointer<vtkCallbackCommand>::New();
	selectionCallback3->SetClientData(this);
	selectionCallback3->SetCallback(SelectionCallbackFunction3);

	this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
	this->keyPress->SetCallback(HandleKeyPress);
	this->keyPress->SetClientData(this);

	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->view->GetInteractor()->AddObserver(vtkCommand::RightButtonPressEvent, selectionCallback2);
	this->view->GetInteractor()->AddObserver(vtkCommand::RightButtonReleaseEvent, selectionCallback3);
	this->view->GetInteractor()->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);

	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
}

rgb Heatmap::GetRGBValue(double val)
{
	int index = 64 * (val+1) - 1;   // when val = 1; index should be the max index
	if( index >= COLOR_MAP_SIZE)
	{
		index = COLOR_MAP_SIZE - 1;
	}
	else if( index < 0)
	{
		index = 0;
	}
	return COLORMAP[index];
}

void Heatmap::setDataForDendrograms(double** treedata1, double** treedata2)
{
	if(treedata1 != NULL)
	{
		this->connect_Data_Tree1 = new double*[this->num_samples-1];
		this->rowMapForTreeData.clear();
		for(int i = 0; i<this->num_samples - 1; i++)
		{
			this->connect_Data_Tree1[i] = new double[4];
			for(int j = 0; j<4; j++)
				this->connect_Data_Tree1[i][j] = treedata1[i][j];

			this->rowMapForTreeData.insert( std::pair< int, int>(treedata1[i][3], i));

		}
	}

	if( treedata2 != NULL)
	{
		this->connect_Data_Tree2 = new double*[this->num_features-1];
		this->columnMapForTreeData.clear();
		for(int i = 0; i<this->num_features - 1; i++)
		{
			this->connect_Data_Tree2[i] = new double[4];
			for(int j = 0; j<4; j++)
				this->connect_Data_Tree2[i][j] = treedata2[i][j];

			this->columnMapForTreeData.insert( std::pair< int, int>(treedata2[i][3], i));
		}
	}
}

void Heatmap::createDataForDendogram1(double powCof)
{
	this->Processed_Coordinate_Data_Tree1.resize(2*(this->num_samples) - 1);
	for(int i = 0; i < 2*(this->num_samples) - 1; i++)
	{
		this->Processed_Coordinate_Data_Tree1[i].resize(4);
	}

	for(int i = 0; i < num_samples; i++)
	{
		Processed_Coordinate_Data_Tree1[i][0] = i;
		int k = rowMapFromOriginalToReorder.find(i)->second;
		Processed_Coordinate_Data_Tree1[i][2] = (k + 0.5)/(double)this->num_samples - 0.5;
		Processed_Coordinate_Data_Tree1[i][1] = -0.5;
		Processed_Coordinate_Data_Tree1[i][3] = 0; 
	}

	std::cout<<std::endl;
	std::cout<< "Max value:"<<connect_Data_Tree1[num_samples-2][2]<<std::endl;
	std::cout<< "devided by:"<<2 * pow(connect_Data_Tree1[num_samples - 2][2], powCof)<<std::endl;

	for(int i = 0; i < num_samples-1; i++)
	{
		connect_Data_Tree1[i][2] = pow(connect_Data_Tree1[i][2], powCof);
		connect_Data_Tree1[i][2] /= pow(connect_Data_Tree1[num_samples-2][2], powCof);
		connect_Data_Tree1[i][2] /= 2;
	}
	connect_Data_Tree1[num_samples-2][2] = 0.5;

	for(int i = num_samples ; i < 2*num_samples - 1; i++)
	{
		Processed_Coordinate_Data_Tree1[i][0] = i;

		for(int k = 0; k < num_samples -1 ; k++)
		{
			if(i == connect_Data_Tree1[k][3])
			{
				double temp1, temp2;
				temp1 = connect_Data_Tree1[k][0];
				temp2 = connect_Data_Tree1[k][1];
				Processed_Coordinate_Data_Tree1[i][2] = (Processed_Coordinate_Data_Tree1[temp1][2] + Processed_Coordinate_Data_Tree1[temp2][2])/2;
				Processed_Coordinate_Data_Tree1[i][1] = -connect_Data_Tree1[k][2] - 0.5;
			}
		}

		Processed_Coordinate_Data_Tree1[i][3] = 0; 
	}
}

void Heatmap::createDataForDendogram2()
{
	this->Processed_Coordinate_Data_Tree2.resize(this->num_features);
	for(int i = 0; i < this->num_features; i++)
	{
		this->Processed_Coordinate_Data_Tree2[i].resize(4);
	}

	for(int i = 0; i < num_features; i++)
	{
		Processed_Coordinate_Data_Tree2[i][0] = i;
		int k = columnMapFromOriginalToReorder.find(i)->second;
		Processed_Coordinate_Data_Tree2[i][1] = (k+0.5)/(double)this->num_features - 0.5;
		Processed_Coordinate_Data_Tree2[i][2] = 0.5;
		Processed_Coordinate_Data_Tree2[i][3] = 0; 
	}
}


void Heatmap::createDataForDendogram2(double powCof)
{
	this->Processed_Coordinate_Data_Tree2.resize(2*(this->num_features) - 1);
	for(int i = 0; i < 2*(this->num_features) - 1; i++)
	{
		this->Processed_Coordinate_Data_Tree2[i].resize(4);
	}

	for(int i = 0; i < num_features; i++)
	{
		Processed_Coordinate_Data_Tree2[i][0] = i;
		int k = columnMapFromOriginalToReorder.find(i)->second;
		Processed_Coordinate_Data_Tree2[i][1] = (k+0.5)/(double)this->num_features - 0.5;
		Processed_Coordinate_Data_Tree2[i][2] = 0.5;
		Processed_Coordinate_Data_Tree2[i][3] = 0; 
	}

	for(int i = 0; i < num_features-1; i++)
	{
		connect_Data_Tree2[i][2] = pow(connect_Data_Tree2[i][2], powCof);
		connect_Data_Tree2[i][2] /= pow(connect_Data_Tree2[num_features - 2][2], powCof);
		connect_Data_Tree2[i][2] /= 2;
	}

	connect_Data_Tree2[num_features - 2][2] = 0.5;

	for(int i = num_features ; i < 2*num_features - 1; i++)
	{
		Processed_Coordinate_Data_Tree2[i][0] = i;

		for(int k = 0; k < num_features -1 ; k++)
		{
			if(i == connect_Data_Tree2[k][3])
			{
				double temp1, temp2;
				temp1 = connect_Data_Tree2[k][0];
				temp2 = connect_Data_Tree2[k][1];
				Processed_Coordinate_Data_Tree2[i][1] = (Processed_Coordinate_Data_Tree2[temp1][1] + Processed_Coordinate_Data_Tree2[temp2][1])/2;
				Processed_Coordinate_Data_Tree2[i][2] = connect_Data_Tree2[k][2] + 0.5;
			}
		}

		Processed_Coordinate_Data_Tree2[i][3] = 0; 
	}
}
void Heatmap::showDendrogram1()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	
	this->dencolors1->SetNumberOfComponents(3);
	this->dencolors1->SetName("denColors1");
	unsigned char color[3] = {0, 0, 0};
	for(int i=0; i<3*(this->num_samples - 1);i++)
		this->dencolors1->InsertNextTupleValue(color);

	for(int i=0; i<this->num_samples-1;i++)
	{
		double temp1 = this->connect_Data_Tree1[i][0];
        double temp2 = this->connect_Data_Tree1[i][1];

		for(int j=0; j<(2*(this->num_samples))-1; j++)
        {
            if(this->Processed_Coordinate_Data_Tree1[j][0]==temp1)
			{
				p1[0]=this->Processed_Coordinate_Data_Tree1[j][1];
				p1[1]=this->Processed_Coordinate_Data_Tree1[j][2];
				p1[2]=this->Processed_Coordinate_Data_Tree1[j][3];
            }   
            if(this->Processed_Coordinate_Data_Tree1[j][0]==temp2)
			{
                p2[0]=this->Processed_Coordinate_Data_Tree1[j][1];
                p2[1]=this->Processed_Coordinate_Data_Tree1[j][2];
                p2[2]=this->Processed_Coordinate_Data_Tree1[j][3];
			}                             
        }
	   
        p3[0]=-connect_Data_Tree1[i][2] - 0.5;
        p3[1]=p1[1];
        p3[2]=p1[2];

        p4[0]=-connect_Data_Tree1[i][2] - 0.5;
        p4[1]=p2[1];
        p4[2]=p2[2];


		this->denpoints1->InsertNextPoint(p1);
		this->denpoints1->InsertNextPoint(p2);
		this->denpoints1->InsertNextPoint(p3);
		this->denpoints1->InsertNextPoint(p4);

		
		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0,0 + i*4);
		line0->GetPointIds()->SetId(1,2 + i*4);
		this->denlines1->InsertNextCell(line0);

		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0,1 + i*4);
		line1->GetPointIds()->SetId(1,3 + i*4);
		this->denlines1->InsertNextCell(line1);

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0,2 + i*4);
		line2->GetPointIds()->SetId(1,3 + i*4);
		this->denlines1->InsertNextCell(line2);
	}

	
	this->denlinesPolyData1->SetPoints(denpoints1);
	this->denlinesPolyData1->SetLines(denlines1);
	this->denlinesPolyData1->GetCellData()->SetScalars(dencolors1);

	
	this->denmapper1->SetInput(denlinesPolyData1);
	this->denmapper1->SetScalarRange(0, 3*this->num_samples-1);

	this->denactor1 = vtkSmartPointer<vtkActor>::New();
	this->denactor1->SetMapper(denmapper1);
	
	this->view->GetRenderer()->AddActor(denactor1);
}

void Heatmap::showDendrogram2()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	this->dencolors2->SetNumberOfComponents(3);
	this->dencolors2->SetName("denColors2");
	unsigned char color[3] = {0, 0, 0};
	for(int i=0; i<3*(this->num_features - 1);i++)
		this->dencolors2->InsertNextTupleValue(color);

	for(int i=0; i<this->num_features-1; i++)
	{
		double temp1 = this->connect_Data_Tree2[i][0];
        double temp2 = this->connect_Data_Tree2[i][1];

		for(int j=0; j<(2*(this->num_features))-1; j++)
        {
            if(this->Processed_Coordinate_Data_Tree2[j][0]==temp1)
			{
				p1[0]=this->Processed_Coordinate_Data_Tree2[j][1];
				p1[1]=this->Processed_Coordinate_Data_Tree2[j][2];
				p1[2]=this->Processed_Coordinate_Data_Tree2[j][3];
            }   
            if(this->Processed_Coordinate_Data_Tree2[j][0]==temp2)
			{
                p2[0]=this->Processed_Coordinate_Data_Tree2[j][1];
                p2[1]=this->Processed_Coordinate_Data_Tree2[j][2];
                p2[2]=this->Processed_Coordinate_Data_Tree2[j][3];
			}                             
        }
	   
        p3[0]=p1[0];
        p3[1]=this->connect_Data_Tree2[i][2] + 0.5;
        p3[2]=p1[2];

        p4[0]=p2[0];
        p4[1]=this->connect_Data_Tree2[i][2] + 0.5;
        p4[2]=p2[2];

		this->denpoints2->InsertNextPoint(p1);
		this->denpoints2->InsertNextPoint(p2);
		this->denpoints2->InsertNextPoint(p3);
		this->denpoints2->InsertNextPoint(p4);

		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0,0 + i*4);
		line0->GetPointIds()->SetId(1,2 + i*4);
		this->denlines2->InsertNextCell(line0);

		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0,1 + i*4);
		line1->GetPointIds()->SetId(1,3 + i*4);
		this->denlines2->InsertNextCell(line1);

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0,2 + i*4);
		line2->GetPointIds()->SetId(1,3 + i*4);
		this->denlines2->InsertNextCell(line2);
	}

	this->denlinesPolyData2->SetPoints(denpoints2);
	this->denlinesPolyData2->SetLines(denlines2);
	this->denlinesPolyData2->GetCellData()->SetScalars(dencolors2);

	this->denmapper2->SetInput(denlinesPolyData2);
	this->denmapper2->SetScalarRange(0, 3*this->num_features-1);
	
	this->denactor2->SetMapper(denmapper2);
	this->view->GetRenderer()->AddActor(denactor2);
}

void Heatmap::GetSelecectedIDs()
{
	cout<<"get selected"<<endl;
	std::set<long int> selectedIDs2 = this->Selection2->getSelections();
	std::set<long int> selectedIDs1 = this->Selection->getSelections();	
	std::set<long int>::iterator iter1 = selectedIDs1.begin();
	std::set<long int>::iterator iter2 = selectedIDs2.begin();
	vtkSmartPointer<vtkIdTypeArray> cellids = vtkSmartPointer<vtkIdTypeArray>::New();
	cellids->SetNumberOfComponents(1);

	int num1 = selectedIDs1.size();
	int num2 = selectedIDs2.size();

	std::vector<int > IDs1;
	std::vector<int > IDs2;
	IDs1.resize(num1);
	IDs2.resize(num2);

	int count1 = 0;
	int count2 = 0;

	#pragma omp parallel sections
	{
		#pragma omp section
		while(iter1 != selectedIDs1.end())
		{
			int index1 = *iter1;
			int var = indMapFromVertexToInd.find(index1)->second;
			int id1 = rowMapFromOriginalToReorder.find(var)->second;
			IDs1[count1++] = id1;
			iter1++;
		}

		#pragma omp section 
		while(iter2 != selectedIDs2.end())
		{
			int index2 = *iter2;
			int id2 = columnMapFromOriginalToReorder.find(index2)->second;
			IDs2[count2++] = id2;
			iter2++;
		}
	}

	if( num1 == 0 && num2 != 0)
	{
		num1 = this->num_samples;
		IDs1.resize(num1);
		for( int i = 0; i < this->num_samples; i++)
		{
			IDs1[i] = i;
		}
	}

	if( num2 == 0 && num1 != 0)
	{
		num2 = this->num_features;
		IDs2.resize(num2);
		for( int i = 0; i < this->num_features; i++)
		{
			IDs2[i] = i;
		}
	}

	for(int i = 0; i<num1; i++)
		for(int j = 0; j<num2; j++)
			cellids->InsertNextValue( (IDs1[i])*(this->num_features) + IDs2[j]);

	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(cellids);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);
	 
	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInput(0, this->aPlane->GetOutput());
	extractSelection->SetInput(1, selection);
	extractSelection->Update();
	 
	vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());
	
	vtkSmartPointer<vtkDataSetMapper> selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selectedMapper->SetInputConnection(selected->GetProducerPort());
	 
	vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->GetProperty()->EdgeVisibilityOn();
	selectedActor->GetProperty()->SetEdgeColor(1,1,1);
	selectedActor->GetProperty()->SetLineWidth(0.5);

	if(this->dragLineFlag)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor);
	try
	{
		if(continueselect == false)
		{
			if(continueselectnum > 0)
			{
				cout<<"I'm here "<<continueselectnum<<endl;
				for(int i = 0; i<continueselectnum + 1; i++)
					this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
				continueselectnum = 0;
			}
			else
			{
				if (this->removeActorflag != 0)
					this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
			}

			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
			this->removeActorflag += 1;
		}
		else
		{
			this->continueselectnum += 1;
			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);

		}
		if(this->dragLineFlag)
			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor);
		this->view->Render();
	}
	catch(...)
	{
		cout<<"GetSelecectedIDs failed, please try it again!"<<endl;
	}
}

void Heatmap::drawPoints1()
{
	int max_table_values = 2*this->num_samples + 2*this->num_features-2;

	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->v = vtkSmartPointer<vtkIdTypeArray>::New();
	v->SetNumberOfValues (max_table_values);

	for(int i=0; i<2*this->num_samples-1;i++)
    {
		v->SetValue (i,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree1[i][1],this->Processed_Coordinate_Data_Tree1[i][2],this->Processed_Coordinate_Data_Tree1[i][3]);
	}
	
	for(int i=0; i<(2*this->num_features-1);i++)
    {
		v->SetValue (i+2*this->num_samples-1,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1],this->Processed_Coordinate_Data_Tree2[i][2],this->Processed_Coordinate_Data_Tree2[i][3]);
	}
    this->graph_Layout->SetPoints(this->points);

	this->vertexColors = vtkSmartPointer<vtkIntArray>::New();   
    vertexColors->SetNumberOfComponents(1);
    vertexColors->SetName("Color1");	

	this->vetexlut = vtkSmartPointer<vtkLookupTable>::New();
	vetexlut->SetNumberOfTableValues(max_table_values);
    for(int i=0; i<max_table_values;i++)
    {
		vetexlut->SetTableValue(i, 0.5, 0.5,0.5); // color the vertices- blue
    }  
    vetexlut->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(int j=0;j<max_table_values;j++)
    {
		vertexColors->InsertNextValue(j);
		scales1->InsertNextValue(1);
    }

	this->graph_Layout->GetVertexData()->AddArray(vertexColors);
	this->graph_Layout->GetVertexData()->AddArray(scales1);

	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
    this->view->AddRepresentationFromInput(graph_Layout);
    this->view->SetLayoutStrategy("Pass Through");	
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
	this->view->ColorVerticesOn();
	this->view->SetVertexColorArrayName("Color1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->theme = vtkSmartPointer<vtkViewTheme>::New(); 
	this->theme->SetPointLookupTable(vetexlut);

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	for(int i=0; i<max_table_values;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
	}
	for(int i = 0; i < this->num_samples; i++)
	{
		vtkVariant v = i;
		label->InsertNextValue(v.ToString());
		orient->InsertNextValue(0.0);
	}
	for(int i = this->num_samples; i <2*(this->num_samples)-1; i++)
	{
		label->InsertNextValue("");
		orient->InsertNextValue(0.0);
	}
	for(int i = 2*(this->num_samples)-1; i <2*(this->num_samples)-1 + this->num_features; i++)
	{	
		vtkIdType id = i - (2*this->num_samples-1) + 1 ;
		label->InsertNextValue(this->table->GetColumn(id)->GetName ());
		orient->InsertNextValue(45.0);
	}	
	for(int i = 2*(this->num_samples)-1 + this->num_features; i <2*(this->num_samples)-1 + 2*(this->num_features)-1; i++)
	{
		label->InsertNextValue("");
		orient->InsertNextValue(0.0);
	}

	pd->SetPoints(points);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);
    this->selectionCallback1->SetCallback ( SelectionCallbackFunction1);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void Heatmap::drawPoints3()
{
	int max_table_values = 2*this->num_samples-1;

	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->v = vtkSmartPointer<vtkIdTypeArray>::New();
	v->SetNumberOfValues (max_table_values);

	for(int i=0; i<max_table_values;i++)
    {
		v->SetValue (i,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree1[i][1],this->Processed_Coordinate_Data_Tree1[i][2],this->Processed_Coordinate_Data_Tree1[i][3]);
	}
	
    this->graph_Layout->SetPoints(this->points);
     
	this->vertexColors = vtkSmartPointer<vtkIntArray>::New();
    vertexColors->SetNumberOfComponents(1);
    vertexColors->SetName("Color1");
	
	this->vetexlut = vtkSmartPointer<vtkLookupTable>::New();
	vetexlut->SetNumberOfTableValues(max_table_values);
    for(int i=0; i<max_table_values;i++)
    {
		vetexlut->SetTableValue(i, 0.5, 0.5,0.5); 
    }  
    vetexlut->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(int j=0; j<max_table_values;j++)
    {
		vertexColors->InsertNextValue(j);
		scales1->InsertNextValue(1);
    }

	this->graph_Layout->GetVertexData()->AddArray(scales1);
	this->graph_Layout->GetVertexData()->AddArray(vertexColors);
	
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
    this->view->AddRepresentationFromInput(graph_Layout);
    this->view->SetLayoutStrategy("Pass Through");
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
	this->view->ColorVerticesOn();
	this->view->SetVertexColorArrayName("Color1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	this->theme->SetPointLookupTable(vetexlut);

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for(int i=0; i<this->num_features;i++)
    {
		pts->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1], 0.5, 0);
	}

	for(int i=0; i<this->num_features;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(45.0);
		vtkIdType id = i+1  ;
		label->InsertNextValue(this->table->GetColumn(id)->GetName());
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);
	this->selectionCallback1->SetCallback (SelectionCallbackFunction1);

    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void Heatmap::drawPointsForSPD()
{
	int max_table_values = 2 * this->num_samples-1 + this->num_features;

	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->v = vtkSmartPointer<vtkIdTypeArray>::New();
	v->SetNumberOfValues (max_table_values);

	for( int i = 0; i < 2 * this->num_samples - 1; i++)
    {
		v->SetValue(i, graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree1[i][1],this->Processed_Coordinate_Data_Tree1[i][2],this->Processed_Coordinate_Data_Tree1[i][3]);
	}
	for( int i = 0; i < this->num_features; i++)
	{
		v->SetValue(i + 2 * this->num_samples - 1, graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1],this->Processed_Coordinate_Data_Tree2[i][2],this->Processed_Coordinate_Data_Tree2[i][3]);
	}
	
    this->graph_Layout->SetPoints(this->points);
     
	this->vertexColors = vtkSmartPointer<vtkIntArray>::New();
    vertexColors->SetNumberOfComponents(1);
    vertexColors->SetName("Color1");
	
	this->vetexlut = vtkSmartPointer<vtkLookupTable>::New();
	vetexlut->SetNumberOfTableValues(max_table_values);
    for(int i=0; i<max_table_values;i++)
    {
		vetexlut->SetTableValue(i, 0.5, 0.5,0.5); 
    }  
    vetexlut->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(int j=0; j<max_table_values;j++)
    {
		vertexColors->InsertNextValue(j);
		scales1->InsertNextValue(1);
    }

	this->graph_Layout->GetVertexData()->AddArray(scales1);
	this->graph_Layout->GetVertexData()->AddArray(vertexColors);
	
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
    this->view->AddRepresentationFromInput(graph_Layout);
    this->view->SetLayoutStrategy("Pass Through");
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
	this->view->ColorVerticesOn();
	this->view->SetVertexColorArrayName("Color1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	this->theme->SetPointLookupTable(vetexlut);

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for(int i=0; i<this->num_features;i++)
    {
		pts->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1], 0.5, 0);
	}

	bool bplus = true;
	if(this->num_features == table->GetNumberOfColumns())
	{
		bplus = false;
	}
	
	for(int i=0; i<this->num_features;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(45.0);
		if(bplus)
		{
			label->InsertNextValue(this->table->GetColumn(i+1)->GetName());
		}
		else
		{
			label->InsertNextValue(this->table->GetColumn(i)->GetName());
		}
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);

	this->selectionCallback1->SetCallback (SelectionCallbackFunctionForSPD);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void Heatmap::drawPointsForOrderHeatmap()
{
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for(int i=0; i<this->num_features;i++)
    {
		pts->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1], 0.5, 0);
	}

	for(int i=0; i<this->num_features;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(45.0);
		vtkIdType id = i+1  ;
		label->InsertNextValue(this->table->GetColumn(id)->GetName());
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);

    //this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    //this->selectionCallback1->SetClientData(this);
    //this->selectionCallback1->SetCallback ( SelectionCallbackFunctionForSPD);
    //this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void Heatmap::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	Heatmap* heatmapWin = (Heatmap*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;
	vtkSelectionNode* cells = NULL;

    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(0);
        }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(0);
        }
 
    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(1);
        }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(1);
        }

	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
	
		std::set<long int> IDs;
		if(vertexList->GetNumberOfTuples() > 0)
		{

			for( vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
			{
				long int value = vertexList->GetValue(i);
				IDs.insert(value);
			}
		}
		try
		{
			heatmapWin->SetdenSelectedIds1( IDs, false);
		}
		catch(...)
		{
			cout<<"SetdenSelectedIds1 failed, please try again !"<<endl;
		}
	}
}

void Heatmap::SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	Heatmap* heatmapWin = (Heatmap*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;
	vtkSelectionNode* cells = NULL;

    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(0);
        }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(0);
        }
 
    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(1);
        }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(1);
        }

	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
	
		std::set<long int> IDs;
		if(vertexList->GetNumberOfTuples() > 0)
		{

			for( vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
			{
				long int value = vertexList->GetValue(i);
				IDs.insert(value);
			}
		}
		try
		{
			heatmapWin->SetdenSelectedIds1( IDs, true);
		}
		catch(...)
		{
			cout<<"SetdenSelectedIds1 failed, please try again !"<<endl;
		}
	}
}

void Heatmap::SelectionCallbackFunction2(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	Heatmap* heatmapWin = (Heatmap*)clientData;
	int* pos = heatmapWin->view->GetInteractor()->GetEventPosition();
 
	vtkCellPicker *cell_picker = (vtkCellPicker *)heatmapWin->view->GetInteractor()->GetPicker();

	cell_picker->Pick(pos[0], pos[1], 0, heatmapWin->view->GetRenderer());
	double* worldPosition = cell_picker->GetPickPosition();
 
	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, heatmapWin->view->GetRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();

		if(cellpicker->GetCellId() != -1)
			heatmapWin->id1 = cellpicker->GetCellId();
	}
	if(worldPosition[0]<-0.5)
	{
		std::cout<<"world position: "<<worldPosition[0] <<endl;
		std::cout<<"rescaled value: "<< - ( worldPosition[0] + 0.5) <<endl;
		heatmapWin->addDragLineforSPD(worldPosition);
	}	
}

void Heatmap::SelectionCallbackFunction3(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	Heatmap* heatmapWin = (Heatmap*)clientData;
	int* pos = heatmapWin->view->GetInteractor()->GetEventPosition();
 
	vtkCellPicker *cell_picker = (vtkCellPicker *)heatmapWin->view->GetInteractor()->GetPicker();
 
	// Pick from this location.
	cell_picker->Pick(pos[0], pos[1], 0, heatmapWin->view->GetRenderer());
	double* worldPosition = cell_picker->GetPickPosition();

	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, heatmapWin->view->GetRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();
		if(cellpicker->GetCellId() != -1)
		{
			try
			{
			heatmapWin->id2 = cellpicker->GetCellId();
			heatmapWin->ids = vtkSmartPointer<vtkIdTypeArray>::New();
			heatmapWin->ids->SetNumberOfComponents(1);
			heatmapWin->computeselectedcells();
			heatmapWin->setselectedCellIds();
			emit heatmapWin->SelChanged();
			}
			catch(...)
			{
				cout<<"SelectionCallbackFunction3 failed, please try again !"<<endl;
			}
		}
	}
}

void Heatmap::HandleKeyPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	/*Heatmap* heatmapWin = (Heatmap*)clientData;
	char key = heatmapWin->view->GetInteractor()->GetKeyCode();
	int size = 0;
	switch (key)
	{
	case 'c':
		heatmapWin->continueselect = true;
		break;
	case 'i':
		heatmapWin->intersectionselect = true;
		break;
	case 'r':
		heatmapWin->continueselect = false;
		heatmapWin->intersectionselect = false;
		break;
	case 'd':
		size = heatmapWin->Selection->getSelections().size();
		if( size > 0)
		{
			std::cout<< size<< " items deleted! Please rerun SPD analysis!"<<endl;
			heatmapWin->Selection->DeleteCurrentSelectionInTable();
		}
		else
		{
			std::cout<< "No items have been selected!"<<endl;
		}
		break;
	default:
		break;
	}*/
}

void Heatmap::SetdenSelectedIds1(std::set<long int>& IDs, bool bfirst)
{
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;
	std::set<long int>::iterator it;
	long int id;

	if(continueselect == false)
	{
		if(this->denResetflag1 == 1)
		{
			for(int i = 0; i<this->dencolors1->GetSize() ; i++)
				this->dencolors1->SetValue(i, 0);
			denlinesPolyData1->Modified();
			denlinesPolyData1->Update();
			denmapper1->Modified();
			denmapper1->Update();
			denactor1->Modified();
		}
		if(this->denResetflag2 == 1)
		{
			for(int i = 0; i<this->dencolors2->GetSize() ; i++)
			this->dencolors2->SetValue(i, 0);
			denlinesPolyData2->Modified();
			denlinesPolyData2->Update();
			denmapper2->Modified();
			denmapper2->Update();
			denactor2->Modified();
		}
	}

	if( bfirst)
	{
		if( IDs.size() > 0)
		{
			for( it = IDs.begin(); it != IDs.end(); it++ )
			{
				id = *it;
				//cout<<"id ====="<<id<<endl;
				if(id < 2*this->num_samples-1)
				{
					reselectIds1(selectedIDs1, id);
				}
				else
				{
					try
					{
						reselectIds2(selectedIDs2, id - (2*this->num_samples-1));
					}
					catch(...)
					{
						cout<<"reselectIds2 failed, please try again!"<<endl;
					}
				}
			}
		}
	}
	else
	{
		if( IDs.size() > 0)
		{
			for( it = IDs.begin(); it != IDs.end(); it++ )
			{
				id = *it;
				if(id < 2 * this->num_features - 1)
				{
					reselectIds2(selectedIDs2, id);
				}
			}
		}
	}

	if(selectedIDs1.size() > 0)	
	{
		for(int i = 0; i<this->num_features; i++)
		{
			selectedIDs2.insert(i);
		}

		denmapper1->ScalarVisibilityOn();
		denlinesPolyData1->Modified();
		denlinesPolyData1->Update();
		denmapper1->Modified();
		denmapper1->Update();
		denactor1->Modified();
		this->view->Render();
		this->denResetflag1 = 1;
		
		if(intersectionselect == true)
			this->interselectedIDs = selectedIDs1;
		try
		{
			//cout<<"set select============="<<endl;
			this->Selection2->select(selectedIDs2);
			this->Selection->select(selectedIDs1);
		}
		catch(...)
		{
			cout<<"Setcelected id after reselect failed, please try again!"<<endl;
		}
	}
	else if(selectedIDs2.size() > 0)	
	{
		this->Selection2->select(selectedIDs2);

		if(intersectionselect == false)
		{
			for(int i = 0; i<this->num_samples; i++)
				selectedIDs1.insert( indMapFromIndToVertex[i]);
			//cout<<"set select============="<<endl;
			this->Selection->select(selectedIDs1);
		}
		else
		{
			//cout<<"set select============="<<endl;
			this->Selection->select(interselectedIDs);
		}

		denmapper2->ScalarVisibilityOn();
		denlinesPolyData2->Modified();
		denlinesPolyData2->Update();
		denmapper2->Modified();
		denmapper2->Update();
		denactor2->Modified();
		this->view->Render();
		this->denResetflag2 = 1;
	}
	else
	{
		this->Selection2->clear();
		//cout<<"set select============="<<endl;
		this->Selection->clear();
	}
}

void Heatmap::reselectIds1(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_samples)
	{
		selectedIDs.insert( indMapFromIndToVertex[id]);
	}
	else
	{
		this->dencolors1->SetValue((id - this->num_samples)*9, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 1, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 2, 204);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 3, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 4, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 5, 204);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 6, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 7, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 8, 204);
		this->reselectIds1(selectedIDs, connect_Data_Tree1[rowMapForTreeData.find(id)->second][0]);
		this->reselectIds1(selectedIDs, connect_Data_Tree1[rowMapForTreeData.find(id)->second][1]);
	}
}

void Heatmap::reselectSPDIds1(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_samples)
	{
		for( int i = 0; i < indSPDMapFromIndToVertex[id].size(); i++)
		{
			selectedIDs.insert( indSPDMapFromIndToVertex[id][i]);
		}
	}
	else
	{
		this->dencolors1->SetValue((id - this->num_samples)*9, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 1, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 2, 204);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 3, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 4, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 5, 204);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 6, 153);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 7, 50);
		this->dencolors1->SetValue((id - this->num_samples)*9 + 8, 204);
		this->reselectSPDIds1(selectedIDs, connect_Data_Tree1[rowMapForTreeData.find(id)->second][0]);
		this->reselectSPDIds1(selectedIDs, connect_Data_Tree1[rowMapForTreeData.find(id)->second][1]);
	}
}

void Heatmap::reselectIds2(std::set<long int>& selectedIDs2, long int id)
{
	if(id  < this->num_features)
	{
		selectedIDs2.insert( id);                                                 
	}
	else
	{
		this->dencolors2->SetValue((id - this->num_features)*9, 153);
		this->dencolors2->SetValue((id - this->num_features)*9 + 1, 50);
		this->dencolors2->SetValue((id - this->num_features)*9 + 2, 204);
		this->dencolors2->SetValue((id - this->num_features)*9 + 3, 153);
		this->dencolors2->SetValue((id - this->num_features)*9 + 4, 50);
		this->dencolors2->SetValue((id - this->num_features)*9 + 5, 204);
		this->dencolors2->SetValue((id - this->num_features)*9 + 6, 153);
		this->dencolors2->SetValue((id - this->num_features)*9 + 7, 50);
		this->dencolors2->SetValue((id - this->num_features)*9 + 8, 204);
		this->reselectIds2(selectedIDs2, connect_Data_Tree2[columnMapForTreeData.find(id)->second][0]);
		this->reselectIds2(selectedIDs2, connect_Data_Tree2[columnMapForTreeData.find(id)->second][1]);
	}
}

void Heatmap::computeselectedcells()
{
	this->r1 = id1/this->num_features;
	this->r2 = id2/this->num_features;
	this->c1 = id1%this->num_features;
	this->c2 = id2%this->num_features;

	cout<<r1<<endl;
	cout<<r2<<endl;
	for(int i = 0; i <= r1 - r2; i++)
	{
		for(int j = 0; j <= c2 - c1; j++)
		{
			ids->InsertNextValue(id2 - j + this->num_features*i);
		}
	}
}

void Heatmap::setselectedCellIds()
{
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;

	if(continueselect == false)
	{
		if(this->denResetflag1 == 1)
		{
			for(int i = 0; i<this->dencolors1->GetSize() ; i++)
				this->dencolors1->SetValue(i, 0);
			denmapper1->ScalarVisibilityOn();
			denlinesPolyData1->Modified();
			denlinesPolyData1->Update();
			denmapper1->Modified();
			denmapper1->Update();
			denactor1->Modified();
			this->view->Render();
			this->denResetflag1 = 0;
		}

		if(this->denResetflag2 == 1)
		{
			for(int i = 0; i<this->dencolors2->GetSize() ; i++)
				this->dencolors2->SetValue(i, 0);

			denmapper2->ScalarVisibilityOn();
			denlinesPolyData2->Modified();
			denlinesPolyData2->Update();
			denmapper2->Modified();
			denmapper2->Update();
			denactor2->Modified();
			this->view->Render();
			this->denResetflag2 = 0;
		}
	}

	for(int i = r2; i<=r1; i++)
	{
		selectedIDs1.insert( indMapFromIndToVertex[ this->Optimal_Leaf_Order1[i]]);
	}
	for(int j = c1; j<=c2; j++)
	{		
		selectedIDs2.insert(this->Optimal_Leaf_Order2[j]);
	}
	try
	{
	this->Selection2->select(selectedIDs2);
	this->Selection->select(selectedIDs1);
	}
	catch(...)
	{
		cout<<"setselectedCellIds failed, please try again!"<<endl;
	}
}

void Heatmap::GetSelRowCol(int &r1, int &c1, int &r2, int &c2)
{
	r1 = this->r1;
	r2 = this->r2;
	c1 = this->c1;
	c2 = this->c2;
}

void Heatmap::SetSelRowCol(int r1, int c1, int r2, int c2)
{
	this->r1 = r1;
	this->r2 = r2;
	this->c1 = c1;
	this->c2 = c2;
	
	setselectedCellIds();	
}

void Heatmap::SetSPDInteractStyle()
{
	this->theme->SetCellValueRange(0, this->num_samples*this->num_features - 1);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
	this->view->ApplyViewTheme(theme);

	this->myCellPicker = vtkSmartPointer<vtkCellPicker>::New();
	this->view->GetInteractor()->SetPicker(this->myCellPicker);
	this->myCellPicker->SetTolerance(0.004);
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2 =vtkSmartPointer<vtkCallbackCommand>::New();
	selectionCallback2->SetClientData(this);
	selectionCallback2->SetCallback(SelectionCallbackFunction2);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->view->GetInteractor()->AddObserver(vtkCommand::RightButtonReleaseEvent, selectionCallback2);

	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
}

void Heatmap::showGraphforSPD( int selCol, int unselCol, bool bprogressionHeatmap)
{	
	if(bprogressionHeatmap)
	{
		drawPointsForOrderHeatmap();
	}
	else
	{
		drawPointsForSPD();
	}

	this->dragLineFlag = false;

	this->aPlane = vtkSmartPointer<vtkPlaneSource>::New();
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	this->cellData = vtkSmartPointer<vtkFloatArray>::New();

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	this->celllut = vtkSmartPointer<vtkLookupTable>::New();
	this->celllut->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->celllut->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->celllut->Build();

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			rgb rgb = GetRGBValue( mapdata[num_samples - i - 1][j]);
			celllut->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->aPlane->Update();
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(celllut);

	this->actor = vtkSmartPointer<vtkActor>::New();
	this->actor->SetMapper(mapper);

	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetTableRange (-1, 1);
	scalarbarLut->SetNumberOfTableValues(COLOR_MAP_SIZE);
	for(int index = 0; index<COLOR_MAP_SIZE;index++)
	{
		rgb rgbscalar = COLORMAP[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(10);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->SetMaximumHeightInPixels(1000);
	scalarBar->SetMaximumWidthInPixels(100);

	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->AddActor2D(scalarBar);
	this->SetSPDInteractStyle();
	this->view->GetRenderer()->GradientBackgroundOff();
	this->view->GetRenderer()->SetBackground(1,1,1);

	if( selCol > 0 || unselCol > 0)
	{
		double p1[3];
		double p2[3];
		double p3[3];
		double length = 0.3;
		double angle = pi / 4;


		p1[0]= -0.5 + 1.0 / (selCol + unselCol) * selCol;
		p1[1]= -0.5;
		p1[2]=0;
		p2[0]= p1[0];
		p2[1]=0.5;
		p2[2]=0;

		p3[0] = p2[0] + length * cos( angle);
		p3[1] = p2[1] + length * sin( angle);
		p3[2] = 0;

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkLineSource> lineSourceForText = vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkPolyDataMapper> lineMapperForText = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
		vtkSmartPointer<vtkActor> lineActorForText = vtkSmartPointer<vtkActor>::New();

		lineSource->SetPoint1(p1);
		lineSource->SetPoint2(p2);
		lineMapper->SetInputConnection(lineSource->GetOutputPort());
		lineActor->SetMapper(lineMapper);
		lineActor->GetProperty()->SetColor(1,1,1);
		lineActor->GetProperty()->SetLineWidth(2);

		lineSourceForText->SetPoint1(p2);
		lineSourceForText->SetPoint2(p3);
		lineMapperForText->SetInputConnection(lineSourceForText->GetOutputPort());
		lineActorForText->SetMapper(lineMapperForText);
		lineActorForText->GetProperty()->SetColor(1,0,0);
		lineActorForText->GetProperty()->SetLineWidth(2);

		this->view->GetRenderer()->AddActor(lineActor);
		this->view->GetRenderer()->AddActor(lineActorForText);
	}

	try
	{
		if( bprogressionHeatmap == false)
		{
			if(this->clusflag == true)
				this->showDendrogram1();
			else
			{
				this->showDendrogram1();
				this->showDendrogram2();
			}
		}
	}
	catch(...)
	{
		cout<<"Draw dendrogram failed ! Please try again!"<<endl;
	}
	this->view->Render();
	this->view->GetInteractor()->Start();
}

void Heatmap::addDragLineforSPD(double* worldPosition)
{
	if(this->dragLineFlag)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor);

	double p1[3];
	double p2[3];
	p1[0]=worldPosition[0];
	p1[1]=-0.75;
	p1[2]=0;
	p2[0]=worldPosition[0];
	p2[1]=0.75;
	p2[2]=0;
	dragLineSource = vtkSmartPointer<vtkLineSource>::New();
	dragLineSource->SetPoint1(p1);
	dragLineSource->SetPoint2(p2);

	dragLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	dragLineMapper->SetInputConnection(dragLineSource->GetOutputPort());

	dragLineActor = vtkSmartPointer<vtkActor>::New();
	dragLineActor->SetMapper(dragLineMapper);
	dragLineActor->GetProperty()->SetColor(0.5,0.7,0);
	dragLineActor->GetProperty()->SetLineWidth(1.5);
	dragLineActor->DragableOn();
	this->view->GetRenderer()->AddActor(dragLineActor);
	this->dragLineFlag = true;
	this->view->Render();

	this->selectClustersforSPD(worldPosition);
}

void Heatmap::selectClustersforSPD(double* worldPosition)
{
	reselectedClusterSPD.clear();
	clusterNumVec.set_size(parentIndex.size());
	for( int i = 0; i < clusterNumVec.size(); i++)
	{
		clusterNumVec[i] = 0;
	}
	std::set<long int> selectedClusterSPD;
	for( long int index = 0; index < 2 * this->num_samples - 1; index++)
	{
		if(worldPosition[0] > Processed_Coordinate_Data_Tree1[index][1])
		{
			selectedClusterSPD.insert(index);
		}
	}

	std::set<long int>::iterator it;
	for(it = selectedClusterSPD.begin(); it != selectedClusterSPD.end(); it++)
	{
		long int id = *it;
		long int id1 = connect_Data_Tree1[rowMapForTreeData.find(id)->second][0];
		long int id2 = connect_Data_Tree1[rowMapForTreeData.find(id)->second][1];

		std::set<long int>::iterator it1 = selectedClusterSPD.find(id1);
		std::set<long int>::iterator it2 = selectedClusterSPD.find(id2);

		unsigned char num = 0;
		if( it1 == selectedClusterSPD.end() && it2 != selectedClusterSPD.end() || it2 == selectedClusterSPD.end() && it1 != selectedClusterSPD.end())
		{
			num = 1;
		}
		else if( it1 == selectedClusterSPD.end() && it2 == selectedClusterSPD.end())
		{
			num = 2;
		}
		else
		{
			num = 0;
		}

		if( id <= parentIndex[0])
		{
			clusterNumVec[0] += num;
		}
		else
		{
			for( int j = 0; j < parentIndex.size() - 1; j++)
			{
				if( id > parentIndex[j] && id <= parentIndex[j+1])
				{
					
					clusterNumVec[j+1] += num;
					break;
				}
			}
		}
	}

	for( int i = 0; i < clusterNumVec.size(); i++)
	{
		if( clusterNumVec[i] == 0)
		{
			 clusterNumVec[i] = 1;
		}
	}
	
	this->reselectClustersforSPD(selectedClusterSPD);
}

void Heatmap::GetSubTreeClusterNum(std::vector<int> &clusterNum)
{
	clusterNum.resize( clusterNumVec.size());
	for( int i = 0; i < clusterNumVec.size(); i++)
	{
		clusterNum[i] = clusterNumVec[i];
	}
}

void Heatmap::reselectClustersforSPD(std::set<long int>& selectedClusterSPD)
{
	std::set<long int>::iterator it;
	int clusternumber = 0;
	for(it = selectedClusterSPD.begin(); it != selectedClusterSPD.end(); it++)
	{
		long int id = *it;
		long int id1 = connect_Data_Tree1[rowMapForTreeData.find(id)->second][0];
		long int id2 = connect_Data_Tree1[rowMapForTreeData.find(id)->second][1];

		std::set<long int>::iterator it1 = selectedClusterSPD.find(id1);
		std::set<long int>::iterator it2 = selectedClusterSPD.find(id2);

		if(it1 == selectedClusterSPD.end())
		{
			reselectedClusterSPD.insert(id1);
			clusternumber++;
		}
		if(it2 == selectedClusterSPD.end())
		{
			reselectedClusterSPD.insert(id2);
			clusternumber++;
		}	
	}
	//cout<<"cluster number is "<<clusternumber<<endl;

	std::vector< std::vector< long int> > clusIndex;
	std::vector< std::vector< long int> > sampleIndex;

	for(it = reselectedClusterSPD.begin(); it != reselectedClusterSPD.end(); it++)
	{
		std::vector<long int> sampleIdsVec;
		std::set<long int> clusIdsforSPD;
		std::vector<long int> clusIdsVec;
		this->reselectIdsforSPD(*it, &clusIdsforSPD);

		for( int i = 0; i < this->num_samples; i++)
		{
			if(clusIdsforSPD.find(Optimal_Leaf_Order1[i]) != clusIdsforSPD.end())
			{
				for( int k = 0; k < clusIdsforSPD.size(); k++)
				{
					int ind = Optimal_Leaf_Order1[i + k];

					clusIdsVec.push_back(ind);
					for( int j = 0; j < indSPDMapFromIndToVertex[ind].size(); j++)
					{
						sampleIdsVec.push_back( indSPDMapFromIndToVertex[ind][j]);
					}
				}
				break;
			}
		}
		sampleIndex.push_back(sampleIdsVec);
		clusIndex.push_back(clusIdsVec);
	}
	this->Selection->SetClusterIndex(clusIndex);
	this->Selection->SetSampleIndex(sampleIndex);
}

vtkSmartPointer<vtkTable> Heatmap::GetTreeTable()
{
	std::vector<std::string> headers;
	headers.push_back("node1");
	headers.push_back("node2");
	headers.push_back("weight");

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	for(int i = 0; i < headers.size(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( headers[i].c_str());
		table->AddColumn(column);
	}

	/// ordering the tree by weight
	std::multimap< double, int> treeMap;
	for( int i = 0; i < this->table->GetNumberOfRows() - 1; i++)
	{
		treeMap.insert(std::pair<double, int>(ftreedata[i][2], i));
	}

	std::vector<int> orderInd(this->table->GetNumberOfRows() - 1); 
	std::multimap< double, int>::iterator treeMapIter;
	int mt = 0;
	for(treeMapIter = treeMap.begin(); treeMapIter != treeMap.end(); treeMapIter++)
	{
		orderInd[mt++] = treeMapIter->second;
	}

	if( reselectedClusterSPD.size() > 0 && this->ftreedata)
	{
		std::set<long int>::iterator it;
		std::map<long int, int> idMap;
		int ind = 0;
		for(it = reselectedClusterSPD.begin(); it != reselectedClusterSPD.end(); it++)
		{
			idMap.insert(std::pair<long int, int>(*it, ind++));
		}

		for( int i = this->table->GetNumberOfRows() - reselectedClusterSPD.size(); i < this->table->GetNumberOfRows() - 1; i++)
		{
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
			int pos = orderInd[i];
			int node1 = ftreedata[pos][0];
			int node2 = ftreedata[pos][1];
			int parent = ftreedata[pos][3];
			DataRow->InsertNextValue(idMap[ node1]);
			DataRow->InsertNextValue(idMap[ node2]);
			DataRow->InsertNextValue(pow(this->ftreedata[i][2], POWER_PARAM));
			table->InsertNextRow(DataRow);
			idMap.insert(std::pair<long int, int>(parent, idMap[ node2]));
		}
	}
	//ftk::SaveTable("TreeTable.txt", table);
	return table;
}

void Heatmap::reselectIdsforSPD(long int id, std::set<long int> *clusidforSPD)
{
	if(id < this->num_samples)
	{
		//cout<<id<<"\t";
		if( clusidforSPD != NULL)
		{
			clusidforSPD->insert(id);
		}
	}
	else
	{
		this->reselectIdsforSPD(connect_Data_Tree1[rowMapForTreeData.find(id)->second][0], clusidforSPD);
		this->reselectIdsforSPD(connect_Data_Tree1[rowMapForTreeData.find(id)->second][1], clusidforSPD);
	}
}

void Heatmap::setModelsforSPD(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, std::vector< int> selOrder, std::vector< int> unselOrder, std::map< int, int> *indexCluster,
							   std::vector<int> *component, int numOfComponenets, vnl_matrix<double> *subTreeDistance, ObjectSelection * sels2)
{
	this->table = table;

	this->indMapFromVertexToInd.clear();
	this->indMapFromIndToVertex.clear();
	std::cout<< this->table->GetNumberOfRows()<<"\t"<<this->table->GetNumberOfColumns()<<std::endl;
	if( indexCluster) 
	{
		indMapFromVertexToInd = *indexCluster;
		indSPDMapFromIndToVertex.resize( table->GetNumberOfRows());
		std::map<int, int>::iterator iter;
		for( iter = indMapFromVertexToInd.begin(); iter != indMapFromVertexToInd.end(); iter++)
		{
			std::pair<int, int> pair = *iter;
			indSPDMapFromIndToVertex[ pair.second].push_back( pair.first);
		}
	}

	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;

	if(!sels2)
		this->Selection2 = new ObjectSelection();
	else
		this->Selection2 = sels2;

	if(component)
	{
		connectedComponent = *component;
	}

	selectedFeatureIDs.clear();
	for( int i = 0; i < selOrder.size(); i++)
	{
		selectedFeatureIDs.insert( selOrder[i]);
	}
	std::cout<< selectedFeatureIDs.size()<<std::endl;
	//connect(Selection, SIGNAL(thresChanged()), this, SLOT(GetSelecectedIDsforSPD()));
	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDsForSPD()));

	if( subTreeDistance != NULL && numOfComponenets >= 1)
	{
		this->runClusforSPD(selOrder, unselOrder, numOfComponenets, *subTreeDistance);
	}
}

void Heatmap::setModelsforSPD(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, std::vector< int> sampleOrder, std::vector< int> selOrder, std::vector< int> unselOrder, std::map< int, int> *indexCluster, ObjectSelection * sels2)
{
	this->table = table;
	this->indMapFromVertexToInd.clear();
	this->indMapFromIndToVertex.clear();

	if( indexCluster) 
	{
		indMapFromVertexToInd = *indexCluster;
		indSPDMapFromIndToVertex.resize( table->GetNumberOfRows());
		std::map<int, int>::iterator iter;
		for( iter = indMapFromVertexToInd.begin(); iter != indMapFromVertexToInd.end(); iter++)
		{
			std::pair<int, int> pair = *iter;
			indSPDMapFromIndToVertex[ pair.second].push_back( pair.first);
		}
	}

	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;

	if(!sels2)
		this->Selection2 = new ObjectSelection();
	else
		this->Selection2 = sels2;

	selectedFeatureIDs.clear();
	for( int i = 0; i < selOrder.size(); i++)
	{
		selectedFeatureIDs.insert( selOrder[i]);
	}
	//connect(Selection, SIGNAL(thresChanged()), this, SLOT(GetSelecectedIDsforSPD()));
	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDsForSPD()));
	this->runClusforSPD( sampleOrder, selOrder, unselOrder);
}

void Heatmap::runClusforSPD(std::vector< int> selOrder, std::vector< int> unselOrder, int numberofcomponents, vnl_matrix<double> &subTreeDistance)
{
	this->clusflag = true;

	if(connectedComponent.size() == this->table->GetNumberOfRows())
	{
		clock_t start_time = clock();
		std::vector< std::vector<int> > graphVertex;
		graphVertex.resize(numberofcomponents);

		for( int i = 0; i < connectedComponent.size(); i++)
		{
			int n = connectedComponent[i];
			graphVertex[n].push_back(i);
		}

		double** datas = new double*[this->table->GetNumberOfRows()];
		
		for (int i = 0; i < this->table->GetNumberOfRows(); i++)
		{
			datas[i] = new double[this->table->GetNumberOfColumns() - 1 + 2 ];
			for(int j = 1; j < this->table->GetNumberOfColumns(); j++)
			{
				vtkVariant temp = this->table->GetValue(i, j);
				datas[i][j-1] = temp.ToDouble();
			}
		}

		std::vector< std::vector< ClusterTree> > treeVec;
		treeVec.resize( numberofcomponents);

		#pragma omp parallel for
		for( int i = 0; i < graphVertex.size(); i++)
		{
			double** datasforclus = new double*[graphVertex[i].size()];
			for (int j = 0; j < graphVertex[i].size(); j++)
			{
				datasforclus[j] = new double[selOrder.size() + 2 ];
			}

			int index = 0;
			for(std::set<long int >::iterator iter = selectedFeatureIDs.begin(); iter != selectedFeatureIDs.end(); iter++)
			{
				for( int k = 0; k < graphVertex[i].size(); k++)
				{
					int ni = graphVertex[i][k];
					datasforclus[k][index] = datas[ni][*iter];
				}
				index++;
			}

			clusclus *cc = new clusclus(datasforclus, (int)graphVertex[i].size(), (int)selOrder.size());
			cc->RunClusClus();
			
			///save the tree 
			cc->GetTreeStructure(treeVec[i]);
			
			for (int j = 0; j < graphVertex[i].size(); j++)
			{
				delete datasforclus[j];
			}
			delete datasforclus;
			delete cc;
		}

		//std::ofstream ofs("Tree.txt");
		/// adjust the local tree index to fit the global tree index
	 	int clusNo = this->table->GetNumberOfRows();
		double maxdis = 0; // find the maximum distance of all the tree
		std::vector< int> parentIndex;
		parentIndex.resize(treeVec.size());
		for( int i = 0; i < treeVec.size(); i++)
		{
			int maxi = treeVec[i].size() + 1;
			//ofs<<maxi<<std::endl;
			parentIndex[i] = 0;
			for( int j = 0; j < treeVec[i].size(); j++)
			{
				ClusterTree tree = treeVec[i][j];
				tree.first = tree.first < maxi ? graphVertex[i][tree.first] : tree.first - maxi + clusNo;
				tree.second = tree.second < maxi ? graphVertex[i][tree.second] : tree.second - maxi + clusNo;
				tree.parent = tree.parent < maxi ? graphVertex[i][tree.parent] : tree.parent - maxi + clusNo;
				treeVec[i][j] = tree;
				if(tree.dis > maxdis)
				{
					maxdis = tree.dis;
				}
			}
			clusNo += treeVec[i].size();
			parentIndex[i] = clusNo - 1;
		}
		this->parentIndex = parentIndex;
		/// print out the tree structure
		//for( int i = 0; i < treeVec.size(); i++)
		//{
		//	for( int j = 0; j < treeVec[i].size(); j++)
		//	{
		//		ClusterTree tree = treeVec[i][j];
		//		ofs<< tree.first<< "\t"<< tree.second << "\t"<<tree.dis<< "\t"<< tree.parent <<std::endl;
		//	}
		//	ofs<< std::endl;
		//}

		/// construct the overal tree and the overal order
		std::vector< ClusterTree> overallTree;
		while(overallTree.size() < numberofcomponents - 1)  
		{
			unsigned int location = subTreeDistance.arg_min();
			int nrow = location / subTreeDistance.rows();
			int ncol = location % subTreeDistance.rows();
			if( parentIndex[nrow] != parentIndex[ncol])
			{
				ClusterTree tree( parentIndex[nrow], parentIndex[ncol], 2 * maxdis, clusNo);
				overallTree.push_back(tree);
				subTreeDistance(ncol, nrow)= 1e9;
				subTreeDistance(nrow, ncol)= 1e9;
				int ind1 = parentIndex[nrow];
				int ind2 = parentIndex[ncol];

				for( int i = 0; i < parentIndex.size(); i++)
				{
					if( parentIndex[i] == ind1 || parentIndex[i] == ind2)
					{
						parentIndex[i] = clusNo;
					}
				}
				clusNo++;
			}
			else
			{
				subTreeDistance(ncol, nrow)= 1e9;
				subTreeDistance(nrow, ncol)= 1e9;
			}
		}

		//for( int i = 0; i < overallTree.size(); i++)
		//{
		//	ClusterTree tree = overallTree[i];
		//	//ofs<< tree.first<< "\t"<< tree.second << "\t"<<tree.dis<< "\t"<< tree.parent <<std::endl;
		//}
		//ofs.close();

		int featureNum = selOrder.size() + unselOrder.size();
		int* optimalleaforder2 = new int[featureNum];
		int counter = 0;
		for(int i = 0; i < selOrder.size(); i++)
		{
			optimalleaforder2[i] = selOrder[i];
			counter++;
		}
		for(int i = 0; i < unselOrder.size(); i++)
		{
			optimalleaforder2[i + counter] = unselOrder[i];
		}

		if(ftreedata)
		{
			for( int i = 0; i < this->table->GetNumberOfRows() - 1; i++)
			{
				delete ftreedata[i];
			}
			delete ftreedata;
		}
		ftreedata = new double*[ this->table->GetNumberOfRows() - 1];
		int count = 0;
		for( int i = 0; i < treeVec.size(); i++)
		{
			for( int j = 0; j < treeVec[i].size(); j++)
			{
				ClusterTree tree = treeVec[i][j];
				ftreedata[count] = new double[4];
				ftreedata[count][0] = tree.first;
				ftreedata[count][1] = tree.second;
				ftreedata[count][2] = tree.dis;
				ftreedata[count][3] = tree.parent;
				count++;
			}
		}

		for(int i = 0; i < overallTree.size(); i++)
		{	
			ClusterTree tree = overallTree[i];
			ftreedata[count] = new double[4];
			ftreedata[count][0] = tree.first;
			ftreedata[count][1] = tree.second;
			ftreedata[count][2] = tree.dis;
			ftreedata[count][3] = tree.parent;
			count++;
		}

		clusclus *cc1 = new clusclus();
		cc1->Initialize(ftreedata, this->table->GetNumberOfRows());
		cc1->GetOptimalLeafOrderD();

		this->setDataForHeatmap( datas, cc1->optimalleaforder, optimalleaforder2, this->table->GetNumberOfRows(), featureNum);
		this->setDataForDendrograms(cc1->treedata);
		this->creatDataForHeatmap(POWER_PARAM);
		
		for (int i = 0; i < this->table->GetNumberOfRows(); i++)
		{
			delete datas[i];
		}
		delete datas;
		std::cout << "Total time to generate progression heatmap is: " << (clock() - start_time) / (float) CLOCKS_PER_SEC << std::endl;
	}
	else
	{
		std::cout<< "Error in component array size"<<std::endl;
	}
}

void Heatmap::runClusforSPD(std::vector< int> sampleOrder, std::vector< int> selOrder, std::vector< int> unselOrder)
{
	this->clusflag = true;
	double** datas;
	vtkVariant temp; 
	int nFeature = selOrder.size() + unselOrder.size();

	datas = new double*[this->table->GetNumberOfRows()];

	std::cout<<this->table->GetNumberOfRows()<<endl;
	std::cout<<this->table->GetNumberOfColumns()<<endl;

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		datas[i] = new double[nFeature + 2 ];
		for(int j = 1; j <= nFeature; j++)
		{
			vtkVariant temp = this->table->GetValue(i, j);
			datas[i][j-1] = temp.ToDouble();
		}
	}

	int* optimalleaforder1 = new int[this->table->GetNumberOfRows()];
	for(int i = 0; i < sampleOrder.size(); i++)
	{
		optimalleaforder1[i] = sampleOrder[i];
	}
	
	int* optimalleaforder2 = new int[ nFeature];
	int counter = 0;
	for(int i = 0; i < selOrder.size(); i++)
	{
		optimalleaforder2[i] = selOrder[i];
		counter++;
	}
	for(int i = 0; i < unselOrder.size(); i++)
	{
		optimalleaforder2[i + counter] = unselOrder[i];
	}

	this->setDataForHeatmap( datas, optimalleaforder1, optimalleaforder2, this->table->GetNumberOfRows(), nFeature);
	//this->setDataForDendrograms();
	this->creatDataForHeatmap(POWER_PARAM);
	
	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		delete datas[i];
	}
	delete datas;
}

void Heatmap::SelectionCallbackFunctionForSPD(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	Heatmap* heatmapWin = (Heatmap*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;
	vtkSelectionNode* cells = NULL;

    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
        {                                                                            
        vertices = selection->GetNode(0);
        }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(0);
        }
 
    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(1);
        }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(1);
        }

	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
	
		std::set<long int> IDs;
		if(vertexList->GetNumberOfTuples() > 0)
		{
			for( vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
			{
				long int value = vertexList->GetValue(i);
				if( value < 2 * heatmapWin->num_samples - 1)
				{
					if( heatmapWin->reselectedClusterSPD.find(value) != heatmapWin->reselectedClusterSPD.end())
					{
						IDs.insert(value);
					}
				}
				else
				{
					value = value - 2 * heatmapWin->num_samples + 1;
					emit heatmapWin->columnToColorChanged(value);
					break;
				}			
			}
		}
		try
		{
			heatmapWin->SetdenSelectedIdsForSPD( IDs);
		}
		catch(...)
		{
			cout<<"SetdenSelectedIds1 failed, please try again !"<<endl;
		}
	}
}

void Heatmap::SetdenSelectedIdsForSPD(std::set<long int>& IDs)
{
	std::set<long int> selectedIDs1;
	std::set<long int>::iterator it;
	long int id;

	if(continueselect == false)
	{
		if(this->denResetflag1 == 1)
		{
			for(int i = 0; i<this->dencolors1->GetSize() ; i++)
				this->dencolors1->SetValue(i, 0);
			denlinesPolyData1->Modified();
			denlinesPolyData1->Update();
			denmapper1->Modified();
			denmapper1->Update();
			denactor1->Modified();
		}
	}

	if( IDs.size() > 0)
	{
		for( it = IDs.begin(); it != IDs.end(); it++ )
		{
			id = *it;

			if(id < 2*this->num_samples-1)
			{
				reselectSPDIds1(selectedIDs1, id);
			}
		}
	}

	if(selectedIDs1.size() > 0)	
	{
		denmapper1->ScalarVisibilityOn();
		denlinesPolyData1->Modified();
		denlinesPolyData1->Update();
		denmapper1->Modified();
		denmapper1->Update();
		denactor1->Modified();
		this->view->Render();
		this->denResetflag1 = 1;

		//std::set<long int>::iterator iter;
		//std::set<long int> selectedIDs;
		//for( iter = selectedIDs1.begin(); iter != selectedIDs1.end(); iter++)
		//{
		//	long int var = *iter;
		//	for( int i = 0; i < indSPDMapFromIndToVertex[ var].size(); i++)
		//	{
		//		selectedIDs.insert( indSPDMapFromIndToVertex[ var][i]);
		//	}
		//}
		
		try
		{
			this->Selection2->select(selectedFeatureIDs);
			this->Selection->select(selectedIDs1);

		}
		catch(...)
		{
			cout<<"Setcelected id after reselect failed, please try again!"<<endl;
		}
	}
	//else
	//{
	//	this->Selection2->clear();
	//	this->Selection->clear();
	//}
}
void Heatmap::GetSelecectedIDsForSPD()
{
	std::set<long int> selectedIDs2 = selectedFeatureIDs;
	std::set<long int> selectedIDs1 = this->Selection->getSelections();	
	std::set<long int>::iterator iter1 = selectedIDs1.begin();
	std::set<long int>::iterator iter2 = selectedIDs2.begin();
	vtkSmartPointer<vtkIdTypeArray> cellids = vtkSmartPointer<vtkIdTypeArray>::New();
	cellids->SetNumberOfComponents(1);

	int num1 = selectedIDs1.size();
	int num2 = selectedIDs2.size();

	std::vector<int > IDs1;
	std::vector<int > IDs2;
	IDs1.resize(num1);
	IDs2.resize(num2);

	int count1 = 0;
	int count2 = 0;

	#pragma omp parallel sections
	{
		#pragma omp section
		while(iter1 != selectedIDs1.end())
		{
			int index1 = *iter1;
			int var = indMapFromVertexToInd.find(index1)->second;
			int id1 = rowMapFromOriginalToReorder.find(var)->second;
			IDs1[count1++] = id1;
			iter1++;
		}

		#pragma omp section 
		while(iter2 != selectedIDs2.end())
		{
			int index2 = *iter2;
			int id2 = columnMapFromOriginalToReorder.find(index2)->second;
			IDs2[count2++] = id2;
			iter2++;
		}
	}

	for(int i = 0; i<num1; i++)
		for(int j = 0; j<num2; j++)
			cellids->InsertNextValue( (IDs1[i])*(this->num_features) + IDs2[j]);

	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(cellids);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);
	 
	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInput(0, this->aPlane->GetOutput());
	extractSelection->SetInput(1, selection);
	extractSelection->Update();
	 
	vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());
	
	vtkSmartPointer<vtkDataSetMapper> selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selectedMapper->SetInputConnection(selected->GetProducerPort());
	 
	vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->GetProperty()->EdgeVisibilityOn();
	selectedActor->GetProperty()->SetEdgeColor(1,1,1);
	selectedActor->GetProperty()->SetLineWidth(0.5);

	if(this->dragLineFlag)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor);
	try
	{
		if(continueselect == false)
		{
			if(continueselectnum > 0)
			{
				cout<<"I'm here "<<continueselectnum<<endl;
				for(int i = 0; i<continueselectnum + 1; i++)
					this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
				continueselectnum = 0;
			}
			else
			{
				if (this->removeActorflag != 0)
					this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
			}

			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
			this->removeActorflag += 1;
		}
		else
		{
			this->continueselectnum += 1;
			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);

		}
		if(this->dragLineFlag)
			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor);
		this->view->Render();
	}
	catch(...)
	{
		cout<<"GetSelecectedIDs failed, please try it again!"<<endl;
	}
}

void Heatmap::closeEvent(QCloseEvent *event)
{
	mainQTRenderWidget.close();
}

void Heatmap::reRunClus()
{
	double** datas;
	vtkVariant temp; 

	datas = new double*[this->table->GetNumberOfRows()];

	std::cout<<this->table->GetNumberOfRows()<<endl;
	std::cout<<this->table->GetNumberOfColumns()<<endl;

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		datas[i] = new double[this->table->GetNumberOfColumns() - 1 + 2 ];
	}

	for(int i = 0; i < this->table->GetNumberOfRows(); i++)
	{		
		for(int j = 1; j < this->table->GetNumberOfColumns(); j++)
		{
			temp = this->table->GetValue(i, j);
			datas[i][j-1] = temp.ToDouble();
		}
	}
	int* optimalleaforder1 = new int[this->table->GetNumberOfRows()];
	for(int i = 0;i<this->table->GetNumberOfRows(); i++)
		optimalleaforder1[i]=i;
	int* optimalleaforder2 = new int[this->table->GetNumberOfColumns() - 1];
	for(int i = 0;i<this->table->GetNumberOfColumns() - 1; i++)
		optimalleaforder2[i]=i;
	
	this->setDataForHeatmap(datas, optimalleaforder1, optimalleaforder2,this->table->GetNumberOfRows(), this->table->GetNumberOfColumns() - 1);
	this->creatDataForHeatmap(POWER_PARAM);	
}

void Heatmap::showGraphforNe()
{	
	this->drawPointsforNe();


	this->aPlane = vtkSmartPointer<vtkPlaneSource>::New();
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	this->cellData = vtkSmartPointer<vtkFloatArray>::New();

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	this->celllut = vtkSmartPointer<vtkLookupTable>::New();
	this->celllut->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->celllut->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->celllut->Build();

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			rgb rgb = GetRGBValue( mapdata[num_samples - i - 1][j]);
			celllut->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->aPlane->Update();
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(celllut);

	this->actor = vtkSmartPointer<vtkActor>::New();
	this->actor->SetMapper(mapper);

	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetTableRange (-1, 1);
	scalarbarLut->SetNumberOfTableValues(COLOR_MAP_SIZE);
	for(int index = 0; index<COLOR_MAP_SIZE;index++)
	{
		rgb rgbscalar = COLORMAP[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(10);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->SetMaximumHeightInPixels(1000);
	scalarBar->SetMaximumWidthInPixels(100);

	
	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->AddActor2D(scalarBar);
	this->SetInteractStyle();
	this->view->GetRenderer()->GradientBackgroundOff();
	this->view->GetRenderer()->SetBackground(1,1,1);

	//this->showDendrogram2();
	this->view->Render();
	//this->view->GetInteractor()->Start();
	return;
}

void Heatmap::drawPointsforNe()
{
	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
    this->view->AddRepresentationFromInput(graph_Layout);
    this->view->SetLayoutStrategy("Pass Through");
    this->view->ScaledGlyphsOn();

	this->theme = vtkSmartPointer<vtkViewTheme>::New();

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for(int i=0; i<this->num_features;i++)
    {
		pts->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1], 0.5, 0);
	}

	for(int i=0; i<this->num_features;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(45.0);
		vtkIdType id = i+1  ;
		label->InsertNextValue(this->table->GetColumn(id)->GetName());
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);
	this->selectionCallback1->SetCallback (SelectionCallbackFunction1);

    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);

}