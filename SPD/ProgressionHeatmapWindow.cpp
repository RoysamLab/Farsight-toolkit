#include "ProgressionHeatmapWindow.h"
#include <vtkPointData.h>
#include <vtkTextProperty.h>
#include <vtkLabelPlacementMapper.h>
#include <vtkQtLabelRenderStrategy.h>
#include <vtkActor2D.h>
#include <vtkScalarBarActor.h>

ProgressionHeatmap::ProgressionHeatmap(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;
	this->aPlane = vtkSmartPointer<vtkPlaneSource>::New();
	this->cellData = vtkSmartPointer<vtkFloatArray>::New();
	this->lookuptable = vtkSmartPointer<vtkLookupTable>::New();
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->actor = vtkSmartPointer<vtkActor>::New();
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	this->cellColors = vtkSmartPointer<vtkIntArray>::New();
	this->ids1 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->ids2 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->ids1->SetNumberOfComponents(1);
	this->ids2->SetNumberOfComponents(1);	
	this->v1 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->graph_Layout1 = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points1 = vtkSmartPointer<vtkPoints>::New();
	this->vertexColors1 = vtkSmartPointer<vtkIntArray>::New();
	this->lookupTable1 = vtkSmartPointer<vtkLookupTable>::New();
	this->myCellPicker = vtkSmartPointer<vtkCellPicker>::New();

	this->denpoints1 = vtkSmartPointer<vtkPoints>::New();
	this->denlines1 = vtkSmartPointer<vtkCellArray>::New();
	this->denlinesPolyData1 =vtkSmartPointer<vtkPolyData>::New();
	this->denmapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->denactor1 = vtkSmartPointer<vtkActor>::New();
	this->dencolors1 = vtkSmartPointer<vtkUnsignedCharArray>::New();

	this->denpoints2 = vtkSmartPointer<vtkPoints>::New();
	this->denlines2 = vtkSmartPointer<vtkCellArray>::New();
	this->denlinesPolyData2 =vtkSmartPointer<vtkPolyData>::New();
	this->denmapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->denactor2 = vtkSmartPointer<vtkActor>::New();
	this->dencolors2 = vtkSmartPointer<vtkUnsignedCharArray>::New();

	this->removeActorflag = 0;
	this->denResetflag1 = 0;
	this->denResetflag2 = 0;
	this->continueselectnum = 0;
	this->continueselect = false;
	this->intersectionselect = false;

	this->mapdata = NULL;
	this->Optimal_Leaf_Order1 = NULL;
	this->Optimal_Leaf_Order2 = NULL;
	this->connect_Data_Tree1 = NULL;
	this->connect_Data_Tree2 = NULL;

	 //this->hoverWidget = vtkSmartPointer<vtkHoverWidget>::New();
	 //hoverCallback = vtkSmartPointer<vtkHoverCallback>::New();
}

ProgressionHeatmap::~ProgressionHeatmap()
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
}

void ProgressionHeatmap::setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features)
{
	this->num_samples = num_samples;
	this->num_features = num_features;

	this->mapdata = new double*[num_samples];
	for(int i=0; i<num_samples; i++)
	{
		this->mapdata[i] = new double[num_features];
		for(int j = 0 ; j<num_features; j++)
			this->mapdata[i][j] = features[i][j];
	}

	this->Optimal_Leaf_Order1 = new int[num_samples] ;
	for(int i=0; i<num_samples; i++)
		this->Optimal_Leaf_Order1[i] = optimalleaforder1[i];

	this->Optimal_Leaf_Order2 = new int[num_features];
	for(int i=0; i<num_features; i++)
		this->Optimal_Leaf_Order2[i] = optimalleaforder2[i];
}

void ProgressionHeatmap::setDataForSimilarMatrixHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features)
{
	double max = features[0][0];
	for( int i = 0; i < num_samples; i++)
	{
		if( !this->table)    // if the table is null, build the index mapping here
		{
			this->indMapFromVertexToInd.insert( std::pair< int, int>(i, i));
			this->indMapFromIndToVertex.push_back( i);
		}
		for( int j = 0; j < num_features; j++)
		{
			if( features[i][j] > max)
			{
				max = features[i][j];
			}
		}
	}

	for( int i = 0; i < num_samples; i++)
	{
		for( int j = 0; j < num_features; j++)
		{
			features[i][j] = features[i][j] / max;
		}
	}

	this->mapdata = features;
	this->Optimal_Leaf_Order1 = optimalleaforder1;
	this->Optimal_Leaf_Order2 = optimalleaforder2;
	this->num_samples = num_samples;
	this->num_features = num_features;
}

void ProgressionHeatmap::creatDataForHeatmap(double powCof)
{
	this->scaleData();

	vector<double > temp;
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
			mean += temp[j];
		}

		mean = mean / this->num_features;

		for(int j = 0; j < num_features; j++)
			sum += (temp[j] - mean)*(temp[j] - mean);

		std = sqrt(sum/this->num_features);

		if(std)
			for(int j = 0; j < num_features; j++)
				tempdata[i][j] = (temp[j] - mean)/std;
	}
	for(int i = 0; i < this->num_samples; i++)
		mapdata[this->num_samples - i - 1] = tempdata[Optimal_Leaf_Order1[i]]; 

	FILE *fp = fopen("mapdata.txt","w");
	for(int i=0; i<num_samples; i++)
	{
		for(int j=0; j<num_features; j++)
			fprintf(fp,"%f\t",mapdata[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	if( this->connect_Data_Tree1 != NULL)
	{
		this->createDataForDendogram1(powCof);
	}
	if( this->connect_Data_Tree2 != NULL)
	{
		this->createDataForDendogram2(powCof);
	}
}

void ProgressionHeatmap::creatDataForProgressionHeatmap(double powCof)
{
	vector<double > temp;
	temp.resize(num_features);
	FILE *fp = fopen("progressionmapdata.txt","w");

	double** tempdata;
	tempdata = new double*[this->num_samples];
	for(int i = 0; i < this->num_samples; i++)
		tempdata[i] = new double[this->num_features];

	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			tempdata[i][j] = mapdata[i][Optimal_Leaf_Order2[j]];
		}
	}

	for(int i=0; i<num_samples; i++)
	{
		for(int j=0; j<num_features; j++)
			fprintf(fp,"%f\t",tempdata[i][j]);
		fprintf(fp,"\n");
	}

	for(int i = 0; i < this->num_samples; i++)
		mapdata[this->num_samples - i - 1] = tempdata[Optimal_Leaf_Order1[i]]; 

	for(int i=0; i<num_samples; i++)
	{
		for(int j=0; j<num_features; j++)
			fprintf(fp,"%f\t",mapdata[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	if( this->connect_Data_Tree1 != NULL)
	{
		this->createDataForDendogram1(powCof);
	}
	if( this->connect_Data_Tree2 != NULL)
	{
		this->createDataForDendogram2(powCof);
	}
}

void ProgressionHeatmap::creatDataForSimilarMatrixHeatmap()
{
	//const char* filename = "heatmapdata.txt";
	//FILE *fp = fopen(filename,"w");
	//for(int i=0; i<this->num_samples; i++)
	//{
	//	for(int j=0; j<this->num_features; j++)
	//		fprintf(fp,"%.4f\t",mapdata[i][j]);
	//	fprintf(fp,"\n");
	//}
	//fprintf(fp,"\n");

	double** tempdata;
	tempdata = new double*[this->num_samples];
	for(int i = 0; i < this->num_samples; i++)
	{
		tempdata[i] = new double[this->num_features];
	}

	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			tempdata[i][j] = mapdata[i][Optimal_Leaf_Order2[j]];
		}
	}
	for(int i = 0; i < this->num_samples; i++)
	{
		mapdata[i] = tempdata[Optimal_Leaf_Order1[i]]; 
	}

	//for(int i=0; i<this->num_samples; i++)
	//{
	//	for(int j=0; j<this->num_features; j++)
	//		fprintf(fp,"%.4f\t",mapdata[i][j]);
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
}

void ProgressionHeatmap::scaleData()
{
	for(int i = 0; i<this->num_features; i++)
	{
		double sum = 0.0;
		for(int j = 0; j<this->num_samples; j++)
		{
			sum += mapdata[j][i]*mapdata[j][i];
		}
		sum = sqrt(sum);

		if(sum)
			for(int j = 0; j<this->num_samples; j++)
				mapdata[j][i] /= sum;
	}
}

void ProgressionHeatmap::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
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

void ProgressionHeatmap::runClusclus()
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

	cc1 = new clusclus(datas, (int)this->table->GetNumberOfRows(), (int)this->table->GetNumberOfColumns() - 1);
	cc1->RunClusClus();
	cc1->WriteClusteringOutputToFile("mergers.txt","features.txt","progress.txt", "members.txt",
		"gap.txt", "treedata.txt", "Optimalleaforder.txt");

	cc1->Transpose();
	cc2 = new clusclus(cc1->transposefeatures,cc1->num_features, cc1->num_samples);
	cc2->RunClusClus();
	cc2->WriteClusteringOutputToFile("mergers2.txt","features2.txt","progress2.txt", "members2.txt",
		"gap2.txt", "treedata2.txt", "Optimalleaforder2.txt");


	cout<<"finish clusclus....."<<endl;
	this->setDataForHeatmap(cc1->features, cc1->optimalleaforder, cc2->optimalleaforder,cc1->num_samples, cc2->num_samples);
	this->setDataForDendrograms(cc1->treedata, cc2->treedata);
	this->creatDataForHeatmap(0.2);	
	
	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		delete datas[i];
	}
	delete datas;

	delete cc1;
	delete cc2;
}

void ProgressionHeatmap::showGraph()
{	
	this->drawPoints1();
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	
	this->lookuptable->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->lookuptable->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->lookuptable->Build();

	int k = 0;
	boost::math::normal N;
	try
	{
		for(int i = 0; i < this->num_samples; i++)
		{
			for(int j = 0; j < this->num_features; j++)
			{
				if(mapdata[num_samples - i - 1][j] <= 0)
					lookuptable->SetTableValue(k++, 0, 1 - cdf(N,mapdata[num_samples - i - 1][j]), 0);
				else
					lookuptable->SetTableValue(k++, cdf(N,mapdata[num_samples - i - 1][j]), 0, 0);
			}
		}
	}
	catch(...)
	{
		cout<<"Boost call failed! please try again!"<<endl;
	}

	this->aPlane->Update(); // Force an update so we can set cell data
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);

	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(lookuptable);
	this->actor->SetMapper(mapper);

	this->SetInteractStyle();
	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->SetBackground(1,1,1);

	try
	{
		this->showDendrogram1();
		this->showDendrogram2();
	}
	catch(...)
	{
		cout<<"Draw dendrogram failed ! Please try again!"<<endl;
	}
	/////////////////////////////////////////////////////////////////////
	//vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	//sphereSource->SetCenter(-2.0, 0.0, 0.0);
	//sphereSource->SetRadius(0.5);
	//sphereSource->Update();
 //
	//vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
 //   vtkSmartPointer<vtkPolyDataMapper>::New();
	//sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
 //
	//vtkSmartPointer<vtkActor> sphereActor =
 //   vtkSmartPointer<vtkActor>::New();
	//sphereActor->SetMapper(sphereMapper);

 //   vtkSmartPointer<vtkBalloonRepresentation> balloonRep =
 //   vtkSmartPointer<vtkBalloonRepresentation>::New();
	//balloonRep->SetBalloonLayoutToImageRight();
 //
	//vtkSmartPointer<vtkBalloonWidget> balloonWidget =
 //   vtkSmartPointer<vtkBalloonWidget>::New();
	//balloonWidget->SetInteractor(this->view->GetInteractor());
	//balloonWidget->SetRepresentation(balloonRep);
	//balloonWidget->AddBalloon(sphereActor,"This is a sphere",NULL);
	//this->view->GetRenderer()->AddActor(sphereActor);

	
   // hoverWidget->SetInteractor(this->view->GetInteractor());
   //hoverWidget->SetTimerDuration(1000);
 
  // Create a callback to listen to the widget's two VTK events
	//hoverCallback = vtkSmartPointer<vtkHoverCallback>::New();
 /*   hoverWidget->AddObserver(vtkCommand::TimerEvent,hoverCallback);
    hoverWidget->AddObserver(vtkCommand::EndInteractionEvent,hoverCallback);*/

	//this->view->GetInteractor()->Initialize();
	/////////////////////////////////////////////////////////////////////

	this->view->Render();
	//hoverWidget->On();/////////////
	this->view->GetInteractor()->Start();
}

void ProgressionHeatmap::SetInteractStyle()
{
	this->theme->SetCellValueRange(0, this->num_samples*this->num_features - 1);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
	this->view->ApplyViewTheme(theme);

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

	this->view->GetInteractor()->RemoveObservers(vtkCommand::LeftButtonPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::LeftButtonReleaseEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->view->GetInteractor()->AddObserver(vtkCommand::LeftButtonPressEvent, selectionCallback2);
	this->view->GetInteractor()->AddObserver(vtkCommand::LeftButtonReleaseEvent, selectionCallback3);
	this->view->GetInteractor()->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);
	//this->view->GetInteractor()->AddObserver(vtkCommand::TimerEvent, hoverCallback);
	//this->view->GetInteractor()->AddObserver(vtkCommand::EndInteractionEvent, hoverCallback);

	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
}
void ProgressionHeatmap::showSimilarMatrixGraph()
{	
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	
	this->lookuptable->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->lookuptable->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->lookuptable->Build();

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			rgb rgb = GetRGBValue( mapdata[num_samples - i - 1][j]);
			lookuptable->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->aPlane->Update(); // Force an update so we can set cell data
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);


	this->theme->SetCellValueRange(0, this->num_samples*this->num_features - 1);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
	this->view->ApplyViewTheme(theme);

	this->view->GetInteractor()->SetPicker(this->myCellPicker);
	this->myCellPicker->SetTolerance(0.004);
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2 =vtkSmartPointer<vtkCallbackCommand>::New();
	selectionCallback2->SetClientData(this);
	selectionCallback2->SetCallback(SelectionCallbackFunction2 );

	vtkSmartPointer<vtkCallbackCommand> selectionCallback3 =vtkSmartPointer<vtkCallbackCommand>::New();
	selectionCallback3->SetClientData(this);
	selectionCallback3->SetCallback(SelectionCallbackFunction3);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::LeftButtonPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::LeftButtonReleaseEvent);
	this->view->GetInteractor()->AddObserver(vtkCommand::LeftButtonPressEvent, selectionCallback2);
	this->view->GetInteractor()->AddObserver(vtkCommand::LeftButtonReleaseEvent, selectionCallback3);

	// Setup actor and mapper
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(lookuptable);
	this->actor->SetMapper(mapper);
	
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();

	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->SetBackground(1,1,1);

	this->view->Render();
	this->view->GetInteractor()->Start();
}

void ProgressionHeatmap::showSPDGraph()
{	
	drawPoints3();
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	int index = 0;

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
		}
    }
	
	this->lookuptable->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->lookuptable->SetTableRange(0, this->num_samples*this->num_features - 1);   
	this->lookuptable->Build();

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			rgb rgb = GetRGBValue( mapdata[num_samples - i - 1][j]);
			lookuptable->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->aPlane->Update(); // Force an update so we can set cell data
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);

	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(lookuptable);
	this->actor->SetMapper(mapper);

	this->SetInteractStyle();
	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->SetBackground(1,1,1);
	this->view->GetRenderer()->GradientBackgroundOff();

	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetNumberOfTableValues( COLOR_MAP_SIZE);
	//scalarbarLut->SetTableRange (-1, 1);
	//scalarbarLut->SetHueRange (0.67, 0);
	//scalarbarLut->SetSaturationRange (1, 1);
	//scalarbarLut->SetValueRange (1, 1);

	for(int index = 0; index < COLOR_MAP_SIZE; index++)
	{
		rgb rgbscalar = COLORMAP[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetMaximumNumberOfColors( COLOR_MAP_SIZE);
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(10);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->SetMaximumWidthInPixels(100);
	scalarBar->SetMaximumHeightInPixels(1000);

	//scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	//scalarBar->GetPositionCoordinate()->SetValue(1, 0.0);
	//scalarBar->SetOrientationToVertical();
	//scalarBar->SetWidth(0.5);
	//scalarBar->SetHeight(2);
 
	this->view->GetRenderer()->AddActor2D(scalarBar);

	this->view->Render();
	this->view->GetInteractor()->Start();
}

void ProgressionHeatmap::setDataForDendrograms(double** treedata1, double** treedata2)
{
	if(treedata1 != NULL)
	{
		this->connect_Data_Tree1 = new double*[this->num_samples-1];
		for(int i = 0; i<this->num_samples - 1; i++)
		{
			this->connect_Data_Tree1[i] = new double[4];
			for(int j = 0; j<4; j++)
				this->connect_Data_Tree1[i][j] = treedata1[i][j];
		}
	}

	if( treedata2 != NULL)
	{
		this->connect_Data_Tree2 = new double*[this->num_features-1];
		for(int i = 0; i<this->num_features - 1; i++)
		{
			this->connect_Data_Tree2[i] = new double[4];
			for(int j = 0; j<4; j++)
				this->connect_Data_Tree2[i][j] = treedata2[i][j];
		}
	}
}

void ProgressionHeatmap::createDataForDendogram1(double powCof)
{
	this->Processed_Coordinate_Data_Tree1.resize(2*(this->num_samples) - 1);
	for(int i = 0; i < 2*(this->num_samples) - 1; i++)
	{
		this->Processed_Coordinate_Data_Tree1[i].resize(4);
	}

	for(int i = 0; i < num_samples; i++)
	{
		Processed_Coordinate_Data_Tree1[i][0] = i;

		for(int k = 0; k < num_samples; k++)
		{
			if(Optimal_Leaf_Order1[k] == i)
			{
				Processed_Coordinate_Data_Tree1[i][2] = (k+0.5)/(double)this->num_samples - 0.5;
			}
		}

		Processed_Coordinate_Data_Tree1[i][1] = -0.5;
		Processed_Coordinate_Data_Tree1[i][3] = 0; 
	}

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

void ProgressionHeatmap::createDataForDendogram2(double powCof)
{
	this->Processed_Coordinate_Data_Tree2.resize(2*(this->num_features) - 1);
	for(int i = 0; i < 2*(this->num_features) - 1; i++)
	{
		this->Processed_Coordinate_Data_Tree2[i].resize(4);
	}

	for(int i = 0; i < num_features; i++)
	{
		Processed_Coordinate_Data_Tree2[i][0] = i;
 
		for(int k = 0; k < num_features; k++)
		{
			if(Optimal_Leaf_Order2[k] == i)
			{
				Processed_Coordinate_Data_Tree2[i][1] = (k+0.5)/(double)this->num_features - 0.5;
			}
		}

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

void ProgressionHeatmap::showDendrogram1()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	dencolors1->SetNumberOfComponents(3);
	dencolors1->SetName("denColors1");
	unsigned char color[3] = {0, 0, 0};

	for(int i=0; i<3*(this->num_samples - 1);i++)
		dencolors1->InsertNextTupleValue(color);

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

		denpoints1->InsertNextPoint(p1);
		denpoints1->InsertNextPoint(p2);
		denpoints1->InsertNextPoint(p3);
		denpoints1->InsertNextPoint(p4);

		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0,0 + i*4);
		line0->GetPointIds()->SetId(1,2 + i*4);
		denlines1->InsertNextCell(line0);

		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0,1 + i*4);
		line1->GetPointIds()->SetId(1,3 + i*4);
		denlines1->InsertNextCell(line1);

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0,2 + i*4);
		line2->GetPointIds()->SetId(1,3 + i*4);
		denlines1->InsertNextCell(line2);
	}

	denlinesPolyData1->SetPoints(denpoints1);
	denlinesPolyData1->SetLines(denlines1);
	denlinesPolyData1->GetCellData()->SetScalars(dencolors1);

	denmapper1->SetInput(denlinesPolyData1);
	denactor1->SetMapper(denmapper1);
	denmapper1->SetScalarRange(0, 3*this->num_samples-1);

	this->view->GetRenderer()->AddActor(denactor1);
	//this->view->Render();
}

void ProgressionHeatmap::showDendrogram2()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	dencolors2->SetNumberOfComponents(3);
	dencolors2->SetName("denColors2");
	unsigned char color[3] = {0, 0, 0};

	for(int i=0; i<3*(this->num_features - 1);i++)
		dencolors2->InsertNextTupleValue(color);

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

		denpoints2->InsertNextPoint(p1);
		denpoints2->InsertNextPoint(p2);
		denpoints2->InsertNextPoint(p3);
		denpoints2->InsertNextPoint(p4);

		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0,0 + i*4);
		line0->GetPointIds()->SetId(1,2 + i*4);
		denlines2->InsertNextCell(line0);

		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0,1 + i*4);
		line1->GetPointIds()->SetId(1,3 + i*4);
		denlines2->InsertNextCell(line1);

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0,2 + i*4);
		line2->GetPointIds()->SetId(1,3 + i*4);
		denlines2->InsertNextCell(line2);
	}

	denlinesPolyData2->SetPoints(denpoints2);
	denlinesPolyData2->SetLines(denlines2);
	denlinesPolyData2->GetCellData()->SetScalars(dencolors2);

	denmapper2->SetInput(denlinesPolyData2);
	denactor2->SetMapper(denmapper2);
	denmapper2->SetScalarRange(0, 3*this->num_features-1);
	this->view->GetRenderer()->AddActor(denactor2);
	//this->view->Render();
}

void ProgressionHeatmap::GetSelecectedIDs()
{
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

	while(iter1 != selectedIDs1.end())
	{
		int index1 = *iter1;

		for(int i = 0; i<this->num_samples; i++)
		{
			if(Optimal_Leaf_Order1[i] == indMapFromVertexToInd.find(index1)->second)
			{
				IDs1[count1++] = i;				
				break;
			}		
		}
		iter1++;
	}

	while(iter2 != selectedIDs2.end())
	{
		int index2 = *iter2;

		for(int i = 0; i<this->num_features; i++)
		{
			if(Optimal_Leaf_Order2[i] == index2)
			{
				IDs2[count2++] = i;				
				break;
			}		
		}
		iter2++;
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
		this->view->Render();
	}
	catch(...)
	{
		cout<<"GetSelecectedIDs failed, please try it again!"<<endl;
	}
}

void ProgressionHeatmap::drawPoints1()
{
	v1->SetNumberOfValues (2*this->num_samples-1 + 2*this->num_features-1);
	for(unsigned int i=0; i<2*this->num_samples-1;i++)
    {
		v1->SetValue (i,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree1[i][1],this->Processed_Coordinate_Data_Tree1[i][2],this->Processed_Coordinate_Data_Tree1[i][3]);
	}
	
	for(unsigned int i=0; i<(2*this->num_features-1);i++)
    {
		v1->SetValue (i+2*this->num_samples-1,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1],this->Processed_Coordinate_Data_Tree2[i][2],this->Processed_Coordinate_Data_Tree2[i][3]);
	}
    this->graph_Layout1->SetPoints(this->points1);
     
    ///////////////coloring/////////////////////
    vertexColors1->SetNumberOfComponents(1);
    vertexColors1->SetName("Color1");
	
	int max_table_values = 2*this->num_samples + 2*this->num_features-2;
	lookupTable1->SetNumberOfTableValues(max_table_values);
    for(unsigned int i=0; i<max_table_values;i++)
    {
		lookupTable1->SetTableValue(i, 0.5, 0.5,0.5); // color the vertices- blue
    }  
    lookupTable1->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(unsigned int j=0;j<2*this->num_samples + 2*this->num_features -2;j++)
    {
		vertexColors1->InsertNextValue(j);
		scales1->InsertNextValue(1);
    }

	this->graph_Layout1->GetVertexData()->AddArray(vertexColors1);
    this->view->AddRepresentationFromInput(graph_Layout1);
    this->view->SetLayoutStrategy("Pass Through");
	this->graph_Layout1->GetVertexData()->AddArray(scales1);
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->view->SetVertexColorArrayName("Color1");
    this->view->ColorVerticesOn();
	this->theme->SetPointLookupTable(lookupTable1);


	vtkSmartPointer<vtkStringArray> vertexIDarrays = vtkSmartPointer<vtkStringArray>::New();
	vertexIDarrays->SetNumberOfComponents(1);
	vertexIDarrays->SetName("vertexIDarrays");

	for( unsigned int i = 0; i < this->num_samples; i++)
	{
		//vertexIDarrays->InsertNextValue("");
		vtkVariant v = i;
		vertexIDarrays->InsertNextValue(v.ToString());
	}
	for( unsigned int i = this->num_samples; i <2*(this->num_samples)-1; i++)
		vertexIDarrays->InsertNextValue("");

	for(unsigned int i = 2*(this->num_samples)-1; i <2*(this->num_samples)-1 + this->num_features; i++)
	{
		//vertexIDarrays->InsertNextValue("");
		vtkIdType id = i - (2*this->num_samples-1) + 1 ;
		vertexIDarrays->InsertNextValue(this->table->GetColumn(id)->GetName ());
	}	

	for(unsigned int i = 2*(this->num_samples)-1 + this->num_features; i <2*(this->num_samples)-1 + 2*(this->num_features)-1; i++)
		vertexIDarrays->InsertNextValue("");

	this->graph_Layout1->GetVertexData()->AddArray(vertexIDarrays);
	this->view->SetVertexLabelVisibility(true);
	this->view->SetVertexLabelArrayName("vertexIDarrays");
	this->view->SetVertexLabelFontSize(25);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);
    this->selectionCallback1->SetCallback ( SelectionCallbackFunction1);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void ProgressionHeatmap::drawPoints2()
{
	v1->SetNumberOfValues (2*this->num_features - 1);

	for(int i = 0; i < 2 * this->num_features - 1; i++)
    {
		v1->SetValue (i,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1],this->Processed_Coordinate_Data_Tree2[i][2],this->Processed_Coordinate_Data_Tree2[i][3]);
	}
    this->graph_Layout1->SetPoints(this->points1);
     
    ///////////////coloring/////////////////////
    vertexColors1->SetNumberOfComponents(1);
    vertexColors1->SetName("Color1");
 
    lookupTable1->SetNumberOfTableValues( 2 * this->num_features - 1);
    for(int i = 0; i< 2 * this->num_features - 1; i++)
    {
		lookupTable1->SetTableValue(i, 0.5, 0.5,0.5); // color the vertices- blue
    }  
    lookupTable1->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(int j = 0; j< 2 * this->num_features - 1; j++)
    {
		vertexColors1->InsertNextValue(j);
		scales1->InsertNextValue(0.5);
    }

	this->graph_Layout1->GetVertexData()->AddArray(vertexColors1);
    this->view->AddRepresentationFromInput(graph_Layout1);
    this->view->SetLayoutStrategy("Pass Through");
	this->graph_Layout1->GetVertexData()->AddArray(scales1);
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->view->SetVertexColorArrayName("Color1");
    this->view->ColorVerticesOn();
	this->theme->SetPointLookupTable(lookupTable1);

    this->selectionCallback1 = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback1->SetClientData(this);
    this->selectionCallback1->SetCallback ( SelectionCallbackFunction);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}

void ProgressionHeatmap::drawPoints3()
{
	v1->SetNumberOfValues (this->num_features);
	for(unsigned int i = 0; i < this->num_features; i++)
    {
		v1->SetValue(i ,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1], this->Processed_Coordinate_Data_Tree2[i][2], this->Processed_Coordinate_Data_Tree2[i][3]);
	}
    this->graph_Layout1->SetPoints(this->points1);
     
    ///////////////coloring/////////////////////
    vertexColors1->SetNumberOfComponents(1);
    vertexColors1->SetName("Color1");
	
	int max_table_values = this->num_features;
	lookupTable1->SetNumberOfTableValues(max_table_values);
    for(unsigned int i = 0; i < max_table_values; i++)
    {
		lookupTable1->SetTableValue(i, 0.5, 0.5, 0.5); 
    }  
    lookupTable1->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(unsigned int j = 0; j < this->num_features; j++)
    {
		vertexColors1->InsertNextValue(j);
		scales1->InsertNextValue(1);
    }

	this->graph_Layout1->GetVertexData()->AddArray(vertexColors1);
    this->view->AddRepresentationFromInput(graph_Layout1);
    this->view->SetLayoutStrategy("Pass Through");
	this->graph_Layout1->GetVertexData()->AddArray(scales1);
    this->view->ScaledGlyphsOn();
    this->view->SetScalingArrayName("Scales1");
    vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->view->SetVertexColorArrayName("Color1");
    this->view->ColorVerticesOn();
	this->theme->SetPointLookupTable(lookupTable1);

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);
	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");
	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");
	for(unsigned int i = 0; i < max_table_values;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(45.0);
	}

	for(unsigned int i = 0; i < max_table_values; i++)
	{
		orient->InsertNextValue(90.0);
		//vtkIdType id = Optimal_Leaf_Order2[i];
		label->InsertNextValue(this->table->GetColumn(i + 1)->GetName());   // discard the ID column
	}	

	pd->SetPoints(points1);
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
	//lmapper->SetShapeToRoundedRect();
	lmapper->SetShapeToNone();
	//lmapper->SetBackgroundColor(1.0, 1.0, 0.7);
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

void ProgressionHeatmap::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	ProgressionHeatmap* heatmapWin = (ProgressionHeatmap*)clientData;

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

void ProgressionHeatmap::SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	ProgressionHeatmap* heatmapWin = (ProgressionHeatmap*)clientData;

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

void ProgressionHeatmap::SelectionCallbackFunction2(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	ProgressionHeatmap* heatmapWin = (ProgressionHeatmap*)clientData;
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
}

void ProgressionHeatmap::SelectionCallbackFunction3(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	ProgressionHeatmap* heatmapWin = (ProgressionHeatmap*)clientData;
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

void ProgressionHeatmap::HandleKeyPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	ProgressionHeatmap* heatmapWin = (ProgressionHeatmap*)clientData;
	char key = heatmapWin->view->GetInteractor()->GetKeyCode();
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
	default:
		break;
	}
}

void ProgressionHeatmap::SetdenSelectedIds1(std::set<long int>& IDs, bool bfirst)
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
			this->Selection->select(selectedIDs1);
		}
		else
			this->Selection->select(interselectedIDs);

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
		this->Selection->clear();
	}
}

void ProgressionHeatmap::reselectIds1(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_samples)
	{
		selectedIDs.insert( indMapFromIndToVertex[id]);
	}
	else
	{
		for (int i = 0; i < this->num_samples - 1; i++)
		{
			if(id == connect_Data_Tree1[i][3])
			{
				this->dencolors1->SetValue((id - this->num_samples)*9, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 1, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 2, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 3, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 4, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 5, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 6, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 7, 255);
				this->dencolors1->SetValue((id - this->num_samples)*9 + 8, 255);
				this->reselectIds1(selectedIDs, connect_Data_Tree1[i][0]);
				this->reselectIds1(selectedIDs, connect_Data_Tree1[i][1]);
			}
		}
	}
}
void ProgressionHeatmap::reselectIds2(std::set<long int>& selectedIDs2, long int id)
{
	if(id  < this->num_features)
	{
		selectedIDs2.insert( id);
	}
	else
	{
		for (int i = 0; i < this->num_features - 1; i++)
		{
			if(id == connect_Data_Tree2[i][3])
			{
				this->dencolors2->SetValue((id - this->num_features)*9, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 1, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 2, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 3, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 4, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 5, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 6, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 7, 255);
				this->dencolors2->SetValue((id - this->num_features)*9 + 8, 255);
				this->reselectIds2(selectedIDs2, connect_Data_Tree2[i][0]);
				this->reselectIds2(selectedIDs2, connect_Data_Tree2[i][1]);
			}
		}
	}
}

void ProgressionHeatmap::computeselectedcells()
{
	this->r1 = id1/this->num_features;
	this->r2 = id2/this->num_features;
	this->c1 = id1%this->num_features;
	this->c2 = id2%this->num_features;

	for(int i = 0; i <= r1 - r2; i++)
	{
		for(int j = 0; j <= c2 - c1; j++)
		{
			ids->InsertNextValue(id2 - j + this->num_features*i);
		}
	}
}

void ProgressionHeatmap::setselectedCellIds()
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

void ProgressionHeatmap::GetSelRowCol(int &r1, int &c1, int &r2, int &c2)
{
	r1 = this->r1;
	r2 = this->r2;
	c1 = this->c1;
	c2 = this->c2;
}

void ProgressionHeatmap::SetSelRowCol(int r1, int c1, int r2, int c2)
{
	this->r1 = r1;
	this->r2 = r2;
	this->c1 = c1;
	this->c2 = c2;
	
	setselectedCellIds();	
}

void ProgressionHeatmap::closeEvent(QCloseEvent *event)
{
	mainQTRenderWidget.close();
}
