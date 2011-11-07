#include "HeatmapWindow.h"
 
vtkStandardNewMacro(MouseInteractorStyle);

Heatmap::Heatmap(QWidget *parent)
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
	///////////////////////////////////////////////
	this->v1 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->v2 = vtkSmartPointer<vtkIdTypeArray>::New();
	this->graph_Layout1 = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->graph_Layout2 = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points1 = vtkSmartPointer<vtkPoints>::New();
	this->points2 = vtkSmartPointer<vtkPoints>::New();
	this->vertexColors1 = vtkSmartPointer<vtkIntArray>::New();
	this->vertexColors2 = vtkSmartPointer<vtkIntArray>::New();
	this->lookupTable1 = vtkSmartPointer<vtkLookupTable>::New();
	this->lookupTable2 = vtkSmartPointer<vtkLookupTable>::New();
}

Heatmap::~Heatmap()
{
}

void Heatmap::setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features)
{
	this->mapdata = features;
	this->Optimal_Leaf_Order1 = optimalleaforder1;
	this->Optimal_Leaf_Order2 = optimalleaforder2;
	this->num_samples = num_samples;
	this->num_features = num_features;
}

void Heatmap::creatDataForHeatmap()
{
	this->scaleData();

	vector<double > temp;
	temp.resize(num_features);
	///////
	double** tempdata;
	tempdata = new double*[this->num_samples];
	for(int i = 0; i < this->num_samples; i++)
		tempdata[i] = new double[this->num_features];
	/////

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

		for(int j = 0; j < num_features; j++)
			tempdata[i][j] = (temp[j] - mean)/std;
	}
	for(int i = 0; i < this->num_samples; i++)
		mapdata[this->num_samples - i - 1] = tempdata[Optimal_Leaf_Order1[i]]; 
//////////////////////////////////////////////////////////for debug
	const char* filename = "heatmapdata";
	FILE *fp = fopen(filename,"w");
	for(int i=0; i<this->num_samples; i++)
	{
		for(int j=0; j<this->num_features; j++)
			fprintf(fp,"%14.7e\t",mapdata[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
//////////////////////////////////////////////
	this->createDataForDendogram1();
	this->createDataForDendogram2();
}

void Heatmap::scaleData()
{
	for(int i = 0; i<this->num_features; i++)
	{
		double sum = 0.0;
		for(int j = 0; j<this->num_samples; j++)
		{
			sum += mapdata[j][i]*mapdata[j][i];
		}
		sum = sqrt(sum);

		for(int j = 0; j<this->num_samples; j++)
		{
			mapdata[j][i] /= sum;
		}
	}
}

void Heatmap::setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL)
{
	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;

	if(!sels2)
		this->Selection2 = new ObjectSelection();
	else
		this->Selection2 = sels2;

	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));
	//connect(Selection2, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));
}

void Heatmap::showGraph()
{	
	this->drawPoints1();
	this->drawPoints2();
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

	this->aPlane->Update(); // Force an update so we can set cell data
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);

 /////////////////////////////////////////////////////////////////////
	this->theme->SetCellValueRange(0, this->num_samples*this->num_features - 1);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
	this->view->ApplyViewTheme(theme);
////////////////////////////////////////////////////////////

	// Setup actor and mapper
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(lookuptable);
	this->actor->SetMapper(mapper);
	
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();

	/*vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	/style->SetDefaultRenderer(this->view->GetRenderer());
	/style->Data = this->aPlane->GetOutput();
	/style->hm = this;
	/this->view->GetInteractor()->SetInteractorStyle(style);*/

	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->SetBackground(1,1,1);

	this->showDendrogram1();
	this->showDendrogram2();
	this->view->Render();
	this->view->GetInteractor()->Start();
}

void Heatmap::setDataForDendrograms(double** treedata1, double** treedata2)
{
	this->connect_Data_Tree1 = treedata1;
	this->connect_Data_Tree2 = treedata2;
}

void Heatmap::createDataForDendogram1()
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
		//Processed_Coordinate_Data_Tree[i][4] = 0; 
	}

	for(int i = 0; i < num_samples-1; i++)
	{
		connect_Data_Tree1[i][2] = pow(connect_Data_Tree1[i][2], 0.2);
		connect_Data_Tree1[i][2] /= pow(connect_Data_Tree1[num_samples-2][2], 0.2);
	}
	connect_Data_Tree1[num_samples-2][2] = 1;

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
		//Processed_Coordinate_Data_Tree[i][4] = 0;
	}
}

void Heatmap::createDataForDendogram2()
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
		//Processed_Coordinate_Data_Tree[i][4] = 0; 
	}

	for(int i = 0; i < num_features-1; i++)
	{
		connect_Data_Tree2[i][2] = pow(connect_Data_Tree2[i][2], 0.2);
		connect_Data_Tree2[i][2] /= pow(connect_Data_Tree2[num_features - 2][2], 0.2);
	}
	connect_Data_Tree2[num_features - 2][2] = 1;

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
		//Processed_Coordinate_Data_Tree[i][4] = 0;
	}
}

void Heatmap::showDendrogram1()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

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

		vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
		lineSource1->SetPoint1(p1);
		lineSource1->SetPoint2(p3);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper1->SetInputConnection(lineSource1->GetOutputPort());

		vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
		actor1->GetProperty()->SetColor(0.5,0.7,0); //set the colours of the line of the dendrogram
		actor1->GetProperty()->SetLineWidth(0.5);
		actor1->SetMapper(mapper1);
		this->view->GetRenderer()->AddActor(actor1);

		vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
		lineSource2->SetPoint1(p2);
		lineSource2->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper2->SetInputConnection(lineSource2->GetOutputPort());

		vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
		actor2->SetMapper(mapper2);
		actor2->GetProperty()->SetColor(0.5,0.7,0);
		actor2->GetProperty()->SetLineWidth(0.5);
		this->view->GetRenderer()->AddActor(actor2);

		vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
		lineSource3->SetPoint1(p3);
		lineSource3->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper3->SetInputConnection(lineSource3->GetOutputPort());

		vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
		actor3->GetProperty()->SetColor(0.5,0.7,0);
		actor3->GetProperty()->SetLineWidth(0.5);
		actor3->SetMapper(mapper3);
		this->view->GetRenderer()->AddActor(actor3);

		//this->view->Render();
	}
}

void Heatmap::showDendrogram2()
{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

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

		vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
		lineSource1->SetPoint1(p1);
		lineSource1->SetPoint2(p3);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper1->SetInputConnection(lineSource1->GetOutputPort());

		vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
		actor1->GetProperty()->SetColor(0.5,0.7,0); //set the colours of the line of the dendrogram
		actor1->GetProperty()->SetLineWidth(0.5);
		actor1->SetMapper(mapper1);
		this->view->GetRenderer()->AddActor(actor1);

		vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
		lineSource2->SetPoint1(p2);
		lineSource2->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper2->SetInputConnection(lineSource2->GetOutputPort());

		vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
		actor2->SetMapper(mapper2);
		actor2->GetProperty()->SetColor(0.5,0.7,0);
		actor2->GetProperty()->SetLineWidth(0.5);
		this->view->GetRenderer()->AddActor(actor2);

		vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
		lineSource3->SetPoint1(p3);
		lineSource3->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper3->SetInputConnection(lineSource3->GetOutputPort());

		vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
		actor3->GetProperty()->SetColor(0.5,0.7,0);
		actor3->GetProperty()->SetLineWidth(0.5);
		actor3->SetMapper(mapper3);
		this->view->GetRenderer()->AddActor(actor3);

		//this->view->Render();

	}
}

void Heatmap::GetSelecectedIDs()
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
			if(Optimal_Leaf_Order1[i] == index1)
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
	selectedActor->GetProperty()->SetEdgeColor(1,0,1);
	selectedActor->GetProperty()->SetLineWidth(0.5);
	
	this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
	this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
}

/*void Heatmap::GetSelecectedIDs2()
{
	std::set<long int> selectedIDs2 = this->Selection->getSelections();
}*/

void Heatmap::drawPoints1()
{
	v1->SetNumberOfValues (2*this->num_samples-1);
	for(int i=0; i<((2*this->num_samples)-1);i++)
    {
		v1->SetValue (i,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree1[i][1],this->Processed_Coordinate_Data_Tree1[i][2],this->Processed_Coordinate_Data_Tree1[i][3]);
	}
	for(int i=0; i<((2*this->num_features)-1);i++)
    {
		v1->SetValue (i,graph_Layout1->AddVertex());
		this->points1->InsertNextPoint(this->Processed_Coordinate_Data_Tree2[i][1],this->Processed_Coordinate_Data_Tree2[i][2],this->Processed_Coordinate_Data_Tree2[i][3]);
	}
    this->graph_Layout1->SetPoints(this->points1);
     
    ///////////////coloring/////////////////////
    vertexColors1->SetNumberOfComponents(1);
    vertexColors1->SetName("Color1");
 
    lookupTable1->SetNumberOfTableValues(2*(this->num_samples)-1 + 2*(this->num_features)-1);
    for(int i=0; i<(2*(this->num_samples)-1) + (2*(this->num_features)-1);i++)
    {
		lookupTable1->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
    }  
    lookupTable1->Build();
   
	vtkSmartPointer<vtkFloatArray> scales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales1->SetNumberOfComponents(1);
	scales1->SetName("Scales1");

    for(int j=0;j<(2*(this->num_samples)-1) + (2*(this->num_features)-1);j++)
    {
		vertexColors1->InsertNextValue(j);
		scales1->InsertNextValue(1.3);
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
    this->selectionCallback1->SetCallback ( SelectionCallbackFunction1);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback1);
}
void Heatmap::drawPoints2()
{
}

void Heatmap::SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	cout<<"vtkselectionhaha.............................."<<endl;
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	Heatmap* heatmapWin = (Heatmap*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;

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

		heatmapWin->SetdenSelectedIds1( IDs);
	}
}


void Heatmap::SetdenSelectedIds1(std::set<long int>& IDs)
{
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;
	std::set<long int>::iterator it;
	long int id;
	if( IDs.size() > 0)
	{
		for(it=IDs.begin(); it != IDs.end(); it++ )
		{
			id = *it;
			cout<<"id ====="<<id<<endl;
			if(id < 2*this->num_samples-1)
			{
				reselectIds1(selectedIDs1, id);
			}
			else
				reselectIds2(selectedIDs2, id - (2*this->num_samples-1));
		}
	}

	if(selectedIDs1.size() > 0)	
	{
		for(int i = 0; i<this->num_features; i++)
		{
			selectedIDs2.insert(i);
		}
		this->Selection2->select(selectedIDs2);
		this->Selection->select(selectedIDs1);

	}

	else if(selectedIDs2.size() > 0)	
	{
		this->Selection2->select(selectedIDs2);
		for(int i = 0; i<this->num_samples; i++)
		{
			selectedIDs1.insert(i);
		}
		this->Selection->select(selectedIDs1);
	}

	else
	{
		this->Selection2->clear();
		this->Selection->clear();
	}
}

void Heatmap::reselectIds1(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_samples)
	{
		selectedIDs.insert(id);
	}
	else
	{
		for (int i = 0; i < this->num_samples - 1; i++)
		{
			if(id == connect_Data_Tree1[i][3])
			{
				this->reselectIds1(selectedIDs, connect_Data_Tree1[i][0]);
				this->reselectIds1(selectedIDs, connect_Data_Tree1[i][1]);
			}
		}
	}
}
void Heatmap::reselectIds2(std::set<long int>& selectedIDs2, long int id)
{
	if(id  < this->num_features)
	{
		selectedIDs2.insert(id);
	}
	else
	{
		for (int i = 0; i < this->num_features - 1; i++)
		{
			if(id == connect_Data_Tree2[i][3])
			{
				this->reselectIds2(selectedIDs2, connect_Data_Tree2[i][0]);
				this->reselectIds2(selectedIDs2, connect_Data_Tree2[i][1]);
			}
		}
	}
}
MouseInteractorStyle::MouseInteractorStyle()
{
	selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selectedActor = vtkSmartPointer<vtkActor>::New();
}


void MouseInteractorStyle::OnLeftButtonDown()
{
	// Get the location of the click (in window coordinates)
	int* pos = this->GetInteractor()->GetEventPosition();
 
	vtkSmartPointer<vtkPicker> picker = vtkSmartPointer<vtkPicker>::New();
	picker->SetTolerance(0.0005);
 
	// Pick from this location.
	picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
	double* worldPosition = picker->GetPickPosition();
 
	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();
		std::cout << "Cell id is: " << cellpicker->GetCellId() << std::endl;

		if(cellpicker->GetCellId() != -1)
		{
			std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1] << " " << worldPosition[2] << endl;
	 
			id1 = cellpicker->GetCellId();
		}
	}

	else if(worldPosition[0]<-0.5)
	{
		vtkSmartPointer<vtkPointPicker> pointpicker = vtkSmartPointer<vtkPointPicker>::New();
		pointpicker->SetTolerance(0.005);
 
		// Pick from this location.
		pointpicker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
 
		double* worldPosition = pointpicker->GetPickPosition();
		std::cout << "Point id is: " << pointpicker->GetPointId() << std::endl;

		if(pointpicker->GetPointId() != -1)
		{
			std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1] << " " << worldPosition[2] << endl;	 
			id1 = pointpicker->GetPointId();
		}
	}
	//cout<<"haha"<<endl;
	vtkInteractorStyleImage::OnLeftButtonDown();
}

void MouseInteractorStyle::OnLeftButtonUp()
{
	// Get the location of the click (in window coordinates)
	int* pos = this->GetInteractor()->GetEventPosition();
 
	vtkSmartPointer<vtkPicker> picker = vtkSmartPointer<vtkPicker>::New();
	picker->SetTolerance(0.0005);
 
	// Pick from this location.
	picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer()); 
	double* worldPosition = picker->GetPickPosition();
 
	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();
		std::cout << "Cell id is: " << cellpicker->GetCellId() << std::endl;
		if(cellpicker->GetCellId() != -1)
		{
			std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1] << " " << worldPosition[2] << endl;
	 
			this->ids = vtkSmartPointer<vtkIdTypeArray>::New();
			this->ids->SetNumberOfComponents(1);
			id2 = cellpicker->GetCellId();
			this->computeselectedcells();
	//////////////////////////////////////////////////////////////////////////
			/*selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
			selectionNode->SetFieldType(vtkSelectionNode::CELL);
			selectionNode->SetContentType(vtkSelectionNode::INDICES);
			selectionNode->SetSelectionList(ids);

			selection = vtkSmartPointer<vtkSelection>::New();
			selection->AddNode(selectionNode);
		 
			extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
			extractSelection->SetInput(0, this->Data);
			extractSelection->SetInput(1, selection);
			extractSelection->Update();
		 
			selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
			selected->ShallowCopy(extractSelection->GetOutput());
		 
			std::cout << "There are " << selected->GetNumberOfPoints()
					  << " points in the selection." << std::endl;
			std::cout << "There are " << selected->GetNumberOfCells()
					  << " cells in the selection." << std::endl; 
		 
			selectedMapper->SetInputConnection(selected->GetProducerPort());
		 
			selectedActor->SetMapper(selectedMapper);
			selectedActor->GetProperty()->EdgeVisibilityOn();
			selectedActor->GetProperty()->SetEdgeColor(1,0,1);
			selectedActor->GetProperty()->SetLineWidth(1);
		 
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);*/

			this->setselectedIds();
		}
	}
	//cout<<"haha"<<endl;
	vtkInteractorStyleImage::OnLeftButtonUp();
}

void MouseInteractorStyle::computeselectedcells()
{
	this->r1 = id1/this->hm->num_features;
	this->r2 = id2/this->hm->num_features;
	this->c1 = id1%this->hm->num_features;
	this->c2 = id2%this->hm->num_features;

	cout<<"r1 = "<<r1<<endl;
	cout<<"r2 = "<<r2<<endl;
	cout<<"c1 = "<<c1<<endl;
	cout<<"c2 = "<<c2<<endl;


	for(int i = 0; i<=r1-r2; i++)
	{
		for(int j = 0; j<=c2-c1; j++)
		{
			ids->InsertNextValue(id2 - j + this->hm->num_features*i);
		}
	}
}

void MouseInteractorStyle::setselectedIds()
{
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;

	for(int i = r2; i<=r1; i++)
	{
		selectedIDs1.insert(this->hm->Optimal_Leaf_Order1[i]);
	}
	for(int j = c1; j<=c2; j++)
	{		
		selectedIDs2.insert(this->hm->Optimal_Leaf_Order2[j]);
	}

	this->hm->Selection2->select(selectedIDs2);
	this->hm->Selection->select(selectedIDs1);
}
