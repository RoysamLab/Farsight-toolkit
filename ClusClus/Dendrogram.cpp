#include"Dendrogram.h"

Dendrogram::Dendrogram(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;
	this->Optimal_Leaf_Order = NULL;
	this->connect_Data_Tree = NULL;
	//this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	//this->theme = vtkSmartPointer<vtkViewTheme>::New();
	//this->graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
	//this->points = vtkSmartPointer<vtkPoints>::New();
	//this->vertexColors = vtkSmartPointer<vtkIntArray>::New();
	//this->lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	//this->v = vtkSmartPointer<vtkIdTypeArray>::New();
}

Dendrogram::~Dendrogram()
{
	if(this->Optimal_Leaf_Order)
	{
		delete this->Optimal_Leaf_Order;
		this->Optimal_Leaf_Order = NULL;
	}
	if(this->connect_Data_Tree)
	{
		for(int i = 0; i<this->num_samples - 1; i++)
			delete this->connect_Data_Tree[i];
		delete this->connect_Data_Tree;
		this->connect_Data_Tree = NULL;
	}
}

void Dendrogram::setTreeData(int numsamples, double** treedata, int* optimalleaforder)
{
	this->num_samples = numsamples;
	this->connect_Data_Tree = new double*[this->num_samples-1];
	for(int i = 0; i<this->num_samples - 1; i++)
	{
		this->connect_Data_Tree[i] = new double[4];
		for(int j = 0; j<4; j++)
			this->connect_Data_Tree[i][j] = treedata[i][j];
	}
	//this->connect_Data_Tree = treedata;
	this->Optimal_Leaf_Order = new int[this->num_samples];
	for(int i = 0; i<this->num_samples; i++)
		this->Optimal_Leaf_Order[i] = optimalleaforder[i];	
}

void Dendrogram::createDataForDendogram()
{
	this->Processed_Coordinate_Data_Tree.resize(2*(this->num_samples) - 1);
	for(int i = 0; i < 2*(this->num_samples) - 1; i++)
	{
		this->Processed_Coordinate_Data_Tree[i].resize(4);
	}

	for(int i = 0; i < num_samples; i++)
	{
		Processed_Coordinate_Data_Tree[i][0] = i;

		for(int k = 0; k < num_samples; k++)
		{
			if(Optimal_Leaf_Order[k] == i)
			{
				Processed_Coordinate_Data_Tree[i][1] = k+1;
			}
		}

		Processed_Coordinate_Data_Tree[i][2] = 0;
		Processed_Coordinate_Data_Tree[i][3] = 0; 
		//Processed_Coordinate_Data_Tree[i][4] = 0; 
	}

	for(int i = 0; i < num_samples-1; i++)
	{
		connect_Data_Tree[i][2] = pow(connect_Data_Tree[i][2], this->powCof);
	}

	for(int i = num_samples ; i < 2*num_samples - 1; i++)
	{
		Processed_Coordinate_Data_Tree[i][0] = i;

		for(int k = 0; k < num_samples -1 ; k++)
		{
			if(i == connect_Data_Tree[k][3])
			{
				double temp1, temp2;
				temp1 = connect_Data_Tree[k][0];
				temp2 = connect_Data_Tree[k][1];
				Processed_Coordinate_Data_Tree[i][1] = (Processed_Coordinate_Data_Tree[temp1][1] + Processed_Coordinate_Data_Tree[temp2][1])/2;
				Processed_Coordinate_Data_Tree[i][2] = connect_Data_Tree[k][2];
			}
		}

		Processed_Coordinate_Data_Tree[i][3] = 0; 
		//Processed_Coordinate_Data_Tree[i][4] = 0;
	}
}
void Dendrogram::showGraph()
{
	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	this->graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->vertexColors = vtkSmartPointer<vtkIntArray>::New();
	this->lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	this->v = vtkSmartPointer<vtkIdTypeArray>::New();
	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	v->SetNumberOfValues (2*this->num_samples-1);
	for(int i=0; i<((2*this->num_samples)-1);i++)
    {
		v->SetValue (i,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree[i][1],this->Processed_Coordinate_Data_Tree[i][2],this->Processed_Coordinate_Data_Tree[i][3]);
	}
    graph_Layout->SetPoints(this->points);
     
    ///////////////coloring/////////////////////
    vertexColors->SetNumberOfComponents(1);
    vertexColors->SetName("Color");
 
    lookupTable->SetNumberOfTableValues(2*(this->num_samples)-1);
    for(int i=0; i<2*(this->num_samples)-1;i++)
    {
		lookupTable->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
    }  
    lookupTable->Build();
   
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales->SetNumberOfComponents(1);
	scales->SetName("Scales");

    for(int j=0;j<2*(this->num_samples)-1;j++)
    {
		vertexColors->InsertNextValue(j);
		scales->InsertNextValue(1.3);
    }

	graph_Layout->GetVertexData()->AddArray(vertexColors);
    this->graphLayoutView->SetRepresentationFromInput(graph_Layout);
    this->graphLayoutView->SetLayoutStrategy("Pass Through");
	graph_Layout->GetVertexData()->AddArray(scales);
    this->graphLayoutView->ScaledGlyphsOn();
    this->graphLayoutView->SetScalingArrayName("Scales");
    vtkRenderedGraphRepresentation::SafeDownCast(this->graphLayoutView->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

    
	//this->graphLayoutView->Render();
	this->graphLayoutView->SetVertexColorArrayName("Color");
    this->graphLayoutView->ColorVerticesOn();
	this->graphLayoutView->ResetCamera();
    //this->theme->SetBackgroundColor(1,1,1);   /// to set the background colour of the dendrogram window
	this->theme->SetPointLookupTable(lookupTable);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1);
    this->graphLayoutView->ApplyViewTheme(theme);

	std::cout<<"the no of vertices"<<graph_Layout->GetNumberOfVertices()<<std::endl;

	this->mainQTRenderWidget.SetRenderWindow(graphLayoutView->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();

    this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback->SetClientData(this);
    this->selectionCallback->SetCallback ( SelectionCallbackFunction);
    graphLayoutView->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback);
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
 
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	for(int i=0; i<num_samples-1;i++)
	{
		double temp1 = this->connect_Data_Tree[i][0];
        double temp2 = this->connect_Data_Tree[i][1];

		for(int j=0; j<(2*(this->num_samples))-1; j++)
        {
            if(this->Processed_Coordinate_Data_Tree[j][0]==temp1)
			{
				p1[0]=this->Processed_Coordinate_Data_Tree[j][1];
				p1[1]=this->Processed_Coordinate_Data_Tree[j][2];
				p1[2]=this->Processed_Coordinate_Data_Tree[j][3];
            }   
            if(this->Processed_Coordinate_Data_Tree[j][0]==temp2)
			{
                p2[0]=this->Processed_Coordinate_Data_Tree[j][1];
                p2[1]=this->Processed_Coordinate_Data_Tree[j][2];
                p2[2]=this->Processed_Coordinate_Data_Tree[j][3];
			}                             
        }
	   
        p3[0]=p1[0];
        p3[1]=this->connect_Data_Tree[i][2];
        p3[2]=p1[2];

        p4[0]=p2[0];
        p4[1]=this->connect_Data_Tree[i][2];
        p4[2]=p2[2];

		vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
		lineSource1->SetPoint1(p1);
		lineSource1->SetPoint2(p3);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper1->SetInputConnection(lineSource1->GetOutputPort());

		vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
		actor1->GetProperty()->SetColor(0.5,0.7,0); //set the colours of the line of the dendrogram
		actor1->GetProperty()->SetLineWidth(1.5);
		actor1->SetMapper(mapper1);
		this->graphLayoutView->GetRenderer()->AddActor(actor1);

		vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
		lineSource2->SetPoint1(p2);
		lineSource2->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper2->SetInputConnection(lineSource2->GetOutputPort());

		vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
		actor2->SetMapper(mapper2);
		actor2->GetProperty()->SetColor(0.5,0.7,0);
		actor2->GetProperty()->SetLineWidth(1.5);
		this->graphLayoutView->GetRenderer()->AddActor(actor2);

		vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
		lineSource3->SetPoint1(p3);
		lineSource3->SetPoint2(p4);
	   
		vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper3->SetInputConnection(lineSource3->GetOutputPort());

		vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
		actor3->GetProperty()->SetColor(0.5,0.7,0);
		actor3->GetProperty()->SetLineWidth(1.5);
		actor3->SetMapper(mapper3);
		/////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////
		this->graphLayoutView->GetRenderer()->AddActor(actor3);


		this->graphLayoutView->Render();

	}
	
    this->graphLayoutView->GetRenderer()->Render();
	this->graphLayoutView->GetInteractor()->Start();
	
}

void Dendrogram::runClusclus()
{
	double** datas;
	vtkVariant temp; 

	datas = new double*[this->table->GetNumberOfRows()];

	std::cout<<this->table->GetNumberOfRows()<<endl;
	std::cout<<this->table->GetNumberOfColumns()<<endl;

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		datas[i] = new double[this->table->GetNumberOfColumns() - 1 + 2];
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

	this->setTreeData(cc1->num_samples, cc1->treedata, cc1->optimalleaforder);
	this->createDataForDendogram();
	//this->dendro1->showGraph();

	for (int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		delete datas[i];
	}
	delete datas;
}
void Dendrogram::SetSelectedIds(std::set<long int>& IDs)
{
	std::set<long int> selectedIDs;
	std::set<long int>::iterator it;
	long int id;
	if( IDs.size() > 0)
	{
		for(it=IDs.begin(); it != IDs.end(); it++ )
		{
			id = *it;
			reselectIds(selectedIDs, id);
		}
		
		this->Selection->select(selectedIDs);
	}
	else
	{
		this->Selection->clear();
	}
}

void Dendrogram::reselectIds(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_samples)
	{
		selectedIDs.insert(id);
	}
	else
	{
		for (int i = 0; i < this->num_samples - 1; i++)
		{
			if(id == connect_Data_Tree[i][3])
			{
				this->reselectIds(selectedIDs, connect_Data_Tree[i][0]);
				this->reselectIds(selectedIDs, connect_Data_Tree[i][1]);
			}
		}
	}
}

void Dendrogram::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	Dendrogram* DendrogramWin = (Dendrogram*)clientData;

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

		DendrogramWin->SetSelectedIds( IDs);
	}
}


void Dendrogram::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, double powCof)
{
	if(table)
		this->table = table;

	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;
	this->powCof = powCof;

	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));

	//this->showGraph();
}

void Dendrogram::GetSelecectedIDs()
{
	std::set<long int> selectedIDs = this->Selection->getSelections();

	for( int i = 0; i < this->num_samples; i++)
	{
		if (selectedIDs.find(i) != selectedIDs.end())
		{
			this->lookupTable->SetTableValue(i, 1.0, 0, 0); // color the vertices- red
		}
		else
		{
			this->lookupTable->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
		}
	}
	this->lookupTable->Build();
		
	this->graphLayoutView->ColorVerticesOn(); 

	this->mainQTRenderWidget.GetRenderWindow()->Render();
	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New() ;
	theme->SetPointLookupTable(lookupTable);
	
	this->graphLayoutView->ApplyViewTheme(theme);
}
