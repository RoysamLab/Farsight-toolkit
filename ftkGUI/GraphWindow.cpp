#include "GraphWindow.h"
#include <vtkAnnotationLink.h>
#include <vtkDataRepresentation.h>
#include <vtkSelectionNode.h>
#include <vtkIdTypeArray.h>
#include <vtkSelection.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkDataSetAttributes.h>

GraphWindow::GraphWindow(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;// = new QVTKWidget;
	this->TTG = vtkSmartPointer<vtkTableToGraph>::New();
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->observerTag = 0;
	this->lookupTable = vtkSmartPointer<vtkLookupTable>::New();
}

GraphWindow::~GraphWindow()
{
}

void GraphWindow::setQtModels(QItemSelectionModel *mod)
{
}

void GraphWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels)
{
	this->dataTable = table;
	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		this->rootID.insert(this->dataTable->GetValue( i, 0).ToInt());
	}
	if(!sels)
		this->selection = new ObjectSelection();
	else
		this->selection = sels;
	connect( this->selection, SIGNAL( changed()), this, SLOT( UpdateGraphView()));
}
	
void GraphWindow::SetGraphTable(vtkSmartPointer<vtkTable> table)
{
	//graphTable->Dump(8);	//debug dump
	this->TTG->ClearLinkVertices();
	this->TTG->SetInput(0, table);
	this->TTG->AddLinkEdge("Source", "Target"); 

	
	this->theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	this->theme->SetLineWidth(5);
	this->theme->SetCellOpacity(0.9);
	this->theme->SetCellAlphaRange(0.5,0.5);
	this->theme->SetPointSize(10);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1); 

	
	this->view->AddRepresentationFromInputConnection(TTG->GetOutputPort());
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName("Distance");
	this->view->SetLayoutStrategyToForceDirected();
	this->view->SetVertexLabelArrayName("label");
	this->view->VertexLabelVisibilityOn();
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2)
{
	//graphTable->Dump(8);	//debug dump
	this->TTG->ClearLinkVertices();
	this->TTG->SetInput(0, table);
	this->TTG->AddLinkEdge(ID1.c_str(), ID2.c_str()); 

	
	this->theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	this->theme->SetLineWidth(5);
	this->theme->SetCellOpacity(0.9);
	this->theme->SetCellAlphaRange(0.5,0.5);
	this->theme->SetPointSize(10);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1); 

	this->view->AddRepresentationFromInputConnection(TTG->GetOutputPort());
	/*this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName("Distance");*/
	this->view->SetLayoutStrategyToForceDirected();
	this->view->SetVertexLabelArrayName("label");
	this->view->VertexLabelVisibilityOn();
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel)
{
	//graphTable->Dump(8);	//debug dump

	this->TTG->ClearLinkVertices();
	this->TTG->SetInput(0, table);
	this->TTG->AddLinkEdge(ID1.c_str(), ID2.c_str()); 

	vtkSmartPointer<vtkMutableDirectedGraph> graph = vtkMutableDirectedGraph::New();
	for( int i = 0; i <  this->dataTable->GetNumberOfRows(); i++)
	{
		int vertexID = graph->AddVertex();
	}

	for( int i = 0; i < table->GetNumberOfRows(); i++)
	{
		graph->AddEdge( table->GetValue(i,0).ToInt(), table->GetValue(i,1).ToInt());
	}

	std::ofstream ofs("mutablegraph.txt", std::ofstream::out);
	graph->Print( ofs);
	ofs.close();
	
	this->theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	this->theme->SetLineWidth(5);
	this->theme->SetCellOpacity(0.9);
	this->theme->SetCellAlphaRange(0.8,0.8);
	this->theme->SetPointSize(8);
	this->theme->SetSelectedCellColor(1,0,0);
	this->theme->SetSelectedPointColor(1,0,0); 

	vtkSmartPointer<vtkIntArray> vertexColors = vtkSmartPointer<vtkIntArray>::New();
	vertexColors->SetNumberOfComponents(table->GetNumberOfRows());
	vertexColors->SetName("Color");
	
	this->lookupTable->SetNumberOfTableValues( this->dataTable->GetNumberOfRows());
	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		vertexColors->InsertNextValue( i);
		this->lookupTable->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
    }
	lookupTable->Build();

	graph->GetVertexData()->AddArray(vertexColors);
	this->view->AddRepresentationFromInput( graph);
	this->view->SetVertexColorArrayName("Color");
	this->view->ColorVerticesOn(); 

    theme->SetPointLookupTable(lookupTable);
    theme->SetBackgroundColor(0,0,0); 
	this->view->ApplyViewTheme(theme);

	this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName(edgeLabel.c_str());
	this->view->SetLayoutStrategyToForceDirected();
	this->view->SetVertexLabelArrayName(ID1.c_str());
	this->view->VertexLabelVisibilityOn();
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::ShowGraphWindow()
{
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
	view->ResetCamera();
	view->Render();

    this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback->SetClientData(this);
    this->selectionCallback->SetCallback ( SelectionCallbackFunction);
	vtkAnnotationLink *link = view->GetRepresentation()->GetAnnotationLink();
	this->observerTag = link->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);

	view->GetInteractor()->Start();
}

void GraphWindow::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	GraphWindow* graphWin = (GraphWindow*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;

	if( selection->GetNode(0))
	{
		if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(0);
		}
		else if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
		{
			edges = selection->GetNode(0);
		}
	}

	if( selection->GetNode(1))
	{
		if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(1);
		}
		else if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
		{
			edges = selection->GetNode(1);
		}
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
		graphWin->SetSelectedIds( IDs);
	}

	graphWin->view->GetRenderer()->Render();
}

ObjectSelection * GraphWindow::GetSelection()
{
	return selection;
}

void GraphWindow::SetSelectedIds(std::set<long int>& IDs)
{
	if( IDs.size() > 0)
	{
		std::set<long int> selectedIDs;
		std::set<long int>::iterator iter = this->rootID.begin();
		std::set<long int>::iterator IdIter;
		long int count = 0;
		for( IdIter = IDs.begin(); IdIter != IDs.end() && iter != this->rootID.end(); IdIter++)
		{
			while( count < *IdIter)
			{
				iter++;
				count++;
			}
			selectedIDs.insert( *iter);
		}
		this->selection->select( selectedIDs);
		UpdataLookupTable( this->selection->getSelections());
	}
	this->view->GetRenderer()->Render();
}

void GraphWindow::UpdataLookupTable( std::set<long int>& IDs)
{
	std::set<long int> selectedIDs; 
	std::set<long int>::iterator iter = IDs.begin();

	while( iter != IDs.end())
	{
		std::set<long int>::iterator iterRoot = this->rootID.begin();
		std::set<long int>::iterator iterFind = this->rootID.find(*iter);
		long int count = 0;

		while( iterRoot != iterFind && iterRoot != this->rootID.end())
		{
			iterRoot++;
			count++;
		}
		selectedIDs.insert( count);

		iter++;
	}

	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
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
}

void GraphWindow::UpdateGraphView()
{
	std::set<long int> IDs = this->selection->getSelections();

	if( this->view->GetRepresentation())
	{
		vtkAnnotationLink* annotationLink = this->view->GetRepresentation()->GetAnnotationLink();

		if( this->observerTag != 0)
		{
			annotationLink->RemoveObserver( this->observerTag);
			this->observerTag = 0;
		}

		vtkSelection* selection = annotationLink->GetCurrentSelection();
	
		vtkSmartPointer<vtkIdTypeArray> vertexList = vtkSmartPointer<vtkIdTypeArray>::New();

		std::set<long int>::iterator iter = IDs.begin();
		for( int id = 0; id < IDs.size(); id++, iter++)
		{
			vertexList->InsertNextValue( *iter);
		}

		vtkSmartPointer<vtkSelectionNode> selectNodeList = vtkSmartPointer<vtkSelectionNode>::New();
		selectNodeList->SetSelectionList( vertexList);
		selectNodeList->SetFieldType( vtkSelectionNode::VERTEX);
		selectNodeList->SetContentType( vtkSelectionNode::INDICES);

		selection->RemoveAllNodes();
		selection->AddNode(selectNodeList);
		annotationLink->SetCurrentSelection( selection);
		UpdataLookupTable( this->selection->getSelections());

		this->mainQTRenderWidget.GetRenderWindow()->Render();
		this->view->GetRenderer()->Render();
		this->observerTag = annotationLink->AddObserver("AnnotationChangedEvent", this->selectionCallback);
	}
}
