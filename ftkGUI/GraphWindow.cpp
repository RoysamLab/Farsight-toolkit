#include "GraphWindow.h"
#include <vtkAnnotationLink.h>
#include <vtkCommand.h>
#include <vtkDataRepresentation.h>
#include <vtkSelectionNode.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkSelection.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkViewTheme.h>
#include <vtkDataSetAttributes.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkAbstractArray.h>
#include <vtkVariantArray.h>
#include <QMessageBox>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <math.h>
#include <iomanip>

#define pi 3.1415926
#define ANGLE_STEP 90
#define MAX_HOP 200

double selectColor[3]={0,1,0};

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

void GraphWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, std::vector<int> *indexCluster, ObjectSelection * sels2)
{
	this->dataTable = table;
	this->indMapFromVertexToInd.clear();
	this->indMapFromIndToVertex.clear();

	for( long int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		long int var = this->dataTable->GetValue( i, 0).ToLong();
		this->indMapFromVertexToInd.insert( std::pair< long int, long int>(var, i));
		this->indMapFromIndToVertex.push_back( var);
	}

	this->indMapFromVertexToClusInd.clear();
	for( int i = 0; i < this->indMapFromClusIndToVertex.size(); i++)
	{
		this->indMapFromClusIndToVertex[i].clear();
		this->indMapFromClusIndToInd[i].clear();
	}
	this->indMapFromClusIndToVertex.clear();
	this->indMapFromClusIndToInd.clear();

	if( indexCluster != NULL)
	{
		int clusterSize = 0;
		for( int i = 0; i < indexCluster->size(); i++)
		{
			if( (*indexCluster)[i] + 1 > clusterSize)
			{
				clusterSize = (*indexCluster)[i] + 1;
			}
			this->indMapFromVertexToClusInd.insert( std::pair< int, int>(this->indMapFromIndToVertex[i], (*indexCluster)[i]));
		}
			/// rebuild the datamatrix for MST 
		std::vector<int> clus;
		for( int i = 0; i < clusterSize; i++)
		{
			this->indMapFromClusIndToVertex.push_back(clus);
			this->indMapFromClusIndToInd.push_back(clus);
		}
		for( int i = 0; i < indexCluster->size(); i++)
		{
			int index = (*indexCluster)[i];
			this->indMapFromClusIndToVertex[index].push_back( this->indMapFromIndToVertex[i]);
			this->indMapFromClusIndToInd[index].push_back( i);
		}
	}

	if(!sels)
		this->selection = new ObjectSelection();
	else
		this->selection = sels;

	if(!sels2)
		this->selection2 = new ObjectSelection();
	else
		this->selection2 = sels2;
	connect( this->selection, SIGNAL( changed()), this, SLOT( UpdateGraphView()));
}
	

void GraphWindow::SetGraphTable(vtkSmartPointer<vtkTable> table)
{
	//graphTable->Dump(8);	//debug dump
	this->TTG->ClearLinkVertices();
	this->TTG->SetInput(0, table);
	this->TTG->AddLinkEdge("Source", "Target"); 
	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New();
	
	theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	theme->SetLineWidth(5);
	theme->SetCellOpacity(0.9);
	theme->SetCellAlphaRange(0.5,0.5);
	theme->SetPointSize(10);
	theme->SetSelectedCellColor(1,0,1);
	theme->SetSelectedPointColor(1,0,1); 

	
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
	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New();
	
	theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	theme->SetLineWidth(5);
	theme->SetCellOpacity(0.9);
	theme->SetCellAlphaRange(0.5,0.5);
	theme->SetPointSize(10);
	theme->SetSelectedCellColor(1,0,1);
	theme->SetSelectedPointColor(1,0,1); 

	this->view->AddRepresentationFromInputConnection(TTG->GetOutputPort());
	/*this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName("Distance");*/
	this->view->SetLayoutStrategyToForceDirected();
	this->view->SetVertexLabelArrayName("label");
	this->view->VertexLabelVisibilityOn();
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::string xCol, std::string yCol, std::string zCol)
{
	std::cout<< "SetGraphTable"<<endl;

	vtkAbstractArray *arrayID1 = table->GetColumnByName( ID1.c_str());
	vtkAbstractArray *arrayID2 = table->GetColumnByName( ID2.c_str());
	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkIntArray> vertexIDarrays = vtkSmartPointer<vtkIntArray>::New();
	vertexIDarrays->SetNumberOfComponents(1);
	vertexIDarrays->SetName("vertexIDarrays");

  // Create the edge weight array
	vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
	weights->SetNumberOfComponents(1);
	weights->SetName("edgeLabel");

	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkMutableUndirectedGraph::New();
	for( int i = 0; i <  this->dataTable->GetNumberOfRows(); i++)
	{
		int vertexID = graph->AddVertex();
		//add coord points for pass through here
		double x = this->dataTable->GetValueByName(i, xCol.c_str()).ToDouble();
		double y = this->dataTable->GetValueByName(i, yCol.c_str()).ToDouble();
		double z = this->dataTable->GetValueByName(i, zCol.c_str()).ToDouble();
		points->InsertNextPoint(x,y,z);
		vertexIDarrays->InsertNextValue( this->indMapFromIndToVertex[i]);
	}

	for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
	{
		long int ver1 = arrayID1->GetVariantValue(i).ToLong();
		long int ver2 = arrayID2->GetVariantValue(i).ToLong();
		std::map< long int, long int>::iterator iter1 = this->indMapFromVertexToInd.find( ver1);
		std::map< long int, long int>::iterator iter2 = this->indMapFromVertexToInd.find( ver2);
		if( iter1 != this->indMapFromVertexToInd.end() && iter2 != this->indMapFromVertexToInd.end())
		{
			long int index1 = iter1->second;
			long int index2 = iter2->second;
			graph->AddEdge( index1, index2);
			weights->InsertNextValue(table->GetValueByName(i, edgeLabel.c_str()).ToDouble());
		}
		else
		{
			QMessageBox msg;
			msg.setText("Index Mapping Error!");
			msg.exec();
			exit(-1);
		}
	}
	
	theme = vtkSmartPointer<vtkViewTheme>::New();
	theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	theme->SetLineWidth(5);
	theme->SetCellOpacity(0.9);
	theme->SetCellAlphaRange(0.8,0.8);
	theme->SetPointSize(8);
	theme->SetSelectedCellColor(1,0,0);
	theme->SetSelectedPointColor(1,0,0); 

	vtkSmartPointer<vtkIntArray> vertexColors = vtkSmartPointer<vtkIntArray>::New();
	vertexColors->SetNumberOfComponents(1);
	vertexColors->SetName("Color");
	
	this->lookupTable->SetNumberOfTableValues( this->dataTable->GetNumberOfRows());

	//std::cerr << "Number of Lookup Table values: " << this->dataTable->GetNumberOfRows() << std::endl;

	for( vtkIdType i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		vertexColors->InsertNextValue( i);
		this->lookupTable->SetTableValue(i, 0, 0, 1); // color the vertices- blue
    }
	lookupTable->Build();

	graph->GetVertexData()->AddArray(vertexColors);
	graph->GetVertexData()->AddArray(vertexIDarrays);
	graph->GetEdgeData()->AddArray(weights);
	graph->SetPoints(points);

	this->view->AddRepresentationFromInput( graph);
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetColorVertices(true); 
	this->view->SetVertexLabelVisibility(true);

	this->view->SetVertexColorArrayName("Color");

    theme->SetPointLookupTable(lookupTable);
    theme->SetBackgroundColor(0,0,0); 
	this->view->ApplyViewTheme(theme);

	this->view->SetEdgeLabelArrayName("edgeLabel");

	//this->view->SetLayoutStrategyToForceDirected();
	this->view->SetLayoutStrategyToPassThrough();

	this->view->SetVertexLabelArrayName("vertexIDarrays");
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::set<long int>* colSels, QString filename)
{
	this->fileName = filename;
	vtkAbstractArray *arrayID1 = table->GetColumnByName( ID1.c_str());
	vtkAbstractArray *arrayID2 = table->GetColumnByName( ID2.c_str());
	if( colSels)
	{
		this->colSelectIDs = *colSels;
	}

	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	

	vtkSmartPointer<vtkIntArray> vertexIDarrays = vtkSmartPointer<vtkIntArray>::New();
	vertexIDarrays->SetNumberOfComponents(1);
	vertexIDarrays->SetName("vertexIDarrays");

  // Create the edge weight array
	vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
	weights->SetNumberOfComponents(1);
	weights->SetName("edgeLabel");

	std::cout<< "construct graph"<<endl;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkMutableUndirectedGraph::New();
	vnl_matrix<long int> adj_matrix( table->GetNumberOfRows() + 1, table->GetNumberOfRows() + 1);
	
	vertextList.set_size( table->GetNumberOfRows(), 3);
	edgeWeights.set_size( table->GetNumberOfRows());

	for( int i = 0; i < table->GetNumberOfRows() + 1; i++)
	{
		int vertexID = graph->AddVertex();
		if (this->indMapFromClusIndToVertex.size() <= 0)
		{
			vertexIDarrays->InsertNextValue( this->indMapFromIndToVertex[i]);
		}
		else
		{
			vertexIDarrays->InsertNextValue( this->indMapFromClusIndToVertex[i].size());   // which should be the cluster index
		}
	}

	for( int i = 0; i < table->GetNumberOfRows(); i++)
	{

		long int ver1 = arrayID1->GetVariantValue(i).ToLong();
		long int ver2 = arrayID2->GetVariantValue(i).ToLong();

		if (this->indMapFromClusIndToVertex.size() <= 0)
		{
			std::map< long int, long int>::iterator iter1 = this->indMapFromVertexToInd.find( ver1);
			std::map< long int, long int>::iterator iter2 = this->indMapFromVertexToInd.find( ver2);
			if( iter1 != this->indMapFromVertexToInd.end() && iter2 != this->indMapFromVertexToInd.end())
			{
				long int index1 = iter1->second;
				long int index2 = iter2->second;
				graph->AddEdge( index1, index2);
				adj_matrix( index1, index2) = 1;
				adj_matrix( index2, index1) = 1;
				vertextList( i, 0) = index1;
				vertextList( i, 1) = index2;
				edgeWeights[i] = table->GetValueByName(vtkIdType(i), edgeLabel.c_str()).ToDouble();
				weights->InsertNextValue(table->GetValueByName(vtkIdType(i), edgeLabel.c_str()).ToDouble());
			}
			else
			{
				QMessageBox msg;
				msg.setText("Index Mapping Error!");
				msg.exec();
				exit(-1);
			}
		}
		else    // cluster index, no mapping
		{
			graph->AddEdge( ver1, ver2);
			adj_matrix( ver1, ver2) = 1;
			adj_matrix( ver2, ver1) = 1;
			vertextList( i, 0) = ver1;
			vertextList( i, 1) = ver2;
			edgeWeights[i] = table->GetValueByName(vtkIdType(i), edgeLabel.c_str()).ToDouble();
			weights->InsertNextValue(table->GetValueByName(vtkIdType(i), edgeLabel.c_str()).ToDouble());
		}
	}

	//edgeWeights = edgeWeights - edgeWeights.mean();
	//double std = edgeWeights.two_norm();
	//if( std != 0)
	//{
	//	edgeWeights = edgeWeights / std;
	//	edgeWeights = ( edgeWeights - edgeWeights.min_value()) / ( edgeWeights.max_value() - edgeWeights.min_value()) + 0.1;
	//	edgeWeights = 2.0 * edgeWeights;
	//}
	edgeWeights = edgeWeights / edgeWeights.max_value() * 5;
	for(int i = 0; i < vertextList.rows(); i++)
	{
		vertextList( i, 2) = edgeWeights[i];
	}

	std::ofstream verofs("vertextList.txt");
	verofs<< vertextList<<endl;
	verofs.close();

	std::cout<< "calculate coordinates"<<endl;
	std::vector<Point> oldPointList;
	std::vector<Point> newPointList;
	CalculateCoordinates(adj_matrix, oldPointList);
	newPointList = oldPointList;
	UpdateCoordinatesByEdgeWeights( oldPointList, vertextList, newPointList);

	//std::vector<Point> newPointList;
	//CalculateCoordinates(adj_matrix, newPointList);

	if( newPointList.size() > 0)
	{
		std::cout<< newPointList.size()<<endl;
		for( int i = 0; i <  newPointList.size(); i++)
		{
			points->InsertNextPoint(newPointList[i].x, newPointList[i].y, 0);
		}
		graph->SetPoints( points);
	}
	else
	{	
		QMessageBox msg;
		msg.setText("No coordinates Input!");
		msg.exec();
		exit(-2);
	}

	vtkSmartPointer<vtkIntArray> vertexColors = vtkSmartPointer<vtkIntArray>::New();
	vertexColors->SetNumberOfComponents(1);
	vertexColors->SetName("Color");
	
	this->lookupTable->SetNumberOfTableValues( table->GetNumberOfRows() + 1);


	for( vtkIdType i = 0; i < table->GetNumberOfRows() + 1; i++)               
	{             
		vertexColors->InsertNextValue( i);
		this->lookupTable->SetTableValue(i, 0, 0, 1); // color the vertices- blue
	}

	lookupTable->Build();

	graph->GetVertexData()->AddArray(vertexColors);
	graph->GetVertexData()->AddArray(vertexIDarrays);
	graph->GetEdgeData()->AddArray(weights);
	
	//theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	theme->SetLineWidth(3);
	theme->SetCellOpacity(0.9);
	theme->SetCellAlphaRange(0.8,0.8);
	theme->SetPointSize(5);
	theme->SetCellColor(0.6,0.6,0.6);
	theme->SetSelectedCellColor(selectColor);
	theme->SetSelectedPointColor(selectColor); 
	theme->SetVertexLabelColor(0.3,0.3,0.3);
	theme->SetBackgroundColor(1,1,1); 
	theme->SetBackgroundColor2(1,1,1);

	this->view->RemoveAllRepresentations();
	this->view->SetRepresentationFromInput( graph);
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetColorVertices(true); 
	this->view->SetVertexLabelVisibility(true);

	this->view->SetVertexColorArrayName("Color");
	this->view->SetEdgeLabelArrayName("edgeLabel");
	this->view->SetVertexLabelArrayName("vertexIDarrays");

    theme->SetPointLookupTable(lookupTable);
	this->view->ApplyViewTheme(theme);

	this->view->SetLayoutStrategyToPassThrough();
	this->view->SetVertexLabelFontSize(15);
}

void GraphWindow::UpdateCoordinatesByEdgeWeights(std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList)
{
	long int oldFirstIndex = backbones[0];
	long int oldSecondIndex = 0;
	newPointList[ oldFirstIndex] = oldPointList[oldFirstIndex];
	Point newFirst = oldPointList[oldFirstIndex];
	Point oldFirst = oldPointList[oldFirstIndex];
	Point oldSecond(0,0);
	for( int i = 1; i < backbones.size(); i++)
	{
		oldSecondIndex = backbones[i];
		oldSecond = oldPointList[oldSecondIndex];
		double weight = GetEdgeWeight( vertexList, oldFirstIndex, oldSecondIndex);
		Point newSecond = GetNewPointFromOldPoint(oldFirst, oldSecond, newFirst, weight);
		newPointList[ oldSecondIndex] = newSecond;
		newFirst = newSecond;
		oldFirst = oldSecond;
		oldFirstIndex = oldSecondIndex;
		UpdateChainPointList(oldSecondIndex, oldPointList, vertexList, newPointList);
	}
}

// attachnode's coordinate has been updated already in pointlist
void GraphWindow::UpdateChainPointList(long int attachnode, std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList)
{
	for( long int i = 0; i < chainList.size(); i++)
	{
		std::pair<long int, std::vector<long int> > pair = chainList[i];
		if( pair.first == attachnode)
		{
 			Point oldFirst = oldPointList[ attachnode];
			Point newFirst = newPointList[ attachnode];
			std::vector<long int> vec = pair.second;	
			Point oldSecond(0,0);
			long int oldFirstIndex = attachnode; 
			long int oldSecondIndex = 0;
			for( long int j = 0; j < vec.size(); j++)
			{
				oldSecondIndex = vec[j];
				oldSecond = oldPointList[ oldSecondIndex];
				double weight = GetEdgeWeight( vertexList, oldFirstIndex, oldSecondIndex);
				Point newSecond = GetNewPointFromOldPoint(oldFirst, oldSecond, newFirst, weight);
				newPointList[oldSecondIndex] = newSecond;
				newFirst = newSecond;
				oldFirst = oldSecond;
				oldFirstIndex = oldSecondIndex;	

				UpdateChainPointList(oldSecondIndex, oldPointList, vertexList, newPointList);
			}
		}
	}
}

double GraphWindow::GetEdgeWeight( vnl_matrix<double>& vertexList, long firstIndex, long secondIndex)
{
	double weight;
	for( int i = 0; i < vertexList.rows(); i++)
	{
		if( vertexList(i, 0) == firstIndex && vertexList( i, 1)== secondIndex
			|| vertexList(i, 0) ==  secondIndex && vertexList(i, 1) ==  firstIndex)
		{
			weight = vertexList(i,2);
			break;
		}
	}
	return weight;
}

/// calculate the new point nodes based on the two old points and their weight
Point GraphWindow::GetNewPointFromOldPoint( Point &oldPointFirst, Point &oldPointSecond, Point &newPointFirst, double weight)
{
	Point newpoint(0,0);
	double xd = (oldPointFirst.x - oldPointSecond.x) * (oldPointFirst.x - oldPointSecond.x);
	double yd = (oldPointFirst.y - oldPointSecond.y) * (oldPointFirst.y - oldPointSecond.y);
	double distance = sqrt( xd + yd);

	newpoint.x = newPointFirst.x + ( oldPointSecond.x - oldPointFirst.x) * weight / distance;
	newpoint.y = newPointFirst.y + ( oldPointSecond.y - oldPointFirst.y) * weight /  distance;
	return newpoint;
}

void GraphWindow::ShowGraphWindow()
{
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
	std::cout<< "view->ResetCamera"<<endl;
	view->ResetCamera();
	std::cout<< "view->Render"<<endl;
	view->Render();
	std::cout<< "Successfully rendered"<<endl;

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
		graphWin->SetSelectedIds2();  // only select the selected feature columns
	}
	//graphWin->mainQTRenderWidget.GetRenderWindow()->Render();
	//graphWin->view->GetRenderer()->Render();
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
		std::set<long int>::iterator iter = IDs.begin();
		while( iter != IDs.end())
		{
			if( this->indMapFromClusIndToVertex.size() <= 0)
			{
				long int var = indMapFromIndToVertex[*iter];
				selectedIDs.insert( var);
			}
			else
			{
				std::vector<int> clusVertex = indMapFromClusIndToVertex[*iter];
				for( int ind = 0; ind < clusVertex.size(); ind++)
				{
					selectedIDs.insert(clusVertex[ind]);
				}
			}
			iter++;
		}
		this->selection->select( selectedIDs);
		UpdataLookupTable( selectedIDs);
	}
	this->view->GetRenderer()->Render();
}

void GraphWindow::SetSelectedIds2()
{
	if( this->colSelectIDs.size() > 0)
	{
		this->selection2->select( this->colSelectIDs);
	}
}

void GraphWindow::UpdataLookupTable( std::set<long int>& IDs)
{
	std::set<long int> selectedIDs; 
	std::set<long int>::iterator iter = IDs.begin();

	while( iter != IDs.end())
	{
		long int var = 0;
		if( this->indMapFromClusIndToVertex.size() <= 0)
		{
			var = this->indMapFromVertexToInd.find( *iter)->second;
		}
		else
		{
			var = this->indMapFromVertexToClusInd.find( *iter)->second;
		}

		if( selectedIDs.find( var) == selectedIDs.end())   // selectIDs doesn't have var
		{
			selectedIDs.insert( var);
		}
		iter++;
	}
	//std::cerr << "Number of Lookup Table values: " << this->dataTable->GetNumberOfRows() << std::endl;
	int vertexnum = 0;
	if( this->indMapFromClusIndToVertex.size() <= 0)
	{
		vertexnum = this->dataTable->GetNumberOfRows();
	}
	else
	{
		vertexnum = this->indMapFromClusIndToVertex.size();  // vertex size equals cluster size
	}

	for( vtkIdType i = 0; i < vertexnum; i++)
	{
		if (selectedIDs.find(i) != selectedIDs.end())
		{
			this->lookupTable->SetTableValue(i, selectColor[0], selectColor[1], selectColor[2]); // color the vertices
		}
		else
		{
			this->lookupTable->SetTableValue(i, 0, 0, 1); // color the vertices- blue
		}
	}

	this->lookupTable->Build();
}

void GraphWindow::UpdateGraphView()
{
	std::set<long int> IDs;
	std::set<long int> IDsnew = this->selection->getSelections();
	if( this->indMapFromClusIndToVertex.size() <= 0)
	{
		IDs = IDsnew;
	}
	else
	{
		std::set<long int>::iterator iter = IDsnew.begin();
		while( iter != IDsnew.end())
		{
			int var = this->indMapFromVertexToClusInd.find( *iter)->second;

			if( IDs.find( var) == IDs.end())   // selectIDs doesn't have var
			{
				IDs.insert( var);
			}
			iter++;
		}
	}

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
		vertexList->SetNumberOfComponents(1);

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

		if( this->colSelectIDs.size() > 0)
		{
			selection2->select( this->colSelectIDs);    // only select the selected feature columns
		}

		std::set<long int> updataIDes =  this->selection->getSelections();
		UpdataLookupTable( updataIDes);

		this->mainQTRenderWidget.GetRenderWindow()->Render();
		this->view->GetRenderer()->Render();
		this->observerTag = annotationLink->AddObserver("AnnotationChangedEvent", this->selectionCallback);
	}
}

void GraphWindow::CalculateCoordinates(vnl_matrix<long int>& adj_matrix, std::vector<Point>& pointList)
{
	vnl_matrix<long int> shortest_hop(adj_matrix.rows(), adj_matrix.cols());
	vnl_vector<long int> mark(adj_matrix.rows());
	shortest_hop.fill(0);
	mark.fill(0);
	mark[0] = 1;
	std::vector< long int> checkNode;
	std::vector<long int> noncheckNode;
	find( mark, 0, noncheckNode, checkNode);

	/// calculate the shortest hop matrix
	while( noncheckNode.size() > 0)
	{
		long int checkNodeInd = -1;
		long int noncheckNodeInd = -1;
		for( long int i = 0; i < checkNode.size(); i++)
		{
			for( long int j = 0; j < noncheckNode.size(); j++)
			{
				if( adj_matrix(checkNode[i], noncheckNode[j]) != 0)
				{
					checkNodeInd = checkNode[i];
					noncheckNodeInd =  noncheckNode[j];
					break;
				}
			}
			if( checkNodeInd != -1 && noncheckNodeInd != -1)
			{
				break;
			}
		}

		for( long int i = 0; i < checkNode.size(); i++)
		{
			if( indMapFromClusIndToInd.size() <= 0)
			{
				shortest_hop( checkNode[i], noncheckNodeInd) =  shortest_hop( checkNode[i], checkNodeInd) + 1;
				shortest_hop( noncheckNodeInd, checkNode[i]) =  shortest_hop( checkNode[i], checkNodeInd) + 1;
			}
			else
			{
				int size = ( this->indMapFromClusIndToInd[noncheckNodeInd].size() + this->indMapFromClusIndToInd[checkNodeInd].size()) / 2;
				shortest_hop( checkNode[i], noncheckNodeInd) =  shortest_hop( checkNode[i], checkNodeInd) + size;
				shortest_hop( noncheckNodeInd, checkNode[i]) =  shortest_hop( checkNode[i], checkNodeInd) + size;
			}
		}

		mark[ noncheckNodeInd] = 1;

		checkNode.clear();
		noncheckNode.clear();
		find( mark, 0, noncheckNode, checkNode);
	}

	QString hopStr = this->fileName + "shortest_hop.txt";
	ofstream ofs(hopStr.toStdString().c_str());
	ofs<< shortest_hop<<endl;
	ofs.close();

	/// find the root and chains of the tree
	long int maxhop = shortest_hop.max_value();
	unsigned int maxId = shortest_hop.arg_max();
	unsigned int coln = maxId / shortest_hop.rows();
	unsigned int rown = maxId  - coln * shortest_hop.rows();

	backbones.clear();
	chainList.clear();  // store the chains, no repeat

	std::queue<long int> tmpbackbones;   // for calculating all the sidechains
	std::vector<long int> debugbackbones;  // for debugging
	std::map< long int, long int> tmpChain;  // for ordering 
	vnl_vector< int> tag( shortest_hop.rows());     // whether the node has been included in the chain
	tag.fill( 0);

	std::cout<<"build chains"<<endl;
	QString chainStr = this->fileName + "chains.txt";
	ofstream ofChains( chainStr.toStdString().c_str());
	/// find the backbone
	for( long int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( coln, i) + shortest_hop( rown, i) == maxhop)
		{
			tmpChain.insert( std::pair< long int, long int>( shortest_hop( coln, i), i));
		}
	}

	std::map< long int, long int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		backbones.push_back((*iter).second);
		tmpbackbones.push((*iter).second);
		debugbackbones.push_back( (*iter).second);
		tag[ (*iter).second] = 1;
	}

	/// find the branches' backbones
	std::vector< long int> branchnodes;
	std::vector< long int> chains;
	while( !tmpbackbones.empty())
	{
		long int ind = tmpbackbones.front(); 
		for( long int i = 0; i < adj_matrix.cols(); i++)
		{
			if( adj_matrix( ind, i) != 0 && tag[i] == 0)    // ind's  neighbours that haven't been checked
			{
				branchnodes.clear();
				branchnodes.push_back( ind);
				for( long int j = 0; j < shortest_hop.cols(); j++)
				{
					if( shortest_hop( ind, j) > shortest_hop( i, j))  // find neighbour i's branch nodes including i
					{
						if( tag[j] == 0)
						{
							branchnodes.push_back( j);
						}
						else
						{
							ofChains<< "already searched branchnode: "<< ind<<"\t"<< i <<"\t"<< j<<endl;
							for( long int st = 0; st < debugbackbones.size(); st++)
							{
								ofChains<< debugbackbones[st] + 1<<endl;
							}
							ofChains<< tag<<endl;
						}
					}
				}

				if( branchnodes.size() > 1)
				{
					getBackBones( shortest_hop, branchnodes, chains);
					chainList.push_back( std::pair< long int, std::vector<long int> >(ind, chains));
					for( long int k = 0; k < chains.size(); k++)
					{
						tmpbackbones.push( chains[k]);
						debugbackbones.push_back( chains[k]);
						tag[ chains[k]] = 1;
					}
				}
			}
		}
		tmpbackbones.pop();
	}

	ofChains << "backbones:" <<endl;
	for( long int i = 0; i < backbones.size(); i++)
	{
		ofChains << backbones[i]<<"\t";
	}
	ofChains <<endl;
	ofChains << "branch chains:"<<endl;
	std::vector< std::pair< long int, std::vector<long int> > >::iterator chainIter;
	for( chainIter = chainList.begin(); chainIter != chainList.end(); chainIter++)
	{
		std::pair< long int, std::vector< long int> >tmp = *chainIter;
		ofChains << tmp.first;
		std::vector< long int> branchlist = tmp.second;
		for( long int i = 0; i < branchlist.size(); i++)
		{
			ofChains << "\t"<< branchlist[i];
		}
		ofChains <<endl;
	}

	SortChainList( shortest_hop, backbones, chainList);   // the order to draw the subbones
	//for( chainIter = chainList.begin(); chainIter != chainList.end(); chainIter++)
	//{
	//	std::pair< long int, std::vector< long int> >tmp = *chainIter;
	//	ofChains << tmp.first;
	//	std::vector< long int> branchlist = tmp.second;
	//	for( long int i = 0; i < branchlist.size(); i++)
	//	{
	//		ofChains << "\t"<< branchlist[i];
	//	}
	//	ofChains <<endl;
	//}
	ofChains.close();

	/// calculate the coordinates of the nodes
	std::cout<<"calculate nodes position"<<endl;
	QString corStr = this->fileName + "coordinates.txt";
	ofstream ofCoordinate( corStr.toStdString().c_str());
	ofCoordinate.precision(4);

	vnl_matrix< double> nodePos( 2, adj_matrix.cols());
	vnl_vector< int> nodePosAssigned( adj_matrix.cols());
	nodePos.fill(0);
	nodePosAssigned.fill(0);

	/// calculate backbone node position
	vnl_vector< double> backboneAngle( backbones.size());
	for( long int i = 0; i < backboneAngle.size(); i++)
	{
		backboneAngle[i] = i + 1;
	}
	backboneAngle = backboneAngle -  backboneAngle.mean();
	backboneAngle = backboneAngle / backboneAngle.max_value() * pi / ANGLE_STEP * 25;
	
	double powx = pow( sin(backboneAngle[0]) - sin(backboneAngle[1]), 2);
	double powy = pow( cos(backboneAngle[1]) - cos(backboneAngle[0]), 2);
	double norm = sqrt( powx + powy);

	if( backbones.size() > 500)
	{
		long int size = backbones.size();
		for( long int i = 0; i < backbones.size(); i++)
		{
			nodePos(0, backbones[i]) = sin( backboneAngle[i]) / norm * 500 / size;
			nodePos(1, backbones[i]) = -cos( backboneAngle[i]) / norm * 500 / size;
			nodePosAssigned[ backbones[i]] = 1;
		}
	}
	else
	{
		for( long int i = 0; i < backbones.size(); i++)
		{
			nodePos(0, backbones[i]) = sin( backboneAngle[i]) / norm;
			nodePos(1, backbones[i]) = -cos( backboneAngle[i]) / norm;
			nodePosAssigned[ backbones[i]] = 1;
		}
	}

	/// cacluate the branch backbone nodes position
	for( long int i = 0; i < chainList.size(); i++)
	{
		if( i % 100 == 0)
		{
			std::cout<<i<<endl;
		}
		std::pair< long int, std::vector<long int> > branch = chainList[i];
		long int attachNode = branch.first;
		std::vector< long int> branchNode = branch.second;
		for( long int j = 0; j < branchNode.size(); j++)
		{
			long int newNode = branchNode[j];
			double bestR = 0.3;
			double bestTheta = 0;
			vnl_vector< double> minForce( 7);
			vnl_vector< double> mintheta( 7);
			minForce.fill(0);
			mintheta.fill(0);

			int count = 0;
			for( double r = 0.3; r <= 0.9; r = r + 0.1, count++)
			{
				bool bfirst = true;
				double minVar = 0;
				for( double theta = 0; theta <= 2 * pi; theta += 2 * pi / ANGLE_STEP)
				{
					Point newNodePoint;
					vnl_matrix<double> repel_mat;
					GetElementsIndexInMatrix( shortest_hop, newNode, MAX_HOP, nodePos, repel_mat, nodePosAssigned);

					newNodePoint.x = nodePos( 0, attachNode) + r * cos( theta);
					newNodePoint.y = nodePos( 1, attachNode) + r * sin( theta);
					for( long int k = 0; k < repel_mat.cols(); k++)
					{
						repel_mat(0, k) = newNodePoint.x - repel_mat( 0, k);
						repel_mat(1, k) = newNodePoint.y - repel_mat( 1, k);
					}

					vnl_vector< double> repel_tmp_xvector(repel_mat.cols());
					vnl_vector< double> repel_tmp_yvector(repel_mat.cols());
					for( long int k = 0; k < repel_mat.cols(); k++)
					{
						double force = sqrt( pow(repel_mat(0, k),2) + pow(repel_mat(1, k),2));
						force = pow( force, 5);
						repel_tmp_xvector[ k] = repel_mat(0, k) / force;
						repel_tmp_yvector[ k] = repel_mat(1, k) / force;
					}

					vnl_vector<double> repel_force(2);
					repel_force[ 0] = repel_tmp_xvector.sum();
					repel_force[ 1] = repel_tmp_yvector.sum();

					double cosalfa = (repel_force[0] * cos( theta) + repel_force[1] * sin( theta)) / repel_force.magnitude();
					double sinalfa = sqrt( 1 - pow( cosalfa, 2));
		
					if( bfirst)
					{
						if( repel_force.magnitude() * cosalfa >= 0)
						{
							minForce[ count] = repel_force.magnitude() * cosalfa;
							mintheta[ count] = theta;
							minVar = repel_force.magnitude() * sinalfa;
							bfirst = false;
						}
					}
					else
					{
						if( repel_force.magnitude() * cosalfa >= 0 && repel_force.magnitude() * sinalfa < minVar )
						{
							minForce[ count] = repel_force.magnitude() * cosalfa;
							mintheta[ count] = theta;
							minVar = repel_force.magnitude() * sinalfa;
						}
					}
				}
			}
			
			int min = minForce.arg_min();
			bestR = 0.3 + 0.1 * min;
			bestTheta = mintheta[ min];

			nodePos(0, newNode) = nodePos(0, attachNode) + bestR * cos( bestTheta);
			nodePos(1, newNode) = nodePos(1, attachNode) + bestR * sin( bestTheta);
			nodePosAssigned[ newNode] = 1;
			attachNode = newNode;
		}
	}

	double medianX = Median( nodePos.get_row(0));
	double medianY = Median( nodePos.get_row(1));
	double medianXY[2] = {medianX, medianY};

	for( int i = 0; i < 2; i++)
	{
		vnl_vector<double> tmpNodePos = nodePos.get_row(i) - medianXY[i];
		vnl_vector<double> tmp = tmpNodePos;
		for( long int j = 0; j < tmpNodePos.size(); j++)
		{
			tmp[j] = abs( tmp[j]);
		}
		
		tmpNodePos = tmpNodePos / tmp.max_value() * 50;
		nodePos.set_row(i, tmpNodePos);
	}

	ofCoordinate << setiosflags(ios::fixed)<< nodePos.transpose()<<endl;
	ofCoordinate.close();

	for( long int i = 0; i < nodePos.cols(); i++)
	{
		Point pt( nodePos(0, i), nodePos(1, i));
		pointList.push_back( pt);
	}
}

/// branchnodes first element is the chain's attach node
void GraphWindow::getBackBones(vnl_matrix< long int>& shortest_hop, std::vector< long int>& branchnodes, std::vector< long int>& chains)
{
	chains.clear();
	vnl_vector< long int> branchShortestHop( branchnodes.size());

	long int attachNode = branchnodes.front();
	std::map< long int, long int> tmpChain;  // for ordering 

	for( long int i = 0; i < branchShortestHop.size(); i++)
	{
		branchShortestHop[i] = shortest_hop( attachNode, branchnodes[i]);
	}

	long int maxHop = branchShortestHop.max_value();
	long int endNodeIndex = branchShortestHop.arg_max();
	long int endNode = branchnodes[ endNodeIndex];
	for( long int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( attachNode, i) + shortest_hop( endNode, i) == maxHop)
		{
			tmpChain.insert( std::pair< long int, long int>( shortest_hop( attachNode, i), i));
		}
	}

	std::map< long int, long int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		if( (*iter).second != attachNode)
		{
			chains.push_back((*iter).second);
		}
	}
}

void GraphWindow::find(vnl_vector<long int>& vec, long int val, std::vector<long int>& equal, std::vector<long int>& nonequal)
{
	for( unsigned int i = 0; i < vec.size(); i++)
	{
		if( vec[i] == val)
		{
			equal.push_back(i);
		}
		else
		{
			nonequal.push_back(i);
		}
	}
}

void GraphWindow::GetElementsIndexInMatrix(vnl_matrix<long int>& mat, long int rownum, long int max, vnl_matrix<double>& oldmat, vnl_matrix<double>& newmat, vnl_vector< int>& tag)
{
	std::vector< long int> index;
	for( long int i = 0; i < mat.cols(); i++)
	{
		if( mat( rownum, i) < max && tag[i] == 1)
		{
			index.push_back( i);
		}
	}

	newmat.set_size( oldmat.rows(), index.size());
	for( long int i = 0; i < index.size(); i++)
	{
		vnl_vector< double> tmpcol = oldmat.get_column( index[i]);
		newmat.set_column( i, tmpcol);
	}
}

double GraphWindow::Median( vnl_vector<double> vec)
{
	vnl_vector<double> vect = vec;
	vnl_vector<double> tmp( vec.size());
	double max = vect.max_value() + 1;
	double med;
	for( long int i = 0; i < vect.size(); i++)
	{
		tmp[i] = vect.min_value();
		long int ind = vect.arg_min();
		vect[ind] = max;
	}

	long int size = tmp.size();
	if( size % 2 == 1)
	{
		med = tmp[ size / 2];
	}
	else
	{
		med = 1.0 / 2.0 * ( tmp[ size / 2] + tmp[ size / 2 - 1]);
	}
	return med;
}

void GraphWindow::SortChainList( vnl_matrix<long int>& shortest_hop, std::vector<long int>& backbones, 
				   std::vector< std::pair<long int, std::vector<long int> > >& chainList)
{
	std::vector< std::pair<long int, std::vector<long int> > > tmpchainList;
	std::multimap< long int, long int> sortMap;
	long int startNode = backbones[0];
	long int endNode = backbones[ backbones.size() - 1];
	for( long int i = 0; i < chainList.size(); i++)
	{
		std::pair<long int, std::vector<long int> > chainPair = chainList[i];
		if( IsExist( backbones, chainPair.first))
		{
			long int abHop = abs( shortest_hop(chainPair.first, startNode) - shortest_hop( chainPair.first, endNode));
			sortMap.insert( std::pair< long int, long int>(abHop, i));
		}
		else
		{
			std::multimap< long int, long int>::iterator iter;
			for( iter = sortMap.begin(); iter != sortMap.end(); iter++)
			{
				std::pair<long int, std::vector<long int> > pair = chainList[ (*iter).second];
				tmpchainList.push_back( pair);
			}
			for( long int k = i; k < chainList.size(); k++)
			{
				tmpchainList.push_back( chainList[k]);
			}
			chainList = tmpchainList;
			break;
		}
	}
}

bool GraphWindow::IsExist(std::vector<long int>& vec, long int value)
{
	for( long int i = 0; i < vec.size(); i++)
	{
		if( value == vec[i])
		{
			return true;
		}
	}
	return false;
}

void GraphWindow::GetProgressionTreeOrder(std::vector<long int> &order)
{
	order.clear();
	std::vector<long int> clusterOrder;
	for( long int i = 0; i < backbones.size(); i++)
	{
		clusterOrder.push_back(backbones[i]);
		GetOrder(backbones[i], clusterOrder);
	}

	if( this->indMapFromClusIndToInd.size() > 0 )
	{
		for( int i = 0; i < clusterOrder.size(); i++)
		{
			std::vector<int> clusterNodes = this->indMapFromClusIndToInd[ clusterOrder[i]];
			for( int j = 0; j < clusterNodes.size(); j++)
			{
				order.push_back( clusterNodes[j]);
			}
		}
	}
	else
	{
		order = clusterOrder;
	}
}

void GraphWindow::GetOrder(long int node, std::vector<long int> &order)
{
	std::vector<long int> vec;
	for( long int i = 0; i < chainList.size(); i++)
	{
		std::pair<long int, std::vector<long int> > pair = chainList[i];
		if( pair.first == node)
		{
			vec = pair.second;
			for( long int j = 0; j < vec.size(); j++)
			{
				order.push_back( vec[j]);
				GetOrder( vec[j], order);
			}
		}
	}
	vec.clear();
}

void GraphWindow::closeEvent(QCloseEvent *event)
{
	mainQTRenderWidget.close();
}