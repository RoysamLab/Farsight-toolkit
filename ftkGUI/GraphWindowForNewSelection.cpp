#include "GraphWindowForNewSelection.h"
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
#include <vtkStringArray.h>
#include <vtkCellPicker.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <QMessageBox>
#include <QInputDialog>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <math.h>
#include <iomanip>

#define pi 3.1415926
#define ANGLE_STEP 90
#define MAX_HOP 200
#define ANNOTATION_CORNER 0

static double selectColor[3]={0,1,0};
static double progressionPathColor[3] = {1,0,0};
static double edgeDefaultColor[3] = {0.6, 0.6, 0.6};

GraphWindowForNewSelection::GraphWindowForNewSelection(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;// = new QVTKWidget;
	this->TTG = vtkSmartPointer<vtkTableToGraph>::New();
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();

	bProgressionStart = true;
	progressionStartID = -1;
	progressionEndID = -1;
	cornAnnotation = NULL;
}

GraphWindowForNewSelection::~GraphWindowForNewSelection()
{
}

void GraphWindowForNewSelection::setModels(vtkSmartPointer<vtkTable> table, SelectiveClustering * clusterSelection)
{
	this->dataTable = table;
	this->clusterIds = vtkSmartPointer<vtkIdTypeArray>::New();
	clusterIds->SetName("PedigreeIDs");
	vtkAbstractArray *idArray = table->GetColumn(0);
	indVectorFromIndToVertex.resize(idArray->GetNumberOfTuples());
	std::cout<< idArray->GetNumberOfTuples()<<endl;
	for( int i = 0; i < idArray->GetNumberOfTuples(); i++)
	{
		indVectorFromIndToVertex[i] = idArray->GetVariantValue(i).ToInt();
		this->clusterIds->InsertNextValue(indVectorFromIndToVertex[i]);
		indMapFromVertexToInd.insert( std::pair< int, int>(indVectorFromIndToVertex[i], i));
	}

	if( clusterSelection == NULL)
	{
		clusterSelection = new SelectiveClustering();
	}
	else
	{
		ClusterSelections = clusterSelection;
	}

}
	

void GraphWindowForNewSelection::SetGraphTable(vtkSmartPointer<vtkTable> table)
{
	//graphTable->Dump(8);	//debug dump
	this->TTG->ClearLinkVertices();
	this->TTG->SetInputData(0, table);
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

void GraphWindowForNewSelection::SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2)
{
	//graphTable->Dump(8);	//debug dump
	this->TTG->ClearLinkVertices();
	this->TTG->SetInputData(0, table);
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

void GraphWindowForNewSelection::SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::string xCol, std::string yCol, std::string zCol)
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
		vertexIDarrays->InsertNextValue( this->indVectorFromIndToVertex[i]);
	}

	for( int i = 0; i < table->GetNumberOfRows(); i++)
	{
		int ver1 = arrayID1->GetVariantValue(i).ToLong();
		int ver2 = arrayID2->GetVariantValue(i).ToLong();
		std::map< int, int>::iterator iter1 = this->indMapFromVertexToInd.find( ver1);
		std::map< int, int>::iterator iter2 = this->indMapFromVertexToInd.find( ver2);
		if( iter1 != this->indMapFromVertexToInd.end() && iter2 != this->indMapFromVertexToInd.end())
		{
			int index1 = iter1->second;
			int index2 = iter2->second;
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

	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
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

void GraphWindowForNewSelection::SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<int>* clusterSize,
							   std::vector<double> *percentVec, std::vector<double> *disVec, QString filename)
{
	this->fileName = filename;
	vtkAbstractArray *arrayID1 = table->GetColumnByName( ID1.c_str());
	vtkAbstractArray *arrayID2 = table->GetColumnByName( ID2.c_str());

	size.clear();
	if( clusterSize)
	{
		size = *clusterSize;
	}
	else
	{
		for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
		{
			size.push_back(1);
		}
	}

	if( percentVec)
	{
		percentVector.set_size(percentVec->size());
		for( int i = 0; i < percentVec->size(); i++)
		{
			double color = (*percentVec)[i];
			percentVector[i] = color;
		}
	}
	else
	{
		percentVector.set_size(dataTable->GetNumberOfRows());
		for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
		{
			percentVector[i] = 1;
		}
	}
	featureColorVector = percentVector;


	vnl_vector<double> distanceVec;
	if(disVec) 
	{
		distanceVec.set_size( disVec->size());
		for( int i = 0; i < disVec->size(); i++)
		{
			distanceVec[i] = (*disVec)[i];
		}
	}
	else
	{
		distanceVec.set_size(dataTable->GetNumberOfRows());
		for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
		{
			distanceVec[i] = -1;
		}
	}

	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkIntArray> vertexIDarrays = vtkSmartPointer<vtkIntArray>::New();
	vertexIDarrays->SetNumberOfComponents(1);
	vertexIDarrays->SetName("vertexIDarrays");

	vtkSmartPointer<vtkStringArray> vertexLabel = vtkSmartPointer<vtkStringArray>::New();
	vertexLabel->SetNumberOfComponents(1);
	vertexLabel->SetName("vertexLabel");

	// Create the edge weight array
	vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
	weights->SetNumberOfComponents(1);
	weights->SetName("edgeLabel");

	std::cout<< "construct graph"<<endl;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkMutableUndirectedGraph::New();
	graph->GetVertexData()->SetPedigreeIds( this->clusterIds);   // set pedigree ids
	vnl_matrix<int> adj_matrix( dataTable->GetNumberOfRows(), dataTable->GetNumberOfRows());
	
	vnl_matrix<double> vertextList( table->GetNumberOfRows(), 3);
	vnl_vector<double> edgeWeights( table->GetNumberOfRows());
	std::vector< int> size;

	if( clusterSize == NULL)
	{
		for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
		{
			size.push_back(1);
		}
	}
	else
	{
		for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
		{
			size.push_back((*clusterSize)[i]);
		}
	}

	for( int i = 0; i < dataTable->GetNumberOfRows(); i++)
	{
		int vertexID = graph->AddVertex(this->clusterIds->GetValue(i));

		vertexIDarrays->InsertNextValue( size[i]);   // which should be the cluster index

		int per = percentVector[i] * 100 + 0.5;
		int disper = distanceVec[i] * 100 + 0.5;
		if( per > 100)
		{
			per = 100;
		}
		if( disper > 100)
		{
			disper = 100;
		}
		QString strPercent = QString::number(per);

		if( disper > 0 && disper <= 100)
		{
			QString disPercent = QString::number(disper);
			strPercent = "(" + strPercent + "%," + disPercent + "%)";
		}
		else
		{
			strPercent = "(" + strPercent + "%)";
		}
		
		QString str = QString::number(size[i]) + strPercent;
		vertexLabel->InsertNextValue(str.toUtf8().constData());
	}

	edgeMapping.clear();

	for( int i = 0; i < table->GetNumberOfRows(); i++) 
	{
		int ver1 = arrayID1->GetVariantValue(i).ToInt();
		int ver2 = arrayID2->GetVariantValue(i).ToInt();

		graph->AddEdge( ver1, ver2);

		std::pair< int, int> pair1 = std::pair< int, int>( ver1, ver2);
		std::pair< int, int> pair2 = std::pair< int, int>( ver2, ver1);
		edgeMapping.insert( std::pair< std::pair< int, int>, int>(pair1, i));
		edgeMapping.insert( std::pair< std::pair< int, int>, int>(pair2, i));

		adj_matrix( ver1, ver2) = 1;
		adj_matrix( ver2, ver1) = 1;
		vertextList( i, 0) = ver1;
		vertextList( i, 1) = ver2;
		edgeWeights[i] = table->GetValueByName(i, edgeLabel.c_str()).ToDouble();
		weights->InsertNextValue(edgeWeights[i]);
	}

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
	this->lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	this->lookupTable->SetNumberOfTableValues( dataTable->GetNumberOfRows());
	for( int i = 0; i < dataTable->GetNumberOfRows(); i++)               
	{             
		vertexColors->InsertNextValue( i);
		int k = int( featureColorVector[i] * COLOR_MAP2_SIZE + 0.5);
		if( k >= COLOR_MAP2_SIZE)
		{
			k = COLOR_MAP2_SIZE - 1;
		}
		this->lookupTable->SetTableValue(i, COLORMAP2[k].r, COLORMAP2[k].g, COLORMAP2[k].b); 
	}
	lookupTable->Build();

	vtkSmartPointer<vtkIntArray> edgeColors = vtkSmartPointer<vtkIntArray>::New();
	edgeColors->SetNumberOfComponents(1);
	edgeColors->SetName("EdgeColor");
	this->edgeLookupTable = vtkSmartPointer<vtkLookupTable>::New();
	this->edgeLookupTable->SetNumberOfTableValues( table->GetNumberOfRows());
	for( int i = 0; i < table->GetNumberOfRows(); i++)               
	{ 
		edgeColors->InsertNextValue(i); // color the edges by default color
		this->edgeLookupTable->SetTableValue(i, edgeDefaultColor[0], edgeDefaultColor[1], edgeDefaultColor[2]);
	}

	graph->GetVertexData()->AddArray(vertexColors);
	graph->GetVertexData()->AddArray(vertexIDarrays);
	graph->GetVertexData()->AddArray(vertexLabel);
	graph->GetEdgeData()->AddArray(weights);
	graph->GetEdgeData()->AddArray(edgeColors);
	
	theme->SetLineWidth(3);
	theme->SetCellOpacity(0.9);
	theme->SetCellAlphaRange(0.8,0.8);
	theme->SetPointSize(5);
	theme->SetSelectedCellColor(selectColor);
	theme->SetSelectedPointColor(selectColor); 
	theme->SetVertexLabelColor(0.3,0.3,0.3);
	theme->SetBackgroundColor(1,1,1); 
	theme->SetBackgroundColor2(1,1,1);

	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetTableRange (0, 1);
	scalarbarLut->SetNumberOfTableValues(COLOR_MAP2_SIZE);
	for(int index = 0; index < COLOR_MAP2_SIZE; index++)
	{
		rgb rgbscalar = COLORMAP2[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(11);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize(10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->SetMaximumHeightInPixels(1000);
	scalarBar->SetMaximumWidthInPixels(100);
	this->view->GetRenderer()->AddActor2D(scalarBar);

	this->view->RemoveAllRepresentations();
	this->view->SetRepresentationFromInput( graph);
	this->view->GetRepresentation()->SetAnnotationLink( this->ClusterSelections->ClusterAnnotationLink);
	this->view->GetRepresentation()->SetSelectionType( vtkSelectionNode::PEDIGREEIDS);
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetColorVertices(true); 
	this->view->SetColorEdges(true);
	this->view->SetVertexLabelVisibility(true);

	this->view->SetVertexColorArrayName("Color");
	this->view->SetEdgeColorArrayName("EdgeColor");
	this->view->SetEdgeLabelArrayName("edgeLabel");
	this->view->SetVertexLabelArrayName("vertexLabel");
	this->view->SetVertexLabelFontSize(5);

    theme->SetPointLookupTable(lookupTable);
	theme->SetCellLookupTable(edgeLookupTable);
	this->view->ApplyViewTheme(theme);

	this->view->SetLayoutStrategyToPassThrough();
	this->view->SetVertexLabelFontSize(15);
}

void GraphWindowForNewSelection::UpdateCoordinatesByEdgeWeights(std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList)
{
	int oldFirstIndex = backbones[0];
	int oldSecondIndex = 0;
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
void GraphWindowForNewSelection::UpdateChainPointList(int attachnode, std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList)
{
	for( int i = 0; i < chainList.size(); i++)
	{
		std::pair<int, std::vector<int> > pair = chainList[i];
		if( pair.first == attachnode)
		{
 			Point oldFirst = oldPointList[ attachnode];
			Point newFirst = newPointList[ attachnode];
			std::vector<int> vec = pair.second;	
			Point oldSecond(0,0);
			int oldFirstIndex = attachnode; 
			int oldSecondIndex = 0;
			for( int j = 0; j < vec.size(); j++)
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

double GraphWindowForNewSelection::GetEdgeWeight( vnl_matrix<double>& vertexList, long firstIndex, long secondIndex)
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
Point GraphWindowForNewSelection::GetNewPointFromOldPoint( Point &oldPointFirst, Point &oldPointSecond, Point &newPointFirst, double weight)
{
	Point newpoint(0,0);
	double xd = (oldPointFirst.x - oldPointSecond.x) * (oldPointFirst.x - oldPointSecond.x);
	double yd = (oldPointFirst.y - oldPointSecond.y) * (oldPointFirst.y - oldPointSecond.y);
	double distance = sqrt( xd + yd);

	newpoint.x = newPointFirst.x + ( oldPointSecond.x - oldPointFirst.x) * weight / distance;
	newpoint.y = newPointFirst.y + ( oldPointSecond.y - oldPointFirst.y) * weight /  distance;
	return newpoint;
}

void GraphWindowForNewSelection::ShowGraphWindow()
{
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();
	std::cout<< "view->ResetCamera"<<endl;
	view->ResetCamera();

	if( cornAnnotation == NULL)
	{
		cornAnnotation = vtkSmartPointer<vtkCornerAnnotation>::New();
		cornAnnotation->SetText(ANNOTATION_CORNER, "Colored by Device Sample Percentage");
		cornAnnotation->GetTextProperty()->SetFontSize(10);
		cornAnnotation->GetTextProperty()->SetColor(0,0,0);
		this->view->GetRenderer()->AddActor2D(cornAnnotation);
	}
	else
	{
		cornAnnotation->SetText(ANNOTATION_CORNER, "Colored by Device Sample Percentage");
		cornAnnotation->GetTextProperty()->SetFontSize(10);
		cornAnnotation->GetTextProperty()->SetColor(0,0,0);
	}

	std::cout<< "view->Render"<<endl;
	view->Render();
	std::cout<< "Successfully rendered"<<endl;

    this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New(); 
    this->selectionCallback->SetClientData(this);
    this->selectionCallback->SetCallback ( SelectionCallbackFunction);
	vtkAnnotationLink *link = view->GetRepresentation()->GetAnnotationLink();
	link->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);

	this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
	this->keyPress->SetCallback(HandleKeyPress);
	this->keyPress->SetClientData(this);
	this->view->GetInteractor()->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);

	view->GetInteractor()->Start();
}

void GraphWindowForNewSelection::SelectionCallbackFunction(vtkObject* caller, long unsigned eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	GraphWindowForNewSelection* graphWin = (GraphWindowForNewSelection*)clientData;

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
		if( vertexList != NULL && vertexList->GetNumberOfTuples() > 0)
		{
			std::vector<int> ClusterIDs;
			int key = graphWin->view->GetInteractor()->GetControlKey();
			if( 0 == key)
			{
				for( int i = 0; i < vertexList->GetNumberOfTuples(); i++)
				{
					int value = vertexList->GetValue(i);
					ClusterIDs.push_back(value);
				}
				graphWin->SetProgressionStartTag(true);
			}
			else
			{
				int value = vertexList->GetValue(0);
				ClusterIDs.push_back(value);
				graphWin->SetUserDefineProgression(value);
			}
			graphWin->UpdataLookupTable(ClusterIDs);
		}
		else
		{
			graphWin->SetProgressionStartTag(true);
		}
	}
	graphWin->mainQTRenderWidget.GetRenderWindow()->Render();
	graphWin->view->GetRenderer()->Render();
}

void GraphWindowForNewSelection::HandleKeyPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	GraphWindowForNewSelection* graphWin = (GraphWindowForNewSelection*)clientData;
	char key = graphWin->view->GetInteractor()->GetKeyCode();
	switch (key)
	{
	case 'r':
		graphWin->RestoreLookupTable();
		break;
	default:
		break;
	}
}

void GraphWindowForNewSelection::RestoreLookupTable()
{
	for( int i = 0; i < percentVector.size(); i++)               
	{
		featureColorVector[i] = percentVector[i];
		int k = int( featureColorVector[i] * COLOR_MAP2_SIZE + 0.5);
		if( k >= COLOR_MAP2_SIZE)
		{
			k = COLOR_MAP2_SIZE - 1;
		}
		this->lookupTable->SetTableValue(i, COLORMAP2[k].r, COLORMAP2[k].g, COLORMAP2[k].b); 
	}
	lookupTable->Build();
	//this->view->GetRenderer()->RemoveActor2D(cornAnnotation);
	cornAnnotation->SetText(ANNOTATION_CORNER, "Colored by Device Sample Percentage");
	cornAnnotation->GetTextProperty()->SetFontSize(10);
	cornAnnotation->GetTextProperty()->SetColor(0,0,0);
	//this->view->GetRenderer()->AddActor2D(cornAnnotation);
	this->view->GetRenderer()->Render();
}

void GraphWindowForNewSelection::SetProgressionStartTag(bool bstart)
{
	bProgressionStart = bstart;
}

void GraphWindowForNewSelection::SetUserDefineProgression(int nodeID)
{
	if( bProgressionStart)
	{
		progressionStartID = nodeID;
		bProgressionStart = false;
		std::cout<< "progression start: "<<progressionStartID<<endl;
	}
	else
	{
		progressionEndID = nodeID;
		bProgressionStart = true;
		std::cout<< "progression end: "<<progressionEndID<<endl;
		UpdateProgressionPath();
	}
}

void GraphWindowForNewSelection::ResetLookupTable(vtkSmartPointer<vtkLookupTable> lookuptable, double* color)
{
	for( int i = 0; i < lookuptable->GetNumberOfTableValues(); i++)
	{
		lookuptable->SetTableValue( i, color[0], color[1], color[2]);
	}
}

void GraphWindowForNewSelection::UpdataLookupTable( std::vector<int>& IDs)
{
	std::set<int> selectedIDs; 

	for( int i = 0; i < IDs.size(); i++)
	{
		std::map<int, int>::iterator iter2 = this->indMapFromVertexToInd.find( IDs[i]);
		if( iter2 != this->indMapFromVertexToInd.end())
		{
			selectedIDs.insert( iter2->second);
		}
	}

	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		if (selectedIDs.find(i) != selectedIDs.end())
		{
			this->lookupTable->SetTableValue(i, selectColor[0], selectColor[1], selectColor[2]); // color the vertices
		}
		else
		{
			int k = int( featureColorVector[i] * COLOR_MAP2_SIZE + 0.5);
			if( k >= COLOR_MAP2_SIZE)
			{
				k = COLOR_MAP2_SIZE - 1;
			}
			this->lookupTable->SetTableValue(i, COLORMAP2[k].r, COLORMAP2[k].g, COLORMAP2[k].b); // color the vertices- blue
		}
	}
	this->lookupTable->Build();
}

void GraphWindowForNewSelection::UpdateProgressionPath()
{
	progressionPath.clear();
	ResetLookupTable( edgeLookupTable, edgeDefaultColor);
	GetProgressionPath(shortest_hop, progressionStartID, progressionEndID, progressionPath);

	for( int i = 0; i < progressionPath.size() - 1; i++)
	{
		std::pair<int, int> pair = std::pair<int, int>( progressionPath[i], progressionPath[i + 1]);
		std::map< std::pair< int, int>, int>::iterator iter = edgeMapping.find( pair);
		if( iter != edgeMapping.end())
		{
			int index = (int)iter->second;
			edgeLookupTable->SetTableValue( index, 0, 0 ,0);
		}
	}
	//this->mainQTRenderWidget.GetRenderWindow()->Render();
	this->view->GetRenderer()->Render();
}

void GraphWindowForNewSelection::GetProgressionPath(vnl_matrix<int> &hopMat, int startNode, int endNode, std::vector< int> &path)
{
	std::map< int, int> tmpChain;  // for ordering 
	int maxhop = hopMat( startNode, endNode);

	for( int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( startNode, i) + shortest_hop( endNode, i) == maxhop)
		{
			tmpChain.insert( std::pair< int, int>( shortest_hop( startNode, i), i));
		}
	}

	path.clear();
	std::map< int, int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		path.push_back((*iter).second);
	}
}

void GraphWindowForNewSelection::CalculateCoordinates(vnl_matrix<int>& adj_matrix, std::vector<Point>& pointList)
{
	shortest_hop.set_size(adj_matrix.rows(), adj_matrix.cols());
	vnl_vector<int> mark(adj_matrix.rows());
	shortest_hop.fill(0);
	mark.fill(0);
	mark[0] = 1;
	std::vector< int> checkNode;
	std::vector<int> noncheckNode;
	find( mark, 0, noncheckNode, checkNode);

	/// calculate the shortest hop matrix
	while( noncheckNode.size() > 0)
	{
		int checkNodeInd = -1;
		int noncheckNodeInd = -1;
		for( int i = 0; i < checkNode.size(); i++)
		{
			for( int j = 0; j < noncheckNode.size(); j++)
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

		for( int i = 0; i < checkNode.size(); i++)
		{
			int size = (this->size[noncheckNodeInd] + this->size[checkNodeInd]) / 2;
			shortest_hop( checkNode[i], noncheckNodeInd) =  shortest_hop( checkNode[i], checkNodeInd) + size;
			shortest_hop( noncheckNodeInd, checkNode[i]) =  shortest_hop( checkNode[i], checkNodeInd) + size;
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
	int maxhop = shortest_hop.max_value();
	unsigned int maxId = shortest_hop.arg_max();
	unsigned int coln = maxId / shortest_hop.rows();
	unsigned int rown = maxId  - coln * shortest_hop.rows();

	backbones.clear();
	chainList.clear();  // store the chains, no repeat

	std::queue<int> tmpbackbones;   // for calculating all the sidechains
	std::vector<int> debugbackbones;  // for debugging
	std::map< int, int> tmpChain;  // for ordering 
	vnl_vector< int> tag( shortest_hop.rows());     // whether the node has been included in the chain
	tag.fill( 0);

	std::cout<<"build chains"<<endl;
	QString chainStr = this->fileName + "chains.txt";
	ofstream ofChains( chainStr.toStdString().c_str());
	/// find the backbone
	for( int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( coln, i) + shortest_hop( rown, i) == maxhop)
		{
			tmpChain.insert( std::pair< int, int>( shortest_hop( coln, i), i));
		}
	}

	std::map< int, int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		backbones.push_back((*iter).second);
		tmpbackbones.push((*iter).second);
		debugbackbones.push_back( (*iter).second);
		tag[ (*iter).second] = 1;
	}

	/// update the progression path
	progressionPath = backbones;   
	progressionStartID = backbones[0];
	progressionEndID = backbones[ backbones.size() - 1];

	/// find the branches' backbones
	std::vector< int> branchnodes;
	std::vector< int> chains;
	while( !tmpbackbones.empty())
	{
		int ind = tmpbackbones.front(); 
		for( int i = 0; i < adj_matrix.cols(); i++)
		{
			if( adj_matrix( ind, i) != 0 && tag[i] == 0)    // ind's  neighbours that haven't been checked
			{
				branchnodes.clear();
				branchnodes.push_back( ind);
				for( int j = 0; j < shortest_hop.cols(); j++)
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
							for( int st = 0; st < debugbackbones.size(); st++)
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
					chainList.push_back( std::pair< int, std::vector<int> >(ind, chains));
					for( int k = 0; k < chains.size(); k++)
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
	for( int i = 0; i < backbones.size(); i++)
	{
		ofChains << backbones[i]<<"\t";
	}
	ofChains <<endl;
	ofChains << "branch chains:"<<endl;
	std::vector< std::pair< int, std::vector<int> > >::iterator chainIter;
	for( chainIter = chainList.begin(); chainIter != chainList.end(); chainIter++)
	{
		std::pair< int, std::vector< int> >tmp = *chainIter;
		ofChains << tmp.first;
		std::vector< int> branchlist = tmp.second;
		for( int i = 0; i < branchlist.size(); i++)
		{
			ofChains << "\t"<< branchlist[i];
		}
		ofChains <<endl;
	}

	SortChainList( shortest_hop, backbones, chainList);   // the order to draw the subbones
	//for( chainIter = chainList.begin(); chainIter != chainList.end(); chainIter++)
	//{
	//	std::pair< int, std::vector< int> >tmp = *chainIter;
	//	ofChains << tmp.first;
	//	std::vector< int> branchlist = tmp.second;
	//	for( int i = 0; i < branchlist.size(); i++)
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
	for( int i = 0; i < backboneAngle.size(); i++)
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
		int size = backbones.size();
		for( int i = 0; i < backbones.size(); i++)
		{
			nodePos(0, backbones[i]) = sin( backboneAngle[i]) / norm * 500 / size;
			nodePos(1, backbones[i]) = -cos( backboneAngle[i]) / norm * 500 / size;
			nodePosAssigned[ backbones[i]] = 1;
		}
	}
	else
	{
		for( int i = 0; i < backbones.size(); i++)
		{
			nodePos(0, backbones[i]) = sin( backboneAngle[i]) / norm;
			nodePos(1, backbones[i]) = -cos( backboneAngle[i]) / norm;
			nodePosAssigned[ backbones[i]] = 1;
		}
	}

	/// cacluate the branch backbone nodes position
	for( int i = 0; i < chainList.size(); i++)
	{
		if( i % 100 == 0)
		{
			std::cout<<i<<endl;
		}
		std::pair< int, std::vector<int> > branch = chainList[i];
		int attachNode = branch.first;
		std::vector< int> branchNode = branch.second;
		for( int j = 0; j < branchNode.size(); j++)
		{
			int newNode = branchNode[j];
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
					for( int k = 0; k < repel_mat.cols(); k++)
					{
						repel_mat(0, k) = newNodePoint.x - repel_mat( 0, k);
						repel_mat(1, k) = newNodePoint.y - repel_mat( 1, k);
					}

					vnl_vector< double> repel_tmp_xvector(repel_mat.cols());
					vnl_vector< double> repel_tmp_yvector(repel_mat.cols());
					for( int k = 0; k < repel_mat.cols(); k++)
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
		for( int j = 0; j < tmpNodePos.size(); j++)
		{
			tmp[j] = abs( tmp[j]);
		}
		
		tmpNodePos = tmpNodePos / tmp.max_value() * 50;
		nodePos.set_row(i, tmpNodePos);
	}

	ofCoordinate << setiosflags(ios::fixed)<< nodePos.transpose()<<endl;
	ofCoordinate.close();

	for( int i = 0; i < nodePos.cols(); i++)
	{
		Point pt( nodePos(0, i), nodePos(1, i));
		pointList.push_back( pt);
	}
}

/// branchnodes first element is the chain's attach node
void GraphWindowForNewSelection::getBackBones(vnl_matrix< int>& shortest_hop, std::vector< int>& branchnodes, std::vector< int>& chains)
{
	chains.clear();
	vnl_vector< int> branchShortestHop( branchnodes.size());

	int attachNode = branchnodes.front();
	std::map< int, int> tmpChain;  // for ordering 

	for( int i = 0; i < branchShortestHop.size(); i++)
	{
		branchShortestHop[i] = shortest_hop( attachNode, branchnodes[i]);
	}

	int maxHop = branchShortestHop.max_value();
	int endNodeIndex = branchShortestHop.arg_max();
	int endNode = branchnodes[ endNodeIndex];
	for( int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( attachNode, i) + shortest_hop( endNode, i) == maxHop)
		{
			tmpChain.insert( std::pair< int, int>( shortest_hop( attachNode, i), i));
		}
	}

	std::map< int, int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		if( (*iter).second != attachNode)
		{
			chains.push_back((*iter).second);
		}
	}
}

void GraphWindowForNewSelection::find(vnl_vector<int>& vec, int val, std::vector<int>& equal, std::vector<int>& nonequal)
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

void GraphWindowForNewSelection::GetElementsIndexInMatrix(vnl_matrix<int>& mat, int rownum, int max, vnl_matrix<double>& oldmat, vnl_matrix<double>& newmat, vnl_vector< int>& tag)
{
	std::vector< int> index;
	for( int i = 0; i < mat.cols(); i++)
	{
		if( mat( rownum, i) < max && tag[i] == 1)
		{
			index.push_back( i);
		}
	}

	newmat.set_size( oldmat.rows(), index.size());
	for( int i = 0; i < index.size(); i++)
	{
		vnl_vector< double> tmpcol = oldmat.get_column( index[i]);
		newmat.set_column( i, tmpcol);
	}
}

double GraphWindowForNewSelection::Median( vnl_vector<double> vec)
{
	vnl_vector<double> vect = vec;
	vnl_vector<double> tmp( vec.size());
	double max = vect.max_value() + 1;
	double med;
	for( int i = 0; i < vect.size(); i++)
	{
		tmp[i] = vect.min_value();
		int ind = vect.arg_min();
		vect[ind] = max;
	}

	int size = tmp.size();
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

void GraphWindowForNewSelection::SortChainList( vnl_matrix<int>& shortest_hop, std::vector<int>& backbones, 
				   std::vector< std::pair<int, std::vector<int> > >& chainList)
{
	std::vector< std::pair<int, std::vector<int> > > tmpchainList;
	std::multimap< int, int> sortMap;
	int startNode = backbones[0];
	int endNode = backbones[ backbones.size() - 1];
	for( int i = 0; i < chainList.size(); i++)
	{
		std::pair<int, std::vector<int> > chainPair = chainList[i];
		if( IsExist( backbones, chainPair.first))
		{
			int abHop = abs( (int)shortest_hop(chainPair.first, startNode) - (int)shortest_hop( chainPair.first, endNode));
			sortMap.insert( std::pair< int, int>(abHop, i));
		}
		else
		{
			std::multimap< int, int>::iterator iter;
			for( iter = sortMap.begin(); iter != sortMap.end(); iter++)
			{
				std::pair<int, std::vector<int> > pair = chainList[ (*iter).second];
				tmpchainList.push_back( pair);
			}
			for( int k = i; k < chainList.size(); k++)
			{
				tmpchainList.push_back( chainList[k]);
			}
			chainList = tmpchainList;
			break;
		}
	}
}

bool GraphWindowForNewSelection::IsExist(std::vector<int>& vec, int value)
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

void GraphWindowForNewSelection::GetProgressionTreeOrder(std::vector<long int> &order)
{
	order.clear();
	for( int i = 0; i < backbones.size(); i++)
	{
		order.push_back(backbones[i]);
		GetOrder(backbones[i], order);
	}
}

void GraphWindowForNewSelection::GetOrder(int node, std::vector<long int> &order)
{
	std::vector<int> vec;
	for( int i = 0; i < chainList.size(); i++)
	{
		std::pair<int, std::vector<int> > pair = chainList[i];
		if( pair.first == node)
		{
			vec = pair.second;
			for( int j = 0; j < vec.size(); j++)
			{
				order.push_back( vec[j]);
				GetOrder( vec[j], order);
			}
		}
	}
	vec.clear();
}

void GraphWindowForNewSelection::ColorTreeAccordingToFeatures(vnl_vector<double> &feature, const char *featureName)
{
	featureColorVector = feature;
	featureColorVector -= feature.mean();
	double var = featureColorVector.two_norm();
	featureColorVector /= var;

	double max = featureColorVector.max_value();
	double min = featureColorVector.min_value();
	featureColorVector = (featureColorVector - min) / ( max - min);
	
	for( int i = 0; i < featureColorVector.size(); i++)               
	{
		int k = int( featureColorVector[i] * COLOR_MAP2_SIZE + 0.5);
		if( k >= COLOR_MAP2_SIZE)
		{
			k = COLOR_MAP2_SIZE - 1;
		}
		this->lookupTable->SetTableValue(i, COLORMAP2[k].r, COLORMAP2[k].g, COLORMAP2[k].b); 
	}
	lookupTable->Build();

	//this->view->GetRenderer()->RemoveActor2D(cornAnnotation);
	cornAnnotation->SetText(ANNOTATION_CORNER, featureName);
	//this->view->GetRenderer()->AddActor2D(cornAnnotation);
	cornAnnotation->GetTextProperty()->SetFontSize(10);
	cornAnnotation->GetTextProperty()->SetColor(0,0,0);
	this->view->GetRenderer()->Render();
	this->mainQTRenderWidget.GetRenderWindow()->Render();
}

void GraphWindowForNewSelection::closeEvent(QCloseEvent *event)
{
	mainQTRenderWidget.close();
}