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

GraphWindow::GraphWindow(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;// = new QVTKWidget;
	this->TTG = vtkSmartPointer<vtkTableToGraph>::New();
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->observerTag = 0;
	this->lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	this->points =  vtkSmartPointer<vtkPoints>::New();
}

GraphWindow::~GraphWindow()
{
}

void GraphWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels)
{
	this->dataTable = table;
	for( long int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		long int var = this->dataTable->GetValue( i, 0).ToLong();
		this->indMapFromVertexToInd.insert( std::pair< long int, long int>(var, i));
		this->indMapFromIndToVertex.push_back( var);
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

	//this->TTG->ClearLinkVertices();
	//this->TTG->SetInput(0, table);
	//this->TTG->AddLinkEdge(ID1.c_str(), ID2.c_str()); 

	vtkAbstractArray *arrayID1 = table->GetColumnByName( ID1.c_str());
	vtkAbstractArray *arrayID2 = table->GetColumnByName( ID2.c_str());

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
		vertexIDarrays->InsertNextValue( this->indMapFromIndToVertex[i]);
	}

	for( int i = 0; i < table->GetNumberOfRows(); i++)
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
	graph->GetVertexData()->AddArray(vertexIDarrays);
	graph->GetEdgeData()->AddArray(weights);

	this->view->AddRepresentationFromInput( graph);
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetColorVertices(true); 
	this->view->SetVertexLabelVisibility(true);

	this->view->SetVertexColorArrayName("Color");

    theme->SetPointLookupTable(lookupTable);
    theme->SetBackgroundColor(0,0,0); 
	this->view->ApplyViewTheme(theme);

	this->view->SetEdgeLabelArrayName("edgeLabel");

	this->view->SetLayoutStrategyToForceDirected();

	this->view->SetVertexLabelArrayName("vertexIDarrays");
	this->view->SetVertexLabelFontSize(20);
}

void GraphWindow::SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel)
{
	vtkAbstractArray *arrayID1 = table->GetColumnByName( ID1.c_str());
	vtkAbstractArray *arrayID2 = table->GetColumnByName( ID2.c_str());

	vtkSmartPointer<vtkIntArray> vertexIDarrays = vtkSmartPointer<vtkIntArray>::New();
	vertexIDarrays->SetNumberOfComponents(1);
	vertexIDarrays->SetName("vertexIDarrays");

  // Create the edge weight array
	vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
	weights->SetNumberOfComponents(1);
	weights->SetName("edgeLabel");

	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkMutableUndirectedGraph::New();
	vnl_matrix<long int> adj_matrix( this->dataTable->GetNumberOfRows(), this->dataTable->GetNumberOfRows());

	for( int i = 0; i <  this->dataTable->GetNumberOfRows(); i++)
	{
		int vertexID = graph->AddVertex();
		vertexIDarrays->InsertNextValue( this->indMapFromIndToVertex[i]);
	}

	for( int i = 0; i < table->GetNumberOfRows(); i++)
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
			adj_matrix( index1, index2) = 1;
			adj_matrix( index2, index1) = 1;
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

	std::vector<Point> pointList;
	CalculateCoordinates(adj_matrix, pointList);
	if( pointList.size() > 0)
	{
		for( int i = 0; i <  pointList.size(); i++)
		{
			this->points->InsertNextPoint(pointList[i].x, pointList[i].y, 0);
		}
		graph->SetPoints( this->points);
	}

	vtkSmartPointer<vtkIntArray> vertexColors = vtkSmartPointer<vtkIntArray>::New();
	vertexColors->SetNumberOfComponents( this->dataTable->GetNumberOfRows());
	vertexColors->SetName("Color");
	
	this->lookupTable->SetNumberOfTableValues( this->dataTable->GetNumberOfRows());
	for( int i = 0; i < this->dataTable->GetNumberOfRows(); i++)
	{
		vertexColors->InsertNextValue( i);
		this->lookupTable->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
    }
	lookupTable->Build();
	
	this->theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	this->theme->SetLineWidth(5);
	this->theme->SetCellOpacity(0.9);
	this->theme->SetCellAlphaRange(0.8,0.8);
	this->theme->SetPointSize(8);
	this->theme->SetSelectedCellColor(1,0,0);
	this->theme->SetSelectedPointColor(1,0,0); 
    this->theme->SetPointLookupTable(lookupTable);
    this->theme->SetBackgroundColor(0,0,0); 
	this->view->ApplyViewTheme(theme);
	
	graph->GetVertexData()->AddArray(vertexColors);
	graph->GetVertexData()->AddArray(vertexIDarrays);
	graph->GetEdgeData()->AddArray(weights);

	this->view->AddRepresentationFromInput( graph);
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName("edgeLabel");

	this->view->SetVertexLabelVisibility(true);
	this->view->SetColorVertices(true); 
	this->view->SetVertexLabelArrayName("vertexIDarrays");
	this->view->SetVertexColorArrayName("Color");
	this->view->SetVertexLabelFontSize(20);
	this->view->SetLayoutStrategyToPassThrough();
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

	graphWin->mainQTRenderWidget.GetRenderWindow()->Render();
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
			long int var = indMapFromIndToVertex[*iter];
			selectedIDs.insert( var);
			iter++;
		}
		this->selection->select( selectedIDs);
		UpdataLookupTable( selectedIDs);
	}
	this->view->GetRenderer()->Render();
}

void GraphWindow::UpdataLookupTable( std::set<long int>& IDs)
{
	std::set<long int> selectedIDs; 
	std::set<long int>::iterator iter = IDs.begin();

	while( iter != IDs.end())
	{
		long int var = this->indMapFromVertexToInd.find( *iter)->second;
		selectedIDs.insert( var);
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
		std::set<long int> updataIDes =  this->selection->getSelections();
		UpdataLookupTable( updataIDes);

		this->mainQTRenderWidget.GetRenderWindow()->Render();
		//this->view->GetRenderer()->Render();
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
			shortest_hop( checkNode[i], noncheckNodeInd) =  shortest_hop( checkNode[i], checkNodeInd) + 1;
			shortest_hop( noncheckNodeInd, checkNode[i]) =  shortest_hop( checkNode[i], checkNodeInd) + 1;
		}

		mark[ noncheckNodeInd] = 1;

		checkNode.clear();
		noncheckNode.clear();
		find( mark, 0, noncheckNode, checkNode);
	}

	ofstream ofs("shortest_hop.txt");
	ofs<< shortest_hop<<endl;
	ofs.close();

	/// find the root and chains of the tree
	long int maxhop = shortest_hop.max_value();
	unsigned int maxId = shortest_hop.arg_max();
	unsigned int coln = maxId / shortest_hop.rows();
	unsigned int rown = maxId  - coln * shortest_hop.rows();

	std::vector<long int> backbones;
	std::queue<long int> tmpbackbones;   // for calculating all the sidechains
	std::map< long int, long int> tmpChain;  // for ordering 
	std::multimap< long int, std::vector<long int> > chainList;  // store the chains, no repeat
	vnl_vector< int> tag( shortest_hop.rows());     // whether the node has been included in the chain
	tag.fill( 0);

	/// find the backbone
	for( long int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( coln, i) + shortest_hop( rown, i) == maxhop)
		{
			tmpChain.insert( std::pair< long int, long int>( shortest_hop( rown, i), i));
		}
	}

	std::map< long int, long int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		backbones.push_back((*iter).second);
		tmpbackbones.push((*iter).second);
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
				for( long int j = 0; j < shortest_hop.cols(); j++)
				{
					if( shortest_hop( ind, j) > shortest_hop( i, j))  // find neighbour i's branch nodes including i
					{
						if( tag[j] == 0)
						{
							branchnodes.push_back( j);
						}
					}
				}

				if( branchnodes.size() > 0)
				{
					getBackBones( shortest_hop, tag, branchnodes, chains);
					chainList.insert( std::pair< long int, std::vector<long int> >(ind, chains));
					for( long int k = 0; k < chains.size(); k++)
					{
						tmpbackbones.push( chains[k]);
					}
					branchnodes.clear();
				}
			}
		}
		tmpbackbones.pop();
	}

	ofstream ofChains("chains.txt");
	ofChains << "backbones:" <<endl;
	for( long int i = 0; i < backbones.size(); i++)
	{
		ofChains << backbones[i]<<"\t";
	}
	ofChains <<endl;
	ofChains << "branch chains:"<<endl;
	std::multimap< long int, std::vector<long int> >::iterator chainIter;
	for( chainIter = chainList.begin(); chainIter != chainList.end(); chainIter++)
	{
		std::vector< long int> tmp = (*chainIter).second;
		ofChains << (*chainIter).first;
		for( long int i = 0; i < tmp.size(); i++)
		{
			ofChains << "\t"<< tmp[i];
		}
		ofChains <<endl;
	}
	ofChains.close();
}

/// branchnodes first element is the chain's first node
void GraphWindow::getBackBones(vnl_matrix< long int>& shortest_hop, vnl_vector< int>& tag, 
							   std::vector< long int>& branchnodes, std::vector< long int>& chains)
{
	chains.clear();
	vnl_vector< long int> branchShortestHop( branchnodes.size());

	long int startNode = branchnodes.front();
	std::map< long int, long int> tmpChain;  // for ordering 

	for( long int i = 0; i < branchShortestHop.size(); i++)
	{
		branchShortestHop[i] = shortest_hop( startNode, branchnodes[i]);
	}

	long int maxHop = branchShortestHop.max_value();
	long int endNodeIndex = branchShortestHop.arg_max();
	long int endNode = branchnodes[ endNodeIndex];
	for( long int i = 0; i < shortest_hop.cols(); i++)
	{
		if( shortest_hop( startNode, i) + shortest_hop( endNode, i) == maxHop)
		{
			tmpChain.insert( std::pair< long int, long int>( shortest_hop( startNode, i), i));
		}
	}

	std::map< long int, long int>::iterator iter;
	for( iter = tmpChain.begin(); iter != tmpChain.end(); iter++)
	{
		chains.push_back((*iter).second);
		tag[ (*iter).second] = 1;
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
