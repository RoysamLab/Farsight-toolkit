#include "Dendrogram_Selection.h"
#include "vtkSmartPointer.h"
#include "vtkGraphLayoutView.h"



Dendrogram_Selection::Dendrogram_Selection(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{

	dendro= new Dendrogram();
	QVTK = new QVTKWidget(this);
	setCentralWidget(QVTK);
	Rubber_Band = vtkSmartPointer<vtkInteractorStyleRubberBand2D >::New();
	
	//vtkSmartPointer<vtkGraphLayoutView> graphLayoutView= dendro->GetGraphLayoutView();	
	QVTK->SetRenderWindow(dendro->GetGraphLayoutView()->GetRenderWindow());
	//QVTK->GetRenderWindow()->AddRenderer(dendro->GetGraphLayoutView()->GetRenderer());
	//QVTK->GetInteractor()->SetRenderWindow(QVTK->GetRenderWindow());
	//QVTK->GetRenderWindow()->SetInteractor(dendro->GetGraphLayoutView()->GetInteractor());
	QVTK->GetInteractor()->SetInteractorStyle(Rubber_Band);
	//QVTK->GetInteractor()->S
	
	//dendro->GetGraphLayoutView()->GetInteractor()->Start();
	/*selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	this->selectionCallback->SetClientData(dendro);
	this->selectionCallback->SetCallback (this->SelectionCallbackFunction);*/



	//QVTK->G
	//QVTK->GetRenderWindow()->Render();
	//dendro->GetGraphLayoutView()->GetRenderWindow()->Render();
	//QVTK->GetRenderWindow()->GetInteractor()->AddObserver("AnotationChangedEvent", this->selectionCallback);
	//this->dendro->GetGraphLayoutView()->GetRepresentation()->GetAnnotationLink()->
	
	/////////////////////////////////////////////////////////
	/*(dendro->GetGraphLayoutView())->GetInteractor()->Start();*/

	//QVTK->GetInteractor()->Start();
}


//void Dendrogram::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
//{
// 
//	Dendrogram* Dendro = (Dendrogram*)clientData;
//	vtkAnnotationLink* annotationLink =
//    static_cast<vtkAnnotationLink*>(caller);
//    vtkSmartPointer<vtkDataObject> view =
//    vtkSmartPointer<vtkDataObject>::New();
//	vtkSelection* selection = annotationLink->GetCurrentSelection();
//  
//    vtkSelectionNode* vertices;
//	vtkSelectionNode* edges;
//	if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
//		{
//		vertices = selection->GetNode(0);
//		}
//	else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
//		{
//		edges = selection->GetNode(0);
//		}
// 
//	if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
//		{
//		vertices = selection->GetNode(1);
//		}
//	else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
//		{
//		edges = selection->GetNode(1);
//		}
// 
//	vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
//	//std::cout << "There are " << vertexList->GetNumberOfTuples() << " vertices selected." << std::endl;
// 
//	if(vertexList->GetNumberOfTuples() > 0)
//		{
//	    //std::cout << "Vertex Ids: ";
//		}
//
//	if(vertexList->GetNumberOfTuples()==0)
//		{
//			for(int n=0;n<Dendro->num_row;n++)
//			{
//			Dendro->lookupTable->SetTableValue(n, 1.0, 1.0, 1.0);
//			}
//		}
//	else
//		{
//			for(vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
//			{
//			//std::cout << vertexList->GetValue(i) << " ";
//			vtkIdType v=vertexList->GetValue(i);
//				if(v>(Dendro->num_row-1))
//					{
//	 
//						int p=v;
//						if(p>(Dendro->num_row-1))
//							{
//							int j=0;
//							for(int i=0;i<(Dendro->num_row);i++)
//								Dendro->colour_child[i]=0;
//
//				
//							while((Dendro->CharLabel[p-(Dendro->num_row-1)][j])!=-1)
//								{
//								int k;
//								///std::cout << "The childs are "<<(Dendro->CharLabel[p-(Dendro->num_row-1)][j])<<"\t";
//					
//								for(k=0;k<Dendro->num_row;k++)
//									{
//									double rep= (Dendro->Tree3D[0][k][0]);
//									if((rep)==(Dendro->CharLabel[p-(Dendro->num_row-1)][j]))
//										{
//										Dendro->colour_child[k]=1;
//										}
//				
//									}
//								for(int n=0;n<Dendro->num_row;n++)
//									{
//										if(Dendro->colour_child[n]==1)
//											Dendro->lookupTable->SetTableValue(n, 1.0, 0.0, 1.0);
//										else
//											Dendro->lookupTable->SetTableValue(n, 1.0, 1.0, 1.0);
//									}
//			
//								j++;
//							}
//
//			}	
//	
//    }
//	else
//		for(int n=0;n<Dendro->num_row;n++)
//		Dendro->lookupTable->SetTableValue(n, 1.0, 1.0, 1.0);
//  }
// }
//  
//  vtkIdTypeArray* edgeList = vtkIdTypeArray::SafeDownCast(edges->GetSelectionList());
// 
//  
// 
//}



Dendrogram_Selection::~Dendrogram_Selection()
{
}
void Dendrogram_Selection::updateSelection()
{
}


