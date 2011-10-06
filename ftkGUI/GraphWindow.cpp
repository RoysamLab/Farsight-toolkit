#include "GraphWindow.h"
GraphWindow::GraphWindow(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;// = new QVTKWidget;
	this->TTG = vtkSmartPointer<vtkTableToGraph>::New();
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
}
GraphWindow::~GraphWindow()
{
}
void GraphWindow::setQtModels(QItemSelectionModel *mod)
{
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

	
	this->theme.TakeReference(vtkViewTheme::CreateMellowTheme());
	this->theme->SetLineWidth(5);
	this->theme->SetCellOpacity(0.9);
	this->theme->SetCellAlphaRange(0.5,0.5);
	this->theme->SetPointSize(10);
	this->theme->SetSelectedCellColor(1,0,1);
	this->theme->SetSelectedPointColor(1,0,1); 

	
	this->view->AddRepresentationFromInputConnection(TTG->GetOutputPort());
	this->view->SetEdgeLabelVisibility(true);
	this->view->SetEdgeLabelArrayName(edgeLabel.c_str());
	this->view->SetLayoutStrategyToForceDirected();
	this->view->SetVertexLabelArrayName("label");
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
	view->GetInteractor()->Start();
}