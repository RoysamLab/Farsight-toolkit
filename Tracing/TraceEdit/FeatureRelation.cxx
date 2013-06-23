#include "FeatureRelation.h"

FeatureRelation::FeatureRelation(void)
{
	mainQTRenderWidget = new QVTKWidget;
}

FeatureRelation::~FeatureRelation(void)
{
}

void FeatureRelation::FeatureGraph()
{
	vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
	reader->SetFileName("C:/Audrey's lab junk/Lmeasures/LmeasureSemanticgroups4.txt");
	reader->SetFieldDelimiterCharacters("	");
	reader->SetHaveHeaders(true);
	reader->Update();

	vtkTable* table = reader->GetOutput();

	//std::cout << "Table has " << table->GetNumberOfRows() << " rows." << std::endl;
	//std::cout << "Table has " << table->GetNumberOfColumns() << " columns." << std::endl;

	//for (int i = 0; i < 10; i++)
	//	std::cout << table->GetColumnName(i) << std::endl;

	vtkSmartPointer<vtkTableToTreeFilter> tableToTree = vtkSmartPointer<vtkTableToTreeFilter>::New();
	//vtkTableToTree * tableToTree = vtkTableToTreeFilter::New();
	tableToTree->SetInputData(table);
	tableToTree->Update();

	vtkTree * tree = tableToTree->GetOutput();

	vtkSmartPointer<vtkGroupLeafVertices> group = vtkSmartPointer<vtkGroupLeafVertices>::New();
	group->SetInputConnection(tableToTree->GetOutputPort());
	group->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "Features");
	group->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES , "Semantic");

	tree = group->GetOutput();


	//vtkSmartPointer<vtkTreeLayoutStrategy> strategy = vtkSmartPointer<vtkTreeLayoutStrategy>::New();
	//strategy->SetRadial(true);
	//strategy->SetAngle(360);

	/*vtkSmartPointer<vtkGraphLayoutView> view = vtkSmartPointer<vtkGraphLayoutView>::New();
	view->AddRepresentationFromInput(group->GetOutput());
	view->SetLayoutStrategyToTree();*/


	vtkSmartPointer<vtkTableToGraph> tableToGraph = vtkSmartPointer<vtkTableToGraph>::New();
	tableToGraph->SetInputConnection(reader->GetOutputPort());

	//tableToGraph->AddLinkVertex("Features","F");
	//tableToGraph->AddLinkVertex("Arbor_Complexity","A");
	tableToGraph->AddLinkEdge("Features", "Morphology");
	tableToGraph->AddLinkEdge("Features", "Semantic");
	tableToGraph->AddLinkEdge("Morphology", "Cell");
	//tableToGraph->Update();

	vtkSmartPointer<vtkGraphLayoutView> graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
	graphLayoutView->AddRepresentationFromInput(tableToGraph->GetOutput());
	graphLayoutView->SetVertexLabelArrayName("label");
	graphLayoutView->SetVertexLabelFontSize(10);
	graphLayoutView->SetVertexLabelVisibility(true);
	graphLayoutView->SetEdgeLabelArrayName("value");
	graphLayoutView->SetEdgeLabelFontSize(6);
	graphLayoutView->SetEdgeLabelVisibility(true);
	graphLayoutView->ColorVerticesOn();
	//graphLayoutView->SetLayoutStrategyToRandom();
		//SetLayoutStrategyToCircular();
//	graphLayoutView->SetLayoutStrategyToClustering2D();
	//graphLayoutView->SetLayoutStrategyToTree();
	//graphLayoutView->SetLayoutStrategyToCone(); //invalid
	graphLayoutView->SetLayoutStrategyToSpanTree();
	graphLayoutView->Update();
	
	mainQTRenderWidget->SetRenderWindow(graphLayoutView->GetRenderWindow());
	mainQTRenderWidget->resize(1200, 1200);
	mainQTRenderWidget->show();
	graphLayoutView->GetRenderWindow()->SetSize(1200, 1200);
	graphLayoutView->ResetCamera();
	graphLayoutView->GetRenderWindow()->Render();
	graphLayoutView->GetRenderWindow()->GetInteractor()->Start();

	//graphLayoutView->GetRenderWindow()->SetSize(600,600);
	//graphLayoutView->ResetCamera();
	//graphLayoutView->Render();

}