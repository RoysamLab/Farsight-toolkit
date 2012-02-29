#include "QvtkTableView.h"

QvtkTableView::QvtkTableView()
{
	/*!
	* 
	*/	
	//this->setModal(false);
	this->TableView = new QTableView();
	this->TableAdapter = new vtkQtTableModelAdapter();

	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->DataTable->Initialize();

	this->TableSorter = new QSortFilterProxyModel();
	this->TableSorter->setSourceModel(this->TableAdapter);

	this->TableView->setModel(this->TableSorter);
	this->TableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	this->TableView->setAlternatingRowColors(true);
	this->TableView->setSortingEnabled(true);
	this->TableView->resizeColumnToContents(0);
	this->TableView->verticalHeader()->setDefaultSectionSize(25);

	this->QTSelectionBool = false;

	QObject::connect(this->TableView->selectionModel(), 
		SIGNAL(selectionChanged(const QItemSelection&,const QItemSelection&)),
		this, 
		SLOT(slotQtSelectionChanged(const QItemSelection&,const QItemSelection&)));

	this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New(); 
	this->selectionCallback->SetClientData(this);
	this->selectionCallback->SetCallback ( SelectionCallbackFunction);

	QVBoxLayout * MainLayout = new QVBoxLayout();
	MainLayout->addWidget(this->TableView);

	this->setLayout(MainLayout);
	//this->TableView->update();
}
QvtkTableView::~QvtkTableView()
{
	//
}
void QvtkTableView::SetInputLink(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink)
{
	/*!
	* 
	*/
	this->DataTable = InputTable;
	this->AnnotationLink = InputAnnotationLink;
	
	AnnotationLink->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);
	this->TableAdapter->setTable(this->DataTable);
	this->TableView->update();

}

void QvtkTableView::slotQtSelectionChanged(const QItemSelection& vtkNotUsed(s1), 
  const QItemSelection& vtkNotUsed(s2))
{
	/*!
	* 
	*/
	std::cout << "QT Selection Changed\n";
	vtkIdTypeArray * Selected = this->getSelectedObjects();

	vtkSmartPointer<vtkSelectionNode> selectNodeList = vtkSmartPointer<vtkSelectionNode>::New();
	selectNodeList->SetSelectionList( Selected );
	selectNodeList->SetFieldType( vtkSelectionNode::VERTEX );
	selectNodeList->SetContentType( vtkSelectionNode::INDICES );

	vtkSelection * TableRowSelection = vtkSelection::New();
	TableRowSelection->RemoveAllNodes();
	TableRowSelection->AddNode(selectNodeList);

	std::cout << "Annotation Link Update\n";
	this->AnnotationLink->SetCurrentSelection( TableRowSelection );
}

void QvtkTableView::SelectionCallbackFunction(vtkObject *caller, unsigned long eventId, void *clientData, void *callData)
{
	/*!
	* 
	*/
	//std::cout << "VTK Selection Changed\n";
	vtkAnnotationLink * annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection * selection = annotationLink->GetCurrentSelection();
	QvtkTableView * QvtkView = (QvtkTableView*)clientData;
	vtkSelection * TableRowSelection = vtkSelection::New();
	vtkSelectionNode* vertices = NULL;

	if( selection->GetNode(0))
	{
		if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(0);
		}
	}

	if( selection->GetNode(1))
	{
		if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(1);
		}
	}
	//QvtkView->TableView->selectionModel()->select();
}

vtkIdTypeArray * QvtkTableView::getSelectedObjects()
{
	/*!
	* 
	*/
	
	vtkIdTypeArray* ItemSelections = vtkIdTypeArray::New();
	ItemSelections->Initialize();
	ItemSelections->SetNumberOfComponents(1);

    const QModelIndexList selectedQTRows = this->TableView->selectionModel()->selectedRows();

	for(int i = 0; i < selectedQTRows.size(); ++i)
	{
		QModelIndex orig = this->TableSorter->mapToSource(selectedQTRows[i]);
		vtkIdType rowID = vtkIdType(orig.row());
		vtkIdType itemID = this->DataTable->GetValue(rowID, 0).ToTypeInt64(); // future change to pedigree col
		ItemSelections->InsertNextValue(itemID);
	}
	return ItemSelections;
}
