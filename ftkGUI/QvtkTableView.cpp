#include "QvtkTableView.h"

QvtkTableView::QvtkTableView()
{
	/*!
	* 
	*/	
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

	this->IdLookUP.clear();
	for (vtkIdType row = 0; row != this->DataTable->GetNumberOfRows(); row++)
	{
		vtkIdType value = this->DataTable->GetValue(row,0).ToTypeInt64();
		this->IdLookUP[value] = row;
	}

	AnnotationLink->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);
	this->TableAdapter->setTable(this->DataTable);
	this->TableView->update();

}

void QvtkTableView::slotQtSelectionChanged(const QItemSelection& vtkNotUsed(s1), 
  const QItemSelection& vtkNotUsed(s2))
{
	/*!
	* Update the QT selections 
	*/
	//std::cout << "QT Selection Changed\n";
	vtkIdTypeArray * Selected = this->getSelectedObjects();

	vtkSelection * TableRowSelection = this->ConvertIDsToVTKSelection(Selected);

	//std::cout << "Annotation Link Update\n";
	this->QTSelectionBool = true;
	this->AnnotationLink->SetCurrentSelection( TableRowSelection );
	this->QTSelectionBool = false;
}

void QvtkTableView::SelectionCallbackFunction(vtkObject *caller, unsigned long eventId, void *clientData, void *callData)
{
	/*!
	* 
	*/
	//std::cout << "VTK Selection Changed\n";
	QvtkTableView * QvtkView = (QvtkTableView*)clientData;
	if (QvtkView->QTSelectionBool)
	{
		return;
	}
	vtkAnnotationLink * annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection * selection = annotationLink->GetCurrentSelection();
	vtkSelection * TableRowSelection = vtkSelection::New();
	vtkIdTypeArray * TableRowIDs = vtkIdTypeArray::New();
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
	std::map< vtkIdType, vtkIdType>::iterator  idIter;

	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
		vtkIdType numTuples = vertexList->GetNumberOfTuples();

		vtkSmartPointer<vtkIdTypeArray> ConvertedVertexList = vtkSmartPointer<vtkIdTypeArray>::New();
		ConvertedVertexList->SetNumberOfComponents(1);

		if( vertexList != NULL && numTuples > 0)
		{
			//std::cout<< "number of selections: " << numTuples << "\n";
			for( vtkIdType i = 0; i < numTuples; i++)
			{
				vtkIdType value = vertexList->GetValue(i);
				//std::cout<< value << "\n";
				idIter = QvtkView->IdLookUP.find(value);
				if (idIter != QvtkView->IdLookUP.end())
				{
					TableRowIDs->InsertNextValue((*idIter).second);
				}
			}
		}//end vertex list
		TableRowSelection = QvtkView->ConvertIDsToVTKSelection(TableRowIDs);
	}//end vertices null

	QvtkView->setCurrentVTKSelection(TableRowSelection);
}

vtkIdTypeArray * QvtkTableView::getSelectedObjects()
{
	/*!
	* returns the array of selected pedigree ids
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

vtkSelection * QvtkTableView::ConvertIDsToVTKSelection(vtkIdTypeArray *vtkIDs)
{
	/*!
	* create a vtk selection from an id array
	*/
	vtkSmartPointer<vtkSelectionNode> selectNodeList = vtkSmartPointer<vtkSelectionNode>::New();
	selectNodeList->SetSelectionList( vtkIDs );
	selectNodeList->SetFieldType( vtkSelectionNode::VERTEX );
	selectNodeList->SetContentType( vtkSelectionNode::INDICES );

	vtkSelection * TableRowSelection = vtkSelection::New();
	TableRowSelection->RemoveAllNodes();
	TableRowSelection->AddNode(selectNodeList);
	return TableRowSelection;
}

void QvtkTableView::setCurrentVTKSelection(vtkSelection * TableRowSelection)
{
	/*!
	* Map the seledted rows of the origonal table to the current 
	*/
	QItemSelection qisList = this->TableAdapter->VTKIndexSelectionToQItemSelection(TableRowSelection);
	QItemSelection sortedSel = this->TableSorter->mapSelectionFromSource(qisList);

    QObject::disconnect(this->TableView->selectionModel(), 
      SIGNAL(selectionChanged(const QItemSelection&,const QItemSelection&)),
      this, SLOT(slotQtSelectionChanged(const QItemSelection&,const QItemSelection&)));
      
    this->TableView->selectionModel()->select(sortedSel, 
      QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
      
    QObject::connect(this->TableView->selectionModel(), 
     SIGNAL(selectionChanged(const QItemSelection&,const QItemSelection&)),
     this, SLOT(slotQtSelectionChanged(const QItemSelection&,const QItemSelection&)));
}

////////////////////////////////////////////////////////////////////////////////////////

QvtkTableDialog::QvtkTableDialog()
{
	TableView = new QvtkTableView();
	this->setModal(false);
	this->setAttribute( Qt::WA_DeleteOnClose, true );
	QVBoxLayout * MainVBoxLayout = new QVBoxLayout();
	MainVBoxLayout->addWidget(TableView);
	this->setLayout(MainVBoxLayout);
}
QvtkTableDialog::~QvtkTableDialog()
{
}

void QvtkTableDialog::setTitle(std::string title)
{
	this->setWindowTitle( QString(title.c_str()));
}

void QvtkTableDialog::UpdateView(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink)
{
	this->TableView->SetInputLink(InputTable, InputAnnotationLink);
	this->setVisible(true);
}

void QvtkTableDialog::close()
{
	this->setVisible(0);
	this->TableView->~QvtkTableView();
}

void QvtkTableDialog::closeEvent(QCloseEvent *event)
{
	event->accept();
}