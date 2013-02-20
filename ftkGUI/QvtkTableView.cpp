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

	this->AnnotationLink = 0;

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
	if(this->AnnotationLink)
	{
		this->AnnotationLink->RemoveAllInputs();
		this->AnnotationLink->RemoveAllObservers();
	}
	this->TableView->~QTableView();
	delete this->TableAdapter;
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
	vtkSmartPointer<vtkIdTypeArray> Selected = this->getSelectedObjects();

	vtkSelection * TableRowSelection = SelectionUtilities::ConvertIDsToVTKSelection(Selected);

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
	vtkSmartPointer<vtkSelection> TableRowSelection = vtkSmartPointer<vtkSelection>::New();
	vtkSmartPointer<vtkIdTypeArray> TableRowIDs = vtkSmartPointer<vtkIdTypeArray>::New();
	std::map< vtkIdType, vtkIdType>::iterator  idIter;
	vtkIdTypeArray* vertexList = SelectionUtilities::ConvertVTKSelectionToIDArray(selection);

	if( vertexList != NULL)
	{
		vtkIdType numTuples = vertexList->GetNumberOfTuples();

		vtkSmartPointer<vtkIdTypeArray> ConvertedVertexList = vtkSmartPointer<vtkIdTypeArray>::New();
		ConvertedVertexList->SetNumberOfComponents(1);

		if( numTuples > 0)
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
		TableRowSelection = SelectionUtilities::ConvertIDsToVTKSelection(TableRowIDs);
	}//end vertices null

	QvtkView->setCurrentVTKSelection(TableRowSelection);
}

vtkSmartPointer<vtkIdTypeArray> QvtkTableView::getSelectedObjects()
{
	/*!
	* returns the array of selected pedigree ids
	*/
	
	vtkSmartPointer<vtkIdTypeArray> ItemSelections = vtkSmartPointer<vtkIdTypeArray>::New();
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
	/*!
	* Constructor to initalize 
	*/
	TableView = new QvtkTableView();
	this->closed = false; 
	this->setModal(false);
	this->setAttribute( Qt::WA_DeleteOnClose, true );
	QVBoxLayout * MainVBoxLayout = new QVBoxLayout();
	MainVBoxLayout->addWidget(TableView);
	this->setLayout(MainVBoxLayout);
}
QvtkTableDialog::~QvtkTableDialog()
{
	/*!
	* Delete
	*/
}

void QvtkTableDialog::setTitle(std::string title)
{
	/*!
	* Sets Window Title
	*/
	this->setWindowTitle( QString(title.c_str()));
}

void QvtkTableDialog::UpdateView(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink)
{
	/*!
	* Pass inputs into tableView and updates visibility
	*/
	this->TableView->SetInputLink(InputTable, InputAnnotationLink);
	this->setVisible(true);
}

void QvtkTableDialog::close()
{
	/*!
	* Program closed  
	*/
	this->setVisible(0);
	this->TableView->~QvtkTableView();
}

void QvtkTableDialog::closeEvent(QCloseEvent *event)
{
	/*!
	* User Closed 
	*/
	this->closed = true;
	this->close();
	event->accept();
}
