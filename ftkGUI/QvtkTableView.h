
#ifndef QVTKTABLEVIEW_H
#define QVTKTABLEVIEW_H

//#include "vtkQtView.h"
#include <QSortFilterProxyModel>

#include <QHeaderView>
#include <QItemSelection>
#include <QTableView>

#include <QtGui/QVBoxLayout>

#include <QtGui/QCloseEvent>
#include <QtGui/QDialog> 
#include <QString>

#include "vtkAbstractArray.h"
#include "vtkAddMembershipArray.h"
#include "vtkAlgorithm.h"
#include "vtkAlgorithmOutput.h"
#include "vtkAnnotation.h"
#include "vtkAnnotationLayers.h"
#include "vtkAnnotationLink.h"
#include "vtkApplyColors.h"
#include "vtkConvertSelection.h"
#include "vtkDataObjectToTable.h"
#include "vtkDataRepresentation.h"
#include "vtkDataSetAttributes.h"
#include "vtkGraph.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkOutEdgeIterator.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkSmartPointer.h"
#include "vtkTable.h"
#include "vtkViewTheme.h"


#include "vtkQtTableModelAdapter.h"
#include "vtkQtAbstractModelAdapter.h"
#include "vtkSmartPointer.h"
#include "vtkCallbackCommand.h"

#include <string>
#include <map>

#include "SelectionUtilities.h"

class QvtkTableView : public QTableView
{
Q_OBJECT
public:
	QvtkTableView();
	~QvtkTableView();
	//QvtkTableView(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink);

	void SetInputLink(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink);

	vtkSmartPointer<vtkIdTypeArray> getSelectedObjects();
	void setCurrentVTKSelection(vtkSelection * TableRowSelection);

private slots:

	void slotQtSelectionChanged(const QItemSelection&,const QItemSelection&);

private:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );


private:
	QTableView * TableView;
	vtkQtTableModelAdapter* TableAdapter;
	QSortFilterProxyModel* TableSorter;

	std::map< vtkIdType, vtkIdType> IdLookUP; 

	//VTK data and selection
	vtkSmartPointer<vtkTable> DataTable;
	vtkSmartPointer<vtkAnnotationLink> AnnotationLink;
	bool QTSelectionBool;
	
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;

};


class  QvtkTableDialog: public QDialog
{
	Q_OBJECT
public:
	QvtkTableDialog();
	~QvtkTableDialog();

	void setTitle(std::string title);
	void UpdateView(vtkSmartPointer<vtkTable> InputTable, vtkSmartPointer<vtkAnnotationLink> InputAnnotationLink);

	void close();
	void closeEvent(QCloseEvent *event);
	
	QvtkTableView * TableView;

	bool closed;
};

#endif
