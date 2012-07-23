#ifndef NODEMODEL_H
#define NODEMODEL_H

//QT INCLUDES
#include <QtCore>
#include <QtGui>
//vtk table includes
#include "vtkTable.h"
#include "vtkSmartPointer.h"

//selective clustering to replace object selection
#include "SelectiveClustering.h"
#include "QvtkTableView.h"
#include "SelectionUtilities.h"

// 
#include <map>

class NodeModel
{
	Q_OBJECT

public:
	NodeModel();
	~NodeModel();


private:
	void nodeHeaders();

	std::vector<QString> headers;
	vtkSmartPointer<vtkTable> DataTable;

	ObjectSelection * Selection;
};
#endif
