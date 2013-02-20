#ifndef NODEMODEL_H
#define NODEMODEL_H

#include <vector>
#include <string>
//QT INCLUDES
#include <QtCore>
#include <QtGui>
//vtk table includes
#include "vtkTable.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <ftkGUI/ObjectSelection.h>

//selective clustering to replace object selection
#include "SelectiveClustering.h"
#include "QvtkTableView.h"
#include "SelectionUtilities.h"

// 
#include <map>

class QStandardItemModel;
class QItemSelectionModel;
class TraceBit;
class TraceLine;

class NodeModel : public QObject
{
	Q_OBJECT

public:
	NodeModel();
	NodeModel(std::vector<TraceLine*> trace_lines);
	~NodeModel();

	void SetLines(std::vector<TraceLine*> trace_lines);
	//std::vector<TraceBit> GetNodes()
	//{
	//	return this->TraceBits;
	//};
	vtkSmartPointer<vtkTable> getDataTable();
	void setDataTable(vtkSmartPointer<vtkTable> table);
	ObjectSelection * GetObjectSelection();
	void SelectByIDs(std::vector<int> IDs);

	void AddNodeHeader(std::string NewFeatureHeader);

	unsigned int getNodeCount();
	//void SelectByIDs(std::vector<int> IDs); //already used elsewhere
	std::vector<long int> GetSelectedIDs();
	//std::vector<TraceBit> GetSelectedNodes();
	//std::map< int ,TraceBit>::iterator GetNodeiterator();
	//std::map< int ,TraceBit>::iterator GetNodeiteratorEnd();


signals:
	//emit this signal to tell the Qt views to update
	void selectionChanged(void);

private:
	void SetupHeaders();
	void SyncModel();

	std::vector<TraceLine*> TraceLines;
	std::vector<TraceBit> TraceBits;
	std::vector<QString> headers;
	std::vector<std::string> additionalHeaders;

	vtkSmartPointer<vtkTable> DataTable;

	ObjectSelection * Selection;

	//std::map< int ,TraceBit>::iterator NodeIDLookupIter;
	//std::map<long int ,TraceBit> NodeIDLookupMAP;
	std::vector<TraceLine*>::iterator TraceIDLookupIter;
};
#endif
