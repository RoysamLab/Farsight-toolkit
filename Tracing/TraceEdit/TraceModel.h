/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#ifndef TRACEMODEL_H
#define TRACEMODEL_H

#include <vector>
#include <string>
//QT INCLUDES
#include <QtCore>
#include <QtGui>
#include "vtkTable.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <ftkGUI/ObjectSelection.h>

#include "SelectiveClustering.h"
#include "QvtkTableView.h"
#include "SelectionUtilities.h"

// 
#include <map>

class QStandardItemModel;
class QItemSelectionModel;
class TraceLine;

class TraceModel : public QObject
{
	Q_OBJECT

public:
	TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders);
	TraceModel(std::vector<TraceLine*> trace_lines);
	~TraceModel();

	int RowForID(int id);
	int GetNumFeatures()
	{
		return this->NumFeatures;
	};
	int GetNumTraces();
	void SetTraces(std::vector<TraceLine*> trace_lines);
	std::vector<TraceLine*>GetTraces()
	{
		return this->TraceLines;
	};
	void AddFeatureHeader(std::string NewFeatureHeader);
	std::vector<TraceLine*>GetSelectedTraces();
	std::vector<int> GetSelectedIDs();
	std::vector<TraceLine*> GetSelectedRoots();
	void SelectByIDs(int ID);
	void SelectByIDs(std::set<long int> ID);
	void SetSelectionByIDs(std::set<long int> ID);
	void SelectByIDs(std::vector<int> IDs);
	void SetSelectionByIDs(std::vector<int> IDs);
	vtkSmartPointer<vtkTable> getDataTable();
	ObjectSelection * GetObjectSelection();
	void CloseClusterManager();
	double scaleFactor;


signals:
	//emit this signal to tell the Qt views to update
	//void modelChanged(void);
	void selectionChanged(void);

	//void deleteTrigger(void);
private:	
	void stdHeaders();
	const static int IDColumn = 0, RootCol = 4;
	
	std::vector<TraceLine*> TraceLines;
	std::vector<QString> headers;
	int NumTraces;
	int NumFeatures;
	void SetupHeaders();
	void SyncModel();	
	vtkSmartPointer<vtkTable> DataTable;
	ObjectSelection * Selection;
	SelectiveClustering * TraceClusterSelection;
	//ClusterManager * TraceClusterManager;
	
	std::vector<std::string> additionalFeatureHeaders;
	std::map<long int ,TraceLine*> TraceIDLookupMAP;
	std::map<long int ,TraceLine*>::iterator TraceIDLookupIter;
};
#endif
