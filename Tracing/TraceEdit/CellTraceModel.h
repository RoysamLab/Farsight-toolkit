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
#ifndef CELLTRACEMODEL_H
#define CELLTRACEMODEL_H

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
class QStandardItemModel;
class QItemSelectionModel;
class TraceLine;
class CellTrace;
class CellTraceModel : public QObject
{
	Q_OBJECT

public:
	CellTraceModel();
	CellTraceModel(std::vector<CellTrace*> Cells);
	void setCells(std::vector<CellTrace*> Cells);
	vtkSmartPointer<vtkTable> getDataTable();
	ObjectSelection * GetObjectSelection();
	void SelectByRootTrace(std::vector<TraceLine*> roots);
	void SelectByIDs(std::vector<int> IDs);
	std::vector<int> GetSelecectedIDs();
signals:
	void selectionChanged(void);
private:
	std::vector<CellTrace*> Cells;
	std::vector<QString> headers;
	void SetupHeaders();
	void SyncModel();
	vtkSmartPointer<vtkTable> DataTable;
	ObjectSelection * Selection;
};

#endif