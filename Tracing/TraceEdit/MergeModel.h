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
#ifndef MERGEMODEL_H
#define MERGEMODEL_H

#include <vector>

//QT INCLUDES
#include <QtCore>
#include <QtGui>
#include "vtkTable.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <ftkGui/ObjectSelection.h>

class QStandardItemModel;
class QItemSelectionModel;

class TraceGap;

class MergeModel : public QObject
{
	Q_OBJECT

public:
	MergeModel();
  MergeModel(std::vector<TraceGap*> gaps);
	~MergeModel();

	QStandardItemModel *GetModel();
	QItemSelectionModel *GetSelectionModel();

	int RowForID(int id);
	int GetNumFeatures();
	int GetNumGaps();
  void SetTraceGaps(std::vector<TraceGap*> gaps);
  std::vector<TraceGap*>GetTraceGaps();
  std::vector<int> GetSelectedGapIDs();
  void SelectbyTraceID(int id);
	vtkSmartPointer<vtkTable> getDataTable();
	ObjectSelection * GetObjectSelection();

signals:
  //emit this signal to tell the Qt views to update
  void modelChanged(void);
  void selectionChanged(void);

public slots:
	void deleteTrigger(void);
	void mergeTrigger(void);
	void MapGapIDsToRows();
	void MapTracesToRows();

private:
  //the first column of each row is the gap ID, as per the definition
  //of this data model's headers.
	const static int IDColumn = 0;
	const static int Trace1Col = 1;
	const static int Trace2Col = 2;
	int NumFeatures;
	int NumGaps;
  std::vector<TraceGap*> TraceGaps;
	QMap<int, int> IDToRowMap;
	QMap<int, int> Trace1ToRow;
	QMap<int, int> Trace2ToRow;

	QStandardItemModel *Model;
	QItemSelectionModel *SelectionModel; 
	void SetupHeaders();
	void SyncModel();	
	vtkSmartPointer<vtkTable> DataTable;
	ObjectSelection * Selection;
};
#endif
