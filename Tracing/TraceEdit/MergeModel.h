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
#include <ftkGUI/ObjectSelection.h>

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

	int RowForID(int id);
	int GetNumFeatures();
	int GetNumGaps();
	void SetTraceGaps(std::vector<TraceGap*> gaps);
	std::vector<TraceGap*>GetTraceGaps();
	std::vector<int> GetSelectedGapIDs();
	void SelectbyTraceID(int id);
	vtkSmartPointer<vtkTable> getDataTable();
	ObjectSelection * GetObjectSelection();

private:
  //the first column of each row is the gap ID, as per the definition
  //of this data model's headers.
	int NumFeatures;
	int NumGaps;
  std::vector<TraceGap*> TraceGaps;
	void SetupHeaders();
	void SyncModel();	
	vtkSmartPointer<vtkTable> DataTable;
	ObjectSelection * Selection;
};
#endif
