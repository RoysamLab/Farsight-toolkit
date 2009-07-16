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

class TraceGap;

class MergeModel : public QStandardItemModel
{
	Q_OBJECT

public:
	MergeModel();
  MergeModel(std::vector<TraceGap*> gaps);
	~MergeModel();

	//QStandardItemModel *GetModel();
	//QItemSelectionModel *GetSelectionModel();
  //get the column that contains the gap id for each row
	int ColumnForID();
	int RowForID(int id);
	int GetNumFeatures();
	int GetNumGaps();
  void SetTraceGaps(std::vector<TraceGap *> gaps);
  std::vector<TraceGap *>GetTraceGaps();

signals:
	void modelChanged(void);

public slots:
	void deleteTrigger(void);
	void mergeTrigger(void);

private:
  //the first column of each row is the gap ID, as per the definition
  //of this data model's headers.
	const static int IDColumn = 0;
	int NumFeatures;
	int NumGaps;
  std::vector<TraceGap *> TraceGaps;
	QMap<int, int> IDToRowMap;	

	//QStandardItemModel *Model; //not needed since this class is *JUST* a model, right?
	//QItemSelectionModel *SelectionModel; 
    //not needed since this class inherits it from StandardItemModel....
  void SetupHeaders();
	void SyncModel();
	void UpdateMapping();
};
#endif
