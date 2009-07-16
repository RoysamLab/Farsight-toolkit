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

//QT INCLUDES
#include <QtCore>
#include <QtGui>
#include "TraceObject.h"

class mergeModel : public QObject
{
	Q_OBJECT

public:
	mergeModel();
	~mergeModel();

	QStandardItemModel *GetModel(){ return model; };
	QItemSelectionModel *GetSelectionModel();

	int ColumnForID()
	{ 
		return IDColumn; 
	};
	int RowForID(int id);
	int NumFeatures()
	{ 
		return numFeatures; 
	};
	int NumGaps()
	{ 
		return numGaps; 
	};
signals:
	void modelChanged(void);

public slots:
	void deleteTrigger(void);
	void mergeTrigger(void);

private:
	int IDColumn;
	int numFeatures;
	int numGaps;
	
	QMap<int, int> LabelToRowMap;	
	QStandardItemModel *model;
	QItemSelectionModel *selectionModel;
	void SyncModel();
	void updateMapping();
};
#endif