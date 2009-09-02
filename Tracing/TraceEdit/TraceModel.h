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

class QStandardItemModel;
class QItemSelectionModel;
class TraceLine;

class TraceModel : public QObject
{
	Q_OBJECT

public:
	TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders);
	TraceModel(std::vector<TraceLine*> trace_lines);

	QStandardItemModel *GetModel()
	{
		return this->Model;
	};
	QItemSelectionModel *GetSelectionModel()
	{
		return this->SelectionModel;
	};
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
	std::vector<int> GetSelecectedIDs();
	void SelectByIDs(int ID);

signals:
	//emit this signal to tell the Qt views to update
	void modelChanged(void);
	void selectionChanged(void);

public slots:
	//void deleteTrigger(void);
private:	
	const static int IDColumn = 0;
	QMap<int, int> IDToRowMap;
	std::vector<TraceLine*> TraceLines;
	std::vector<QString> headers;
	int NumTraces;
	int NumFeatures;
	QStandardItemModel *Model;
	QItemSelectionModel *SelectionModel; 
	void SetupHeaders();
	void SyncModel();	
	void MapTracesToRows();
};
#endif