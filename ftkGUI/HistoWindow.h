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

#ifndef HISTOWINDOW_H
#define HISTOWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QComboBox>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QCloseEvent>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QActionGroup>
#include <QtCore/QAbstractItemModel>

#include <vtkTable.h>
#include <vtkSmartPointer.h>
//#include <vtkQtBarChartView.h>
#include <vtkDoubleArray.h>
#include <vtkCallbackCommand.h>
#include "ObjectSelection.h"

#include <set>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

#if (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 5)
  #define CHART_IS_WIDGET 1
#else
	#define CHART_IS_WIDGET 0 
#endif

class HistoWindow : public QMainWindow
{
	Q_OBJECT;

public:
	HistoWindow(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL, std::vector< std::vector< long int> > * clusIndex = NULL);

public slots:
	void update(void);
signals:
	void closing(QWidget *widget);

public:
	//vtkSmartPointer<vtkQtBarChartView> chartView;

protected:
	void closeEvent(QCloseEvent *event);

private:
	void SyncModel();
	bool ReadHistogramData(const char* fileName);
	void Normalize();
	void SetNumofBins(int n);
	void SetClusterNumber(int n);
	void SetLogNumber(int n);
	bool findFrequencies();
	void setBucketNames();
	void ConstructBarChart();
	//Helpers:
	void updateOptionMenus();
	void initMap(std::map<int,double> &v, int n);
	double GetMaxNumber(){ std::multiset<double>::iterator pos=std::max_element(data.begin(), data.end());
						   return *pos;};
private slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	void columnChange(QAction *action);
	void binsChange(QAction *action);
	void clusterChange(QAction *action);
	void logChange(QAction *action);
	void saveFrequency(QAction *action);

private:
	int columnNum;
	int clusterNum;
	int logNum;
	std::string columnName;
	int numofbins;				//
	bool status;				//will be false if num of bins is less than 2 or larger than 10	
	bool normalized;			//false if not normalized
	double distanceToUpperBound;//change this if you need more precision for bar charts and frequencies
	std::map<int,double> result_fq;
	std::map<int,double> cluster_result_fq; // of the values in each cluster
	std::vector<std::string> names; //Names of bins will be stored here

	std::multiset<double> data;		//The column of all data
	std::multiset<double> cluster_data;

	std::vector< std::vector< long int> > clusIndex;

	ObjectSelection *selection;
	vtkSmartPointer<vtkTable> m_table;
	vtkSmartPointer<vtkTable> hisTable;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;

	QVBoxLayout *layout;

	QMenu *optionsMenu;
	QMenu *columnMenu;
	QMenu *binsMenu;
	QMenu *clusterNoMenu;
	QMenu *logMenu;
	QMenu *fileMenu;
};

#endif
