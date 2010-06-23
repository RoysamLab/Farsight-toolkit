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
#ifndef STATISTICSTOOLBAR_H
#define STATISTICSTOOLBAR_H
//VTK Includes:
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkVariantArray.h>
#include <vtkDoubleArray.h>

//QT
#include<QtGui>
#include<QtCore>
//STD 
#include <vector>
#include <string>



class QStandardItemModel;
class QItemSelectionModel;

class StatisticsToolbar : public QWidget
{
	Q_OBJECT
public:
	StatisticsToolbar(QWidget *parent );
	~StatisticsToolbar();
	void setTable(vtkSmartPointer<vtkTable> dataTable);
    QDockWidget *statisticsDockWidget;

private:
	vtkSmartPointer<vtkTable> inputDataTable;
	
	QStandardItemModel * StatisticsModel;
	QItemSelectionModel * StatisticsSelectionModel;
	std::vector<QString> colHeaders;
	std::vector<QString> rowHeaders;
	std::vector< std::vector<double>> statistics;
	QList<double> *myList;
	  QTableView * StatisticsTable;
	
    
	double Average(vtkAbstractArray *Column);
	double StDeviation(vtkAbstractArray *Column, double average);
	QList<double> SortColumn(vtkAbstractArray *Column);
	void SetUpHeaders(vtkSmartPointer<vtkTable> dataTable);
	std::vector<double> ComputeStatistics(vtkAbstractArray *Column);
	double Mode(QList<double> dataList);
	
	
};
#endif