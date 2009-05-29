#ifndef HISTOWINDOW_H
#define HISTOWINDOW_H

#define VTK_NIGHTLY 0

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
#include <vtkQtBarChartView.h>
#include <vtkDoubleArray.h>

#include <set>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

class HistoWindow : public QMainWindow
{
	Q_OBJECT;

public:
	HistoWindow(QItemSelectionModel *mod, QWidget *parent = 0);

private slots:
	void columnChange(QAction *action);
	void binsChange(QAction *action);
	void modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight);

private:
	void SyncModel();
	bool ReadHistogramData(const char* fileName);
	void Normalize();
	void SetNumofBins(int n);
	bool findFrequencies();
	void setBucketNames();
	void ConstructBarChart();

	int columnNum;
	std::string columnName;
	int numofbins;				//
	bool status;				//will be false if num of bins is less than 2 or larger than 10	
	bool normalized;			//false if not normalized
	double distanceToUpperBound;//change this if you need more precision for bar charts and frequencies
	std::map<int,int> result_fq;
	std::vector<std::string> names; //Names of bins will be stored here

	std::multiset<double> data;		//The column of all data
	QAbstractItemModel *model;
	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkQtBarChartView> chartView;
	QVBoxLayout *layout;

	QMenu *optionsMenu;
	QMenu *columnMenu;
	QMenu *binsMenu;

	//Helpers:
	void updateOptionMenus();
	void initMap(std::map<int,int> &v, int n);
	double GetMaxNumber(){ std::multiset<double>::iterator pos=std::max_element(data.begin(), data.end());
						   return *pos;};


};

#endif