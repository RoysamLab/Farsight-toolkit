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

#ifndef TABLEWINDOW_H
#define TABLEWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QCloseEvent>
#include <QtGui/QDialog>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QCheckBox>
#include <QtGui/QButtonGroup>
#include <QtGui/QScrollArea>
#include <QtGui/QScrollBar>
#include <QtGui/QDoubleSpinBox>

#include <QtCore/QMap>
#include <QtCore/QSignalMapper>

#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkSmartPointer.h>
#include <vtkQtTableModelAdapter.h>

#include "ObjectSelection.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

class ChooseItemDialog;
class ChooseItemsDialog;
class FilterRowsDialog;
class SelectionAdapter;

class TableWindow : public QMainWindow
{
    Q_OBJECT;

public:
	TableWindow(QWidget * parent = 0);
	~TableWindow();
	void setQtModels(QItemSelectionModel *mod);
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);//////////////////////
	vtkSmartPointer<vtkTable> Pointer2Table;

signals:
	void closing(QWidget *widget);
	void sorted();

protected:
	void closeEvent(QCloseEvent *event);

public slots:
	void update();
	void test();

private slots:
	void setup();
	void createMenus();
	void sortBy();
	void changeColumns();
	void exportTable();
	void showFilters();
	void selectColumns();/////////////////////////////////////////////////////////////////////
    void exportSelectedIDs();

private:
	QTableView *tableView;
	vtkQtTableModelAdapter * modAdapter;
	SelectionAdapter * selAdapter;
	ObjectSelection * selection;
	ObjectSelection * selection2;
	ObjectSelection * selectedRows;

	QMenu *viewMenu;
	QMenu *exportMenu;
	QAction *exportAction;
	QAction *visibleColumnsAction;
	QAction *sortByAction;
	QAction *filterRowsAction;
	QAction *exportIDAction;
	QAction *testAction;

	static const int rowHeight = 18;
	bool rowsSelected;
};


class ChooseItemDialog : public QDialog
{
	Q_OBJECT
public:
	ChooseItemDialog(QStringList items, QWidget *parent = 0);
	QString getSelectedItem(void){return itemCombo->currentText();};

private:
	QComboBox *itemCombo;
	QPushButton *okButton;
};

class ChooseItemsDialog : public QDialog
{
	Q_OBJECT
public:
	ChooseItemsDialog(QStringList items, QList<bool> * _selected, QWidget *parent = 0);
private slots:
	void selectionChanged(int);
	void selectNone();
	void selectAll();
private:
	QButtonGroup *itemGroup;
	QPushButton *okButton;
	QList<bool> *selected;
};


class FilterRowsDialog : public QDialog
{
	Q_OBJECT
public:
	FilterRowsDialog(QTableView *tableView, ObjectSelection * sel = 0, QWidget *parent = 0);

private:
	QTableView *mTable;
	QGridLayout *fLayout;
	QGroupBox *groupBox;
	QPushButton *addButton;
	QPushButton *delButton;
	
	static const int tests = 3;
	QDoubleSpinBox *minVal[tests];
	QDoubleSpinBox *maxVal[tests];
	QPushButton *minComp[tests];
	QPushButton *maxComp[tests];
	QComboBox *feature[tests];
	QComboBox *bools[tests-1];

	QPushButton *updateButton;

	int numEquations;
	
	QString smaller;
	QString bigger;

	QStringList GetVisibleFeatures();
	QComboBox * NewFeatureCombo(int n);
	QComboBox * NewBoolCombo();
	QPushButton * NewCompButton(int n);
	int GetColumnFor(QString headerText);
	void GetMinMaxFor(int c, double *min, double *max);

	ObjectSelection * selection;

	QSignalMapper *buttonSignalMapper;
	QSignalMapper *featurSignalMapper;

private slots:
	void DoFilter();
	void AddEquation();
	void RemoveEquation();
	void RemoveWidget(QWidget *widget);
	void AddWidget(QWidget *widget, int r, int c);
	void ToggleComp(int i);
	void InitRanges();
	void SetRanges(int i);
};

class SelectionAdapter : public QObject
{
	Q_OBJECT
public:
	SelectionAdapter();
	SelectionAdapter(QTableView *tableView);
	void SetPair(ObjectSelection * obj, QItemSelectionModel * qmod);
	ObjectSelection * getSelRows(){ return m_obj; };

signals:
	void newRowsSelected();

protected slots:
	void updateQMOD(void);
	void updateOBJ(const QItemSelection & selected, const QItemSelection & deselected);

private:
	QTableView * m_tableView;	//If I have this I will make sure selected row is visible.
	ObjectSelection * m_obj;
	QItemSelectionModel * m_qmod;
	bool okToChange;	//set to false when I am changing so I don't get into recursion
};

#endif
