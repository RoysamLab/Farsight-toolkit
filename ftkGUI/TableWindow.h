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
#include <QtGui/QDoubleSpinBox>


#include <iostream>

//#include "SegmentationModel.h"

class ChooseItemDialog;
class ChooseItemsDialog;
class FilterRowsDialog;

class TableWindow : public QMainWindow
{
    Q_OBJECT;

public:
	TableWindow(QItemSelectionModel *mod, QWidget *parent = 0);

	//void SetModels(QItemSelectionModel *selectionModel);
	void ResizeToOptimalSize(void);

signals:
	void closing(QWidget *widget);
	void sorted();

protected:
	void closeEvent(QCloseEvent *event);

public slots:
	void update();
	void modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight);

private slots:
	void createMenus();
	void sortBy();
	void changeColumns();
	void showFilters();
    
private:
	QTableView *table;
	FilterRowsDialog *filters;

	QMenu *viewMenu;
	QAction *visibleColumnsAction;
	QAction *sortByAction;
	QAction *filterRowsAction;
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
	FilterRowsDialog(QTableView *table, QWidget *parent = 0);

private:
	QTableView *mTable;
	QGridLayout *fLayout;
	QGroupBox *groupBox;
	QPushButton *addButton;
	QPushButton *delButton;
	QDoubleSpinBox *minVal1;
	QDoubleSpinBox *minVal2;
	QDoubleSpinBox *minVal3;
	QDoubleSpinBox *maxVal1;
	QDoubleSpinBox *maxVal2;
	QDoubleSpinBox *maxVal3;
	QPushButton *minComp1;
	QPushButton *minComp2;
	QPushButton *minComp3;
	QPushButton *maxComp1;
	QPushButton *maxComp2;
	QPushButton *maxComp3;
	QComboBox *feature1;
	QComboBox *feature2;
	QComboBox *feature3;
	QComboBox *bool1;
	QComboBox *bool2;
	QPushButton *updateButton;

	int numEquations;
	
	QString smaller;
	QString bigger;

	QStringList GetVisibleFeatures();
	QComboBox * NewFeatureCombo();
	QComboBox * NewBoolCombo();
	QPushButton * NewCompButton(int n);
	int GetColumnFor(QString headerText);
	void GetMinMaxFor(int c, double *min, double *max);

private slots:
	void DoFilter();
	void AddEquation();
	void RemoveEquation();
	void RemoveWidget(QWidget *widget);
	void AddWidget(QWidget *widget, int r, int c);
	void ToggleComp1(){ ToggleComp(1); };
	void ToggleComp2(){ ToggleComp(2); };
	void ToggleComp3(){ ToggleComp(3); };
	void ToggleComp(int n);
	void InitRanges();
	void SetF1Ranges(QString text);
	void SetF2Ranges(QString text);
	void SetF3Ranges(QString text);
};

#endif
