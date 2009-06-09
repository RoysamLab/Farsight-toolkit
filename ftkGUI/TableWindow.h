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
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QCloseEvent>

#include <iostream>

//#include "SegmentationModel.h"

class TableWindow : public QWidget
{
    Q_OBJECT;

public:
	TableWindow(QItemSelectionModel *mod, QWidget *parent = 0);

	//void SetModels(QItemSelectionModel *selectionModel);
	void ResizeToOptimalSize(void);

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

public slots:
	void update();
	void modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight);
    
private:
	QVBoxLayout *layout;
	QTableView *table;

	int visibleRows;

	//void setup(void);
};

#endif
