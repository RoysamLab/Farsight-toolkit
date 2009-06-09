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

#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

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

#include "ScatterView.h"
#include "LibSVMWidget.h"

//class PlotWindow : public QWidget
class PlotWindow : public QMainWindow
{
    Q_OBJECT;

public:
	//PlotWindow(QWidget *parent = 0);
	PlotWindow(QItemSelectionModel *mod, QWidget *parent = 0); 

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);
	//void keyPressEvent(QKeyEvent *event);

private slots:
	void xChange(QAction *action);
	void yChange(QAction *action);
	void colorChange(QAction *action);
	void modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight);
	void startSVM();
    
private:
	QMenu *optionsMenu;
	QMenu *xMenu;
	QMenu *yMenu;
	QMenu *colorMenu;

	QMenu *toolsMenu;
	QAction *svmAction;

	ScatterView *scatter;
	LibSVMWidget *svmWidget;

	void setupUI(void);	//for initial setup
	void updateOptionMenus(bool first);
 };

#endif
