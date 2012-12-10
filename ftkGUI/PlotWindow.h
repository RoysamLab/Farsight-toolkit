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
#include <QtGui/QVBoxLayout>
#include <QtGui/QCloseEvent>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QActionGroup>

#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include "ScatterView.h"
#include "ObjectSelection.h"

class PlotWindow : public QMainWindow
{
    Q_OBJECT;

public:
	PlotWindow(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels = NULL);

public slots:
	void update(void);

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void xChange(QAction *action);
	void yChange(QAction *action);
	void colorChange(QAction *action);
    
private:
	QMenu *optionsMenu;
	QAction *normalizeAction;	//Normalize the data in the plot
	QMenu *xMenu;				//Choose x column
	QMenu *yMenu;				//Choose y column
	QMenu *colorMenu;			//Choose color column
	QAction *clearAction;		//clear selections
	QAction *windowAction;		//set window size and interval

	ScatterView *scatter;

	void setupUI(void);	//for initial setup
	void updateOptionMenus(vtkSmartPointer<vtkTable> tbl);
 };

#endif
