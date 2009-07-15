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

#ifndef SAMPLE_EDITOR_H
#define SAMPLE_EDITOR_H

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QStandardItemModel>
#include <QtGui/QTableView>
#include <QtCore/QFileInfo>


//Farsight Includes:
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/HistoWindow.h"

#include <vector>
#include <string>

class SampleEditor : public QMainWindow
{
    Q_OBJECT;

public:
	SampleEditor(QWidget * parent = 0, Qt::WindowFlags flags = 0);


private slots:
	void loadFile(void);

signals:
    
private:
	void createMenus();
	void createStatusBar();

	QMenu *fileMenu;
	QAction *loadAction;

	QStandardItemModel *model;
	QItemSelectionModel *selModel;

	QTableView *table;
	PlotWindow *plot;
	HistoWindow *histo;

 };


#endif