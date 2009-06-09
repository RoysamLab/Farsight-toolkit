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

#ifndef RESULTLOADWIDGET_H
#define RESULTLOADWIDGET_H

#include "FarsightConfig.h"

#include <iostream>
#include <fstream>

//QT INCLUDES
#include <QtGui/QWidget>
#include <QtGui/QFileDialog>
#include <QtCore/QString>

//GUI INCLUDES
#include "PlotWindow.h"
#include "TableWindow.h"
#include "SegmentationWindow.h"

//OTHER LOCAL INCLUDES
#include "ftkImage/ftkImage.h"
#include "tinyxml/tinyxml.h"


//***************************************************************************
// THIS CLASS LOADS DIRECTLY FROM AN XML FILE TO DISPLAY RESULTS
// AND ALLOW FOR RESULTS NAVIGATION AND VIEWING.
// EDITING IS NOT ALLOWED BECAUSE THERE IS NO CONNECTION TO THE 
// MODULE/PACKAGE THAT CREATED THESE RESULTS.  
//***************************************************************************
class ResultLoadWidget : public QWidget
{
    Q_OBJECT;

public:
    ResultLoadWidget(FTKItemModel*,QItemSelectionModel*);
	void CreateNewPlotWindow();

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void closeWidget(QWidget *widget);

private:
	SegmentationWindow *segWin;
	TableWindow *tblWin;
	std::vector<PlotWindow *> pltWin;

	FTKItemModel *model;				  //comes from parent
	QItemSelectionModel *selectionModel;  //comes from parent

	ftkImage *dataImg;
	ftkImage *labelImg;

	QString path;
	QString filename;
	QString datafname;
	QString labelfname;

	void loadFromResultImages();
	void loadFromXML();
	void loadXML();
	void loadOutliers();
	void extractPath();
	//void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewSegmentationWindow();
};


#endif