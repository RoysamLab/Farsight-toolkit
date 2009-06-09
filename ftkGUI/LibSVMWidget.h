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

#ifndef LIBSVMWIDGET_H
#define LIBSVMWIDGET_H

#include <PatternAnalysis/libsvm/svm.h>

#include <iostream>
#include <vector>

//#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QAction>
//#include <QtGui/QGridLayout>
//#include <QtGui/QComboBox>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
//#include <QtGui/QCloseEvent>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QScrollArea>
#include <QtGui/QDoubleSpinBox>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QList>


//class PlotWindow : public QWidget
class LibSVMWidget : public QWidget
{
    Q_OBJECT;

public:
	LibSVMWidget(QAbstractItemModel *mod, QWidget *parent = 0); 

private slots:
	void go();
	void selectNone();
	void selectAll();
    
private:

	QGroupBox * initFeatureBox();
	QGroupBox * initOptionBox();

	QButtonGroup *featureGroup;

	QCheckBox *normalizeBox;
	QDoubleSpinBox *nuSpin;

	QPushButton *goButton;

	QAbstractItemModel *model;

	int columnForSVM;

 };

#endif
