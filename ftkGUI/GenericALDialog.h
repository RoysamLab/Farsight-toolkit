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
#ifndef GENERICALDIALOG_H
#define GENERICALDIALOG_H

#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QDialog>
#include <QtGui/QColorDialog>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtCore/QStringList>
#include <QtCore/QSignalMapper>
#include <QtGui/QTableWidget>
#include <QtGui/QTableWidgetItem>
#include <QtGui/QLineEdit>
#include <QtGui/QPlainTextEdit>
#include <QHeaderView>
 #include <QRadioButton>

#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkVariant.h>

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

class GenericALDialog : public QDialog
{
	Q_OBJECT
public:
	GenericALDialog(vtkSmartPointer<vtkTable> table,int num_classes,int active_query,std::vector<int>top_feats,QWidget *parent=0);
	bool finish;
	//QLineEdit *lineEdit;
	int class_selected;
	std::vector<QRadioButton *> button_vector;
	//std::pair<QImage,int> temp_pair;

private slots:
	void finished();
	void Set_Class();
	
private:
	//QMap<QString, QColor> m_snapshot;
	//QSignalMapper *signalMapper;
	QPushButton *nextButton;
	QPushButton *doneButton;
};


#endif
