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
#ifndef TRAININGDIALOG_H
#define TRAININGDIALOG_H

//QT INCLUDES
#include <QtGui/QLabel>
#include <QtGui/QDialog>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QVBoxLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QCheckBox>
#include <QtGui/QScrollArea>
#include <QtGui/QMessageBox>
#include <QtGui/QLineEdit>
#include <QtGui/QFileDialog>

//VTK INCLUDES
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <float.h>

#include "ftkUtils.h"
#include "ftkGUI/LabelImageViewQT.h"

class TrainingDialog : public QDialog
{
    Q_OBJECT;

public:
	TrainingDialog(vtkSmartPointer<vtkTable> table, const char * trainColumn,std::string mode,int firstTableSize, QWidget *parent=0);
	//TrainingDialog(vtkSmartPointer<vtkTable> table, const char * trainColumn,std::string mode, QWidget *parent = 0);
	void loadModelFromFile(std::string file_name);
	int firstTableRows;

public slots:
	  void accept();
	  
protected:

signals:
	void changedTable(void);

private slots:
	void addClass(void);
	void remClass(void);
	void saveModel(void);
	void loadModel(void);
	void updateTable(void);
	void updateInputs(void);
	void parseInputValues(void);
	void parseTableValues(void);
	void Append(void);
	void tableToInput(void);
	void inputToTable(void);
	void class_1_selected(void);
	void class_2_selected(void);
	void class_3_selected(void);
	void class_4_selected(void);
	void class_5_selected(void);
	void class_6_selected(void);

private:
	QVBoxLayout * inputsLayout;
	QVBoxLayout * classesLayout;
	QVector<QHBoxLayout *> iLayouts;
	QVector<QLabel *> inputLabels;
	QVector<QLineEdit *> inputValues;
	QLineEdit *className;
	QPushButton * addButton;
	QPushButton * delButton;
	QPushButton * saveButton;
	QPushButton * loadButton;
	QPushButton * quitButton;
	QPushButton * doneButton;
	QRadioButton *radio1;
	QRadioButton *radio2;
	QRadioButton *radio3;
	QRadioButton *radio4;
	QRadioButton *radio5;
	QRadioButton *radio6;
	void GetTrainingNames(vtkSmartPointer<vtkTable> table);
	QHBoxLayout *createFirstExclusiveGroup();
	
	QString lastPath;							//Last path that has been navigated to
	//std::vector< std::set<int> > training;
	std::map< int, std::set<int> > training;
	std::vector< std::string > class_names;
	std::vector< std::string > training_names;
	vtkSmartPointer<vtkTable> m_table;
	vtkSmartPointer<vtkTable> model_table;
	const char * columnForTraining;
	void differentclassselected(int selected_class);
	bool use_train_gui;
	std::string default_training_name,new_class_name;
};

class TrainingDialogNoGUI{
public:
	TrainingDialogNoGUI(vtkSmartPointer<vtkTable> table);
	void loadModelFromFile1(std::string file_name);
private:
	void Append1(void);
	void GetTrainingNames1(vtkSmartPointer<vtkTable> table);

	std::vector< std::string > class_names;
	std::vector< std::string > training_names;
	vtkSmartPointer<vtkTable> m_table;
	vtkSmartPointer<vtkTable> model_table;
	
};
#endif
