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

#ifndef PATTERNANALYSISWIZARD_H
#define PATTERNANALYSISWIZARD_H

//QT INCLUDES
#include <QtGui/QWizard>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QVBoxLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QCheckBox>
#include <QtGui/QScrollArea>
#include <QtGui/QMessageBox>

//VTK INCLUDES
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

//OTHER FARSIGHT INCLUDES
#include <PatternAnalysis/libsvm/svm.h>
#include <PatternAnalysis/embrex/kpls.h>

#include <iostream>
#include <vector>
#include <float.h>
#include <set>

class PatternAnalysisWizard : public QWizard
{
    Q_OBJECT;

public:
	enum { Page_Start, Page_Features, Page_Training, Page_Parameters, Page_Execute };
	typedef enum { _SVM, _KPLS } Module;

	PatternAnalysisWizard(vtkSmartPointer<vtkTable> table, Module mod,
                        const char * trainColumn, const char * resultColumn,
                        QWidget *parent = 0);
	void KPLSrun(std::vector<int> columnsToUse);

protected:
	//void initializePage(int id);
	//void cleanupPage(int id);
	bool validateCurrentPage();
	int nextId() const;

	void runSVM();
	void runKPLS();

signals:
	void changedTable(void);

private:
	void initFeatureGroup(void);
	QButtonGroup *featureGroup;

	void initOptionGroup(void);
	QButtonGroup *optionGroup;

	vtkSmartPointer<vtkTable> m_table;
	const char * columnForTraining;
	const char * columnForPrediction;
	Module m_module;
};


class StartPage : public QWizardPage
{
	Q_OBJECT;

public:
	StartPage(QButtonGroup *oGroup, QWidget *parent = 0);
};

class FeaturesPage : public QWizardPage
{
	Q_OBJECT;

public:
	FeaturesPage(QButtonGroup *fGroup, QWidget *parent = 0);
	bool isComplete() const;
private:
	QGroupBox * initFeatureBox(QButtonGroup *fGroup);
	QButtonGroup *featureGroup;
private slots:
	void selectNone();
	void selectAll();
};

class PatternAnalysisWizardNoGUI{
public:
	PatternAnalysisWizardNoGUI(vtkSmartPointer<vtkTable> table, const char * trainColumn, const char * resultColumn);
	void KPLSrun1(std::vector<int> columnsToUse);
private:
	vtkSmartPointer<vtkTable> m_table;
	const char * columnForTraining;
	const char * columnForPrediction;
};
/*
class TrainingPage : public QWizardPage
{
	Q_OBJECT;
public:
	TrainingPage(QWidget *parent = 0);
};
*/
/*
class ExecutePage : public QWizardPage
{
	Q_OBJECT;
public:
	ExecutePage(QWidget *parent = 0);
};
*/
#endif
