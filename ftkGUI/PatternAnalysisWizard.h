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
#include <QtGui/QFileDialog>
#include <QtGui/QDialog>

//VTK INCLUDES
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>
#include <vtkVariantArray.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <float.h>

#include <PatternAnalysis/libsvm/svm.h>
#ifdef USE_KPLS
#include <PatternAnalysis/embrex/kpls.h>
#endif


class PatternAnalysisWizard : public QWizard
{
    Q_OBJECT;

public:
	enum { Page_Start, Page_Features, Page_Training, Page_Parameters, Page_Execute };
	typedef enum { _SVM, _KPLS, _SEGMODEL, _APPENDMODEL, _ACTIVE, _ACTIVEMODEL, _CLUS } Module;

	PatternAnalysisWizard(vtkSmartPointer<vtkTable> table, Module mod,
                        const char * trainColumn, const char * resultColumn,
                        QWidget *parent = 0);
	PatternAnalysisWizard(vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> model_table, QString fileName, Module mod,
                        const char * trainColumn, const char * resultColumn,
                        QWidget *parent = 0);
	void KPLSrun(std::vector<int> columnsToUse);
	vtkSmartPointer<vtkTable> mod_table;
	QString filename;
	inline vtkSmartPointer<vtkTable>  getExtractedTable(){ return new_table;} 
	bool extractedTable;
	bool fromModel;
	

protected:
	//void initializePage(int id);
	//void cleanupPage(int id);
	bool validateCurrentPage();
	int nextId() const;

	void runSVM();
	void runKPLS();
	void extractTable(bool flag = false);
	void extractTableFromModel(vtkSmartPointer<vtkTable>);
	void saveModel(void);
	void appendModel(vtkSmartPointer<vtkTable>, QString);


signals:
	void changedTable(void);
	void enableModels(void);
	void start_training(vtkSmartPointer<vtkTable>);

private:
	void initFeatureGroup(void);
	void disabledFeatureGroup(void);
	QButtonGroup *featureGroup;

	void initOptionGroup(void);
	QButtonGroup *optionGroup;
	QString lastPath;

	vtkSmartPointer<vtkTable> m_table;
	const char * columnForTraining;
	const char * columnForPrediction;
	Module m_module;
	vtkSmartPointer<vtkTable> new_table;

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
	FeaturesPage(QButtonGroup *fGroup, int caller = 0, QWidget *parent = 0);
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
