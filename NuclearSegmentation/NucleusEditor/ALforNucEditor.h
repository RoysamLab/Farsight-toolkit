
#ifndef AL_NUC_EDITOR_H
#define AL_NUC_EDITOR_H

#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QToolBar>
#include <QtGui/QProgressBar>
#include <QtGui/QDialog>
#include <QtGui/QVBoxLayout>
#include <QtGui/QRadioButton>
#include <QtCore/QFileInfo>
#include <QtCore/QThread>
#include <QtCore/QSettings>
#include <QtCore/QSignalMapper>
#include <QtGui/QGridLayout>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>

#include "ActiveLearningDialog.h"
#include "ftkCommon/ftkUtils.h"
#include "ftkGUI/TrainingDialog.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/LabelImageViewQT.h"
#include "PatternAnalysis/activeLearning/mclr.h"

//#ifdef USE_Clusclus
#include "ClusClus/HeatmapWindow.h"
#include "ClusClus/Heatmap.h"
//#endif

#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariant.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>


class ALforNucEd : public QWidget
{
	Q_OBJECT;

public:
	ALforNucEd( bool val=false );
	~ALforNucEd();

	inline void SetTablesToClassify(std::vector< vtkSmartPointer< vtkTable > > tbls){ classificationTables = tbls; };
	inline void SetTableForTraining(vtkSmartPointer< vtkTable > tbl){ trainingTable = tbl; };
	void SetLabelView(LabelImageViewQT *view);
	void RunALClassification(bool from_model);
	inline std::vector< vtkSmartPointer< vtkTable > > GetClassificationResult(){ return classificationTables; };
	inline std::vector< std::pair< std::string, vnl_vector<double> > > GetALModel(){ return activeModel; };

signals:
	void Classification_Done();

private slots:
	void Start_Training(vtkSmartPointer<vtkTable> pTable);
	void ALDialogPopUP(bool first_pop, std::vector<std::pair<int,int> > query_labels);
	//void CreateActiveLearningModel(MCLR* mclr_alm,  vtkSmartPointer<vtkTable> pWizard_table);
	void Start_Classification(bool create_model = false);
	void Perform_Classification();

private:
	std::vector< vtkSmartPointer< vtkTable > > classificationTables;
	vtkSmartPointer< vtkTable > trainingTable;
	LabelImageViewQT *labelView;
	PatternAnalysisWizard *pWizard;
	bool pixel_class;

	#ifdef	USE_Clusclus
		Heatmap *HeatmapWin;		
	#endif
	
	MCLR *mclr;
	ActiveLearningDialog *alDialog;
	double confidence_thresh;
	bool from_model;
	std::string classification_name;
	std::vector< std::pair<int,int> > id_time;
	vtkSmartPointer<vtkTable> pawTable;
	std::vector<int> active_queries;
	std::vector<QImage> snapshots;
	vnl_vector<double> std_dev_vec;
	vnl_vector<double> mean_vec; 
	vnl_matrix<double> act_learn_matrix;
	std::vector< std::pair< std::string, vnl_vector<double> > >activeModel;
	std::string prediction_col_name;
	std::string confidence_col_name;

};


class ClassNameConfidenceDialog : public QDialog
{
	Q_OBJECT
public:
	ClassNameConfidenceDialog(QWidget *parent = 0);
	std::string getClassName();
	double getConfThresh();

private:
	QLabel *classNameLabel;
	QLineEdit *class_name;
	QHBoxLayout *classNameLayout;
	QLabel *confidenceLabel;
	QLineEdit *conf_thresh;
	QHBoxLayout *confLayout;
	QPushButton *okButton;
	QHBoxLayout *bLayout;
	QVBoxLayout *layout;
};

#endif
