/*
 *  ftkPreprocessDialog.h
 *  Farsight
 *
 *  Created by RAGHAV on 12/3/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef FTKPREPROCESSDIALOG_H
#define FTKPREPROCESSDIALOG_H


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif


//QT Includes:
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

//Farsight Includes:
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "NuclearSegmentation/CytoplasmSegmentation/CytoplasmSegmentation.h"
#include "NuclearSegmentation/Nuclear_Association/ftkNuclearAssociationRules.h"
#include "ftkCommon/ftkLabelImageToFeatures.h"
#include "ftkCommon/ftkUtils.h"
#include "ftkImage/ftkImage.h"
//#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"

//VTK includes:
#include "vtkQtTableView.h"

// Preprocess
//#include "ftkPreprocess.h"
#include <map>


class ftkPreprocessDialog : public QDialog
{
	Q_OBJECT
public:
	ftkPreprocessDialog(QVector<QString> chs, std::string id,QWidget *parent = 0);
	QLabel *channelLabel;
	ftk::Image::Pointer myImg;
	enum Filters {Filter1,Filter2,Filter3,Filter4,Filter5,Filter6,Filter7,Filter8,Filter9,Filter10,Filter11,Filter12,Filter13,Filter14,Filter15,Filter16,Filter17,Filter18,Filter19}; // Value-Defintions of the different String values
	std::map<std::string, Filters> FilterValue;	// Map to associate the strings with the enum values
	ftk::Image::Pointer getImage(void);
	std::vector<double> getParams(std::string id); 
	std::string filtername;

	int getChannelNumber(void);
	void InitializeFilters(void); 

private:

	QComboBox *channelCombo;
	QLabel *QTParamLabel1;
	QLabel *QTParamLabel2;
	QLabel *QTParamLabel3;
	QLabel *QTParamLabel4;
	QLabel *QTParamLabel5;
	
	QLineEdit *QTParam1;
	QLineEdit *QTParam2;	
	QLineEdit *QTParam3;	
	QLineEdit *QTParam4;		
	
	QCheckBox *QTParam5;	
	QCheckBox *QTParam6;	
	
	QString lastPath;
	QPushButton *okButton;
	
	double paramVal1;
	double paramVal2;	
	double paramVal3;
	double paramVal4;
	double paramVal5;
	std::vector<double> parameters;
};


#endif