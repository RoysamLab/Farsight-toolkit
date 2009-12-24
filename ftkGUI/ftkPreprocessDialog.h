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
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"

//VTK includes:
#include "vtkQtTableView.h"

// Preprocess
//#include "ftkPreprocess.h"


class ftkPreprocessDialog : public QDialog
{
	Q_OBJECT
public:
	ftkPreprocessDialog(QVector<QString> chs, unsigned char id,QWidget *parent = 0);
	int getChannelNumber();
	ftk::Image::Pointer getImage();
	std::vector<double> getParams(unsigned char id); 
	QLabel *channelLabel;
	ftk::Image::Pointer myImg;

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