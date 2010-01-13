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
#ifndef FTKPREPROCESSDIALOG_H
#define FTKPREPROCESSDIALOG_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>

//Farsight Includes:
#include "ftkImage/ftkImage.h"

// Preprocess
//#include "ftkPreprocess.h"

#include <map>
#include <vector>

class ftkPreprocessDialog : public QDialog
{
	Q_OBJECT
public:
	ftkPreprocessDialog(QVector<QString> chs, std::string id, ftk::Image::Pointer img, QWidget *parent = 0);
	
public slots:
	void accept();

private:
	enum Filters {Filter1,Filter2,Filter3,Filter4,Filter5,Filter6,Filter7,Filter8,Filter9,Filter10,Filter11,Filter12,Filter13,Filter14,Filter15,Filter16,Filter17,Filter18,Filter19,Filter20}; // Value-Defintions of the different String values

	//functions:
	std::vector<double> getParams(std::string id); 
	int getChannelNumber(void);
	void InitializeFilters(void);
	void doPreprocess(void);

	//variables:
	std::string filtername;
	ftk::Image::Pointer myImg;
	std::map<std::string, Filters> FilterValue;	// Map to associate the strings with the enum valuesk

	QLabel *channelLabel;
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
	
	QPushButton *okButton;
	QPushButton *cancelButton;
	
	double paramVal1;
	double paramVal2;	
	double paramVal3;
	double paramVal4;
	double paramVal5;
	std::vector<double> parameters;
};

#endif