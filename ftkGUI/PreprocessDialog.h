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
#ifndef FTKPREPROCESSDIALOG2_H
#define FTKPREPROCESSDIALOG2_H

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
#include <QtGui/QPlainTextEdit>
#include <QtGui/QFileDialog>
#include <QtCore/QFile>
#include <QtCore/QTextStream>

#include <map>
#include <vector>

#include "ftkPreprocess2.h"

class PreprocessDialog : public QDialog
{
	Q_OBJECT
public:
	PreprocessDialog(QString lastPath = "", QWidget *parent = 0);
	~PreprocessDialog();
	void SetImage(ftk::Preprocess::ImageType3D::Pointer im);
	ftk::Preprocess::ImageType3D::Pointer GetImage();
	
public slots:
	void insertFilter(const QString & text);
	void loadPipe();
	void savePipe();
	void finalizePreprocessing();
	void Process();

private:

	QPlainTextEdit * textEdit;
	QComboBox * filterCombo;
	QPushButton * loadButton;
	QPushButton * saveButton;
	QPushButton * processButton;

	QString wrapper;
	QString lpath;
	QString filename;

	ftk::Preprocess * prep;
	
};

#endif