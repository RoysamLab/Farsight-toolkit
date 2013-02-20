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
#ifndef PROJECTFILENAMESDIALOG_H
#define PROJECTFILENAMESDIALOG_H

#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QDialog>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QFileDialog>
#include <QtGui/QInputDialog>
#include <QtGui/QMessageBox>
#include <QtCore/QFileInfo>

#include "ftkProjectFiles.h"

class ProjectFilenamesDialog : public QDialog
{
    Q_OBJECT;

public:
	ProjectFilenamesDialog(ftk::ProjectFiles * files, QWidget * parent = 0);

signals:

protected:

public slots:
	void accept();

private slots:
	void changeName(void);
	void changePath(void);
	void changeInput(void);
	void changeOutput(void);
	void changeLog(void);
	void changeDefinition(void);
	void changeTable(void);
	void changeAdjTables(void);
	void setupDefaults(void);
    
private:
	ftk::ProjectFiles * pFiles;
	QString lastPath;

	QLabel *projectLabel;
	QLabel *projectFile;

	QLabel *nameLabel;
	QLabel *nameText;
	QPushButton *nameButton;

	QLabel *pathLabel;
	QLabel *inputPath;
	QPushButton *pathButton;

	QLabel *inputLabel;
	QLabel *inputFile;
	QPushButton *inputButton;

	QLabel *outputLabel;
	QLabel *outputFile;
	QPushButton *outputButton;

	QLabel *logLabel;
	QLabel *logFile;
	QPushButton *logButton;

	QLabel *definitionLabel;
	QLabel *definitionFile;
	QPushButton *definitionButton;

	QLabel *tableLabel;
	QLabel *tableFile;
	QPushButton *tableButton;

	QLabel *adjTablesLabel;
	QLabel *adjTablesFile;
	QPushButton *adjTablesButton;

	QPushButton *quitButton;
	QPushButton *okButton;
};

#endif
