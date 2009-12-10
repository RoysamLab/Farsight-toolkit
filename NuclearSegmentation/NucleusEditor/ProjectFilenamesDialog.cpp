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
#include "ProjectFilenamesDialog.h"

//Constructors:
ProjectFilenamesDialog::ProjectFilenamesDialog(QString path, ftk::ProjectFiles * files, QString projfile, QWidget * parent)
: QDialog(parent)
{
	pFiles = files;
	lastPath = path;

	QString buttonText = tr("Change...");
	int labelWidths = 120;

	QVBoxLayout *masterLayout = new QVBoxLayout;

	QHBoxLayout *pLayout = new QHBoxLayout;
	projectLabel = new QLabel(tr("Project File: "));
	projectLabel->setFixedWidth(labelWidths);
	projectFile = new QLabel(projfile);
	projectFile->setFrameShadow(QFrame::Sunken);
	projectFile->setFrameShape(QFrame::StyledPanel);
	pLayout->addWidget(projectLabel);
	pLayout->addWidget(projectFile);
	//pLayout->addStretch(5);
	masterLayout->addLayout(pLayout);

	QHBoxLayout *nLayout = new QHBoxLayout;
	nameLabel = new QLabel(tr("Project Name: "));
	nameLabel->setFixedWidth(labelWidths);
	nameText = new QLabel(QString::fromStdString(files->name));
	nameText->setFrameShadow(QFrame::Sunken);
	nameText->setFrameShape(QFrame::StyledPanel);
	nameButton = new QPushButton(buttonText);
	connect(nameButton, SIGNAL(clicked()), this, SLOT(changeName()));
	nLayout->addWidget(nameLabel);
	nLayout->addWidget(nameText);
	nLayout->addWidget(nameButton);
	masterLayout->addLayout(nLayout);

	masterLayout->addStretch(5);

	QHBoxLayout *iLayout = new QHBoxLayout;
	inputLabel = new QLabel(tr("Input: "));
	inputLabel->setFixedWidth(labelWidths);
	inputFile = new QLabel(QString::fromStdString(files->input));
	inputFile->setFrameShadow(QFrame::Sunken);
	inputFile->setFrameShape(QFrame::StyledPanel);
	inputButton = new QPushButton(buttonText);
	connect(inputButton, SIGNAL(clicked()), this, SLOT(changeInput()));
	iLayout->addWidget(inputLabel);
	iLayout->addWidget(inputFile);
	iLayout->addWidget(inputButton);
	masterLayout->addLayout(iLayout);

	QHBoxLayout *oLayout = new QHBoxLayout;
	outputLabel = new QLabel(tr("Output: "));
	outputLabel->setFixedWidth(labelWidths);
	outputFile = new QLabel(QString::fromStdString(files->output));
	outputFile->setFrameShadow(QFrame::Sunken);
	outputFile->setFrameShape(QFrame::StyledPanel);
	outputButton = new QPushButton(buttonText);
	connect(outputButton, SIGNAL(clicked()), this, SLOT(changeOutput()));
	oLayout->addWidget(outputLabel);
	oLayout->addWidget(outputFile);
	oLayout->addWidget(outputButton);
	masterLayout->addLayout(oLayout);

	QHBoxLayout *lLayout = new QHBoxLayout;
	logLabel = new QLabel(tr("Edit Log: "));
	logLabel->setFixedWidth(labelWidths);
	logFile = new QLabel(QString::fromStdString(files->log));
	logFile->setFrameShadow(QFrame::Sunken);
	logFile->setFrameShape(QFrame::StyledPanel);
	logButton = new QPushButton(buttonText);
	connect(logButton, SIGNAL(clicked()), this, SLOT(changeLog()));
	lLayout->addWidget(logLabel);
	lLayout->addWidget(logFile);
	lLayout->addWidget(logButton);
	masterLayout->addLayout(lLayout);

	QHBoxLayout *dLayout = new QHBoxLayout;
	definitionLabel = new QLabel(tr("Project Definition: "));
	definitionLabel->setFixedWidth(labelWidths);
	definitionFile = new QLabel(QString::fromStdString(files->definition));
	definitionFile->setFrameShadow(QFrame::Sunken);
	definitionFile->setFrameShape(QFrame::StyledPanel);
	definitionButton = new QPushButton(buttonText);
	connect(definitionButton, SIGNAL(clicked()), this, SLOT(changeDefinition()));
	dLayout->addWidget(definitionLabel);
	dLayout->addWidget(definitionFile);
	dLayout->addWidget(definitionButton);
	masterLayout->addLayout(dLayout);

	QHBoxLayout *tLayout = new QHBoxLayout;
	tableLabel = new QLabel(tr("Table: "));
	tableLabel->setFixedWidth(labelWidths);
	tableFile = new QLabel(QString::fromStdString(files->table));
	tableFile->setFrameShadow(QFrame::Sunken);
	tableFile->setFrameShape(QFrame::StyledPanel);
	tableButton = new QPushButton(buttonText);
	connect(tableButton, SIGNAL(clicked()), this, SLOT(changeTable()));
	tLayout->addWidget(tableLabel);
	tLayout->addWidget(tableFile);
	tLayout->addWidget(tableButton);
	masterLayout->addLayout(tLayout);

	masterLayout->addStretch(5);

	QHBoxLayout *okLayout = new QHBoxLayout;
	okButton = new QPushButton(tr("Done"));
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	okLayout->addStretch(20);
	okLayout->addWidget(okButton);
	masterLayout->addLayout(okLayout);

	this->setLayout(masterLayout);

	//this->setupDefaults();
}

void ProjectFilenamesDialog::setupDefaults()
{
	QString projfile = projectFile->text();
	QString path = QFileInfo(projfile).absolutePath();
	QString name = QFileInfo(projfile).baseName();

	
}

void ProjectFilenamesDialog::changeName(void)
{
	bool ok;
	QString text = QInputDialog::getText(this, tr("Name"),
                                          tr("Enter project name:"), QLineEdit::Normal, nameText->text(), &ok);
     if (ok && !text.isEmpty() && text!=nameText->text())
	 {
         nameText->setText(text);
	 }
}

void ProjectFilenamesDialog::changeInput(void)
{
	QString	filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("Input File(*.*)"));
	
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	if(filename != inputFile->text())
	{
		inputFile->setText(filename);
		pFiles->input = filename.toStdString();
		pFiles->inputSaved = false;
	}
}

void ProjectFilenamesDialog::changeOutput(void)
{
	QString	filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("Output File(*.*)"));
	
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	if(filename != outputFile->text())
	{
		outputFile->setText(filename);
		pFiles->output = filename.toStdString();
		pFiles->outputSaved = false;
	}
}
void ProjectFilenamesDialog::changeLog(void)
{
	QString	filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("Edit Log File(*.txt)"));
	
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	if(filename != logFile->text())
	{
		logFile->setText(filename);
		pFiles->log = filename.toStdString();
		pFiles->logSaved = false;
	}
}
void ProjectFilenamesDialog::changeDefinition(void)
{
	QString	filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("Definition File(*.xml)"));
	
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	if(filename != definitionFile->text())
	{
		definitionFile->setText(filename);
		pFiles->definition = filename.toStdString();
		pFiles->definitionSaved = false;
	}
}
void ProjectFilenamesDialog::changeTable(void)
{
	QString	filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("Table File(*.txt)"));
	
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	if(filename != tableFile->text())
	{
		tableFile->setText(filename);
		pFiles->table = filename.toStdString();
		pFiles->tableSaved = false;
	}
}