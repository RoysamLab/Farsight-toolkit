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
#include "PreprocessDialog.h"

PreprocessDialog::~PreprocessDialog()
{
	if(prep)
		delete prep;
}

PreprocessDialog::PreprocessDialog(QString lastPath, QWidget *parent)
: QDialog(parent)
{
	wrapper = QString("<Preprocess>\n</Preprocess>");
	lpath = lastPath;
	filename = "";
	prep = new ftk::Preprocess();

	QVBoxLayout * masterLayout = new QVBoxLayout();

	filterCombo = new QComboBox();
	filterCombo->clear();
	filterCombo->addItem("Insert Filter...");
	std::map<std::string, std::string>::iterator it;
	for(it=ftk::Preprocess::filterMap.begin() ; it!=ftk::Preprocess::filterMap.end(); it++)
	{
		filterCombo->addItem(QString::fromStdString(((*it).first)));
	}
	connect(filterCombo, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(insertFilter(const QString &)));
	masterLayout->addWidget(filterCombo);

	textEdit = new QPlainTextEdit(wrapper);
	QTextCursor cursor = textEdit->textCursor();
	cursor.movePosition(QTextCursor::Down);
	textEdit->setTextCursor(cursor);

	masterLayout->addWidget(textEdit);

	QHBoxLayout * buttonLayout = new QHBoxLayout();
	buttonLayout->addStretch(10);

	loadButton = new QPushButton(tr("Load"));
	connect(loadButton, SIGNAL(clicked()), this, SLOT(loadPipe()));
	buttonLayout->addWidget(loadButton);

	saveButton = new QPushButton(tr("Save"));
	connect(saveButton, SIGNAL(clicked()), this, SLOT(savePipe()));
	buttonLayout->addWidget(saveButton);

	buttonLayout->addStretch(1);

	processButton = new QPushButton(tr("Exit"));
	connect(processButton, SIGNAL(clicked()), this, SLOT(finalizePreprocessing()));
	buttonLayout->addWidget(processButton);

	masterLayout->addLayout(buttonLayout);

	this->setLayout(masterLayout);

	this->setWindowTitle(tr("Preprocessing"));
}

void PreprocessDialog::SetImage(ftk::Preprocess::ImageType3D::Pointer im)
{
	prep->SetImage(im);
}

void PreprocessDialog::insertFilter(const QString & text)
{
	if( text == "Insert Filter..." )
		return;

	std::string filter = text.toStdString();

	textEdit->insertPlainText(QString::fromStdString(ftk::Preprocess::filterMap[filter]));

	filterCombo->setCurrentIndex(0);
}

void PreprocessDialog::loadPipe()
{
	QString fname = QFileDialog::getOpenFileName(this, tr("Select File..."), lpath, tr("All Files (*.*)"));
	if(fname == "")  return;

	filename = fname;
	lpath = QFileInfo(filename).absolutePath();


	QFile file(filename);
	file.open(QFile::ReadOnly | QFile::Text);
	QTextStream in(&file);
	textEdit->setPlainText(in.readAll());

	QTextCursor cursor = textEdit->textCursor();
	cursor.movePosition(QTextCursor::End);
	cursor.movePosition(QTextCursor::StartOfLine);
	textEdit->setTextCursor(cursor);
}

void PreprocessDialog::savePipe()
{
	QString fname = QFileDialog::getSaveFileName(this, tr("Save As..."),lpath, tr("All Files (*.*)"));
	if(fname == "")  return;

	filename = fname;
	lpath = QFileInfo(filename).absolutePath();

	QFile file(filename);
	file.open(QFile::WriteOnly | QFile::Text);
	QTextStream out(&file);
	out << textEdit->toPlainText();
}

void PreprocessDialog::finalizePreprocessing()
{
	if(!prep)
		this->reject();	

	this->accept();
}

void PreprocessDialog::Process()
{
	if( filename.isEmpty() )
	{
		QFile file("temp.pipe");
		file.open(QFile::WriteOnly | QFile::Text);
		QTextStream out(&file);
		out << textEdit->toPlainText();
		file.close();

		prep->RunPipe("temp.pipe");

		file.remove();
	}
	else
	{
		prep->RunPipe(filename.toStdString());
	}	
}

ftk::Preprocess::ImageType3D::Pointer PreprocessDialog::GetImage()
{
	if(!prep)
		return NULL;
	else
		return prep->GetImage();
}



