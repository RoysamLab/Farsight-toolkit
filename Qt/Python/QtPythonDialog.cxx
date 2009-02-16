/*=========================================================================

   Program: ParaView
   Module:    $RCSfile: QtPythonDialog.cxx,v $

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2. 

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include "QtPythonDialog.h"
#include <QFile>
#include <QFileDialog>
#include <QtDebug>

//////////////////////////////////////////////////////////////////////

QtPythonDialog::QtPythonDialog(QWidget* Parent, const char *c) :
  QDialog(Parent),
  Implementation(new QtImplementation())
{
  this->argv0 = new char[strlen(c)+1];
  strcpy(this->argv0, c);
  this->Implementation->Ui.setupUi(this);
  this->Implementation->Ui.shellWidget->setExecutable(this->argv0);
  this->setObjectName("pythonDialog");
  this->setWindowTitle(tr("Python Shell"));
  
  QObject::connect(
    this->Implementation->Ui.clear,
    SIGNAL(clicked()),
    this,
    SLOT(clearConsole()));
    
  QObject::connect(
    this->Implementation->Ui.runScript,
    SIGNAL(clicked()),
    this,
    SLOT(runScript()));
    
  QObject::connect(
    this->Implementation->Ui.reset,
    SIGNAL(clicked()),
    this,
    SLOT(initializeInterpretor()));                

  QObject::connect(
    this->Implementation->Ui.shellWidget,
    SIGNAL(executing(bool)),
    this->Implementation->Ui.runScript,
    SLOT(setDisabled(bool)));

  QObject::connect(
    this->Implementation->Ui.shellWidget,
    SIGNAL(executing(bool)),
    this->Implementation->Ui.clear,
    SLOT(setDisabled(bool)));

  QObject::connect(
    this->Implementation->Ui.shellWidget,
    SIGNAL(executing(bool)),
    this->Implementation->Ui.close,
    SLOT(setDisabled(bool)));
}

QtPythonDialog::~QtPythonDialog()
{
  delete Implementation;
  delete this->argv0;
}

void QtPythonDialog::runScript()
{
  QFileDialog* const dialog = new QFileDialog(
    this,
    tr("Run Script"),
    QString(),
    QString(tr("Python Script (*.py);;All files (*)")));
    
  dialog->setObjectName("PythonShellRunScriptDialog");
  dialog->setFileMode(QFileDialog::ExistingFiles);
  QObject::connect(
    dialog,
    SIGNAL(filesSelected(const QStringList&)), 
    this,
    SLOT(runScript(const QStringList&)));
  dialog->show(); 
}

void QtPythonDialog::runScript(const QStringList& files)
{
  for(int i = 0; i != files.size(); ++i)
    {
    QFile file(files[i]);
    if(file.open(QIODevice::ReadOnly))
      {
      this->Implementation->Ui.shellWidget->executeScript(
        file.readAll().data());
      }
    else
      {
      qCritical() << "Error opening " << files[i];
      }
    }
}

void QtPythonDialog::print(const QString& msg)
{
  this->Implementation->Ui.shellWidget->printMessage(msg);
}

void QtPythonDialog::runString(const QString& str)
{
  this->Implementation->Ui.shellWidget->executeScript(str);
}

void QtPythonDialog::clearConsole()
{
  this->Implementation->Ui.shellWidget->clear();
}

void QtPythonDialog::initializeInterpretor()
{
  this->Implementation->Ui.shellWidget->initializeInterpretor();
  emit this->interpreterInitialized();
}

