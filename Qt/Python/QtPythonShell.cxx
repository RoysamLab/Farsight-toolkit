/*=========================================================================

   Program: ParaView
   Module:    $RCSfile: QtPythonShell.cxx,v $

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

#include <vtkPython.h> // python first

#include "QtPythonConfig.h"

#include "pqConsoleWidget.h"
#include "QtPythonShell.h"
#include "vtkCommand.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkPVPythonInteractiveInterpretor.h"

#include <QCoreApplication>
#include <QResizeEvent>
#include <QTextCharFormat>
#include <QVBoxLayout>

/////////////////////////////////////////////////////////////////////////
// QtPythonShell::QtImplementation

struct QtPythonShell::QtImplementation
{
  QtImplementation(QWidget* Parent) 
    : Console(Parent), Interpreter(0)
  {
    this->VTKConnect = vtkEventQtSlotConnect::New();
  }

  void initialize(int argc, char* argv[])
  {
    this->destroyInterpretor();
    this->Interpreter = vtkPVPythonInteractiveInterpretor::New();
    this->Interpreter->SetCaptureStreams(true);
    this->Interpreter->InitializeSubInterpretor(argc, argv);
    this->Interpreter->MakeCurrent();
    
    // Setup Python's interactive prompts
    PyObject* ps1 = PySys_GetObject(const_cast<char*>("ps1"));
    if(!ps1)
      {
      PySys_SetObject(const_cast<char*>("ps1"), ps1 = PyString_FromString(">>> "));
      Py_XDECREF(ps1);
      }

    PyObject* ps2 = PySys_GetObject(const_cast<char*>("ps2"));
    if(!ps2)
      {
      PySys_SetObject(const_cast<char*>("ps2"), ps2 = PyString_FromString("... "));
      Py_XDECREF(ps2);
      }
    this->Interpreter->ReleaseControl();
    this->MultilineStatement = false;
  }

  ~QtImplementation()
  {
    this->VTKConnect->Disconnect();
    this->VTKConnect->Delete();
    this->destroyInterpretor();
  }

  void destroyInterpretor()
    {
    if (this->Interpreter)
      {
      QTextCharFormat format = this->Console.getFormat();
      format.setForeground(QColor(255, 0, 0));
      this->Console.setFormat(format);
      this->Console.printString("\n... restarting ...\n");
      format.setForeground(QColor(0, 0, 0));
      this->Console.setFormat(format);

      this->Interpreter->MakeCurrent();

      // Restore Python's original stdout and stderr
      PySys_SetObject(const_cast<char*>("stdout"), PySys_GetObject(const_cast<char*>("__stdout__")));
      PySys_SetObject(const_cast<char*>("stderr"), PySys_GetObject(const_cast<char*>("__stderr__")));
      this->Interpreter->ReleaseControl();
      this->Interpreter->Delete();
      }
    this->Interpreter = 0;
    }

  void executeCommand(const QString& Command)
  {
    this->MultilineStatement = 
      this->Interpreter->Push(Command.toAscii().data());
  }

  void promptForInput()
  {
    QTextCharFormat format = this->Console.getFormat();
    format.setForeground(QColor(0, 0, 0));
    this->Console.setFormat(format);

    this->Interpreter->MakeCurrent();
    
    if(this->MultilineStatement || this->Console.getReadingFromStdin())
      {
      this->Console.prompt(PyString_AsString(PySys_GetObject(const_cast<char*>("ps2"))));
      }
    else
      {
      this->Console.prompt(PyString_AsString(PySys_GetObject(const_cast<char*>("ps1"))));
      }
    this->Interpreter->ReleaseControl();
  }

  /// Provides a console for gathering user input and displaying 
  /// Python output
  pqConsoleWidget Console;

  /// Indicates if the last statement processes was incomplete.
  bool MultilineStatement;

  /// Separate Python interpreter that will be used for this shell
  vtkPVPythonInteractiveInterpretor* Interpreter;

  vtkEventQtSlotConnect *VTKConnect;
};

/////////////////////////////////////////////////////////////////////////
// QtPythonShell

QtPythonShell::QtPythonShell(QWidget* Parent) :
  QWidget(Parent),
  Implementation(new QtImplementation(this))
{
  QVBoxLayout* const boxLayout = new QVBoxLayout(this);
  boxLayout->setMargin(0);
  boxLayout->addWidget(&this->Implementation->Console);

  this->setObjectName("pythonShell");
  
  QObject::connect(
    &this->Implementation->Console, SIGNAL(executeCommand(const QString&)), 
    this, SLOT(onExecuteCommand(const QString&)));
}

QtPythonShell::~QtPythonShell()
{
  delete this->Implementation;
}

void QtPythonShell::setExecutable(char *c)
{
  this->argv0 = c;
}

void QtPythonShell::initializeInterpretor()
{
  this->initializeInterpretor(1, (char **)&this->argv0);
}

void QtPythonShell::initializeInterpretor(int argc, char* argv[])
{
  this->Implementation->VTKConnect->Disconnect();

  this->Implementation->initialize(argc, argv);
  QTextCharFormat format = this->Implementation->Console.getFormat();
  format.setForeground(QColor(0, 0, 255));
  this->Implementation->Console.setFormat(format);
  this->Implementation->Console.printString(
    QString("Python %1 on %2\n").arg(Py_GetVersion()).arg(Py_GetPlatform()));
  this->promptForInput();

  this->Implementation->VTKConnect->Connect(
    this->Implementation->Interpreter, vtkCommand::ErrorEvent, 
    this, SLOT(printStderr(vtkObject*, unsigned long, void*, void*))); 
  this->Implementation->VTKConnect->Connect(
    this->Implementation->Interpreter, vtkCommand::WarningEvent, 
    this, SLOT(printStdout(vtkObject*, unsigned long, void*, void*))); 

  this->Implementation->VTKConnect->Connect(
    this->Implementation->Interpreter, vtkCommand::StartPickEvent, 
    this, SLOT(readStdin(vtkObject*, unsigned long, void*, void*))); 
}

void QtPythonShell::clear()
{
  this->Implementation->Console.clear();
  this->Implementation->promptForInput();
}

void QtPythonShell::executeScript(const QString& script)
{
  this->printStdout("\n");
  emit this->executing(true);  
  this->Implementation->Interpreter->RunSimpleString(
    script.toAscii().data());
  emit this->executing(false);
  this->Implementation->promptForInput();
}

void QtPythonShell::printStdout(vtkObject* o, unsigned long l,
                                void* v, void* calldata)
{
  const char* text = reinterpret_cast<const char*>(calldata);
  this->printStdout(text);
  this->Implementation->Interpreter->ClearMessages();
}

void QtPythonShell::printStderr(vtkObject* o, unsigned long l,
                                void* v, void* calldata)
{
  const char* text = reinterpret_cast<const char*>(calldata);
  this->printStderr(text);
  this->Implementation->Interpreter->ClearMessages();
}

void QtPythonShell::readStdin(vtkObject* o, unsigned long l,
                                void* v, void* calldata)
{
  this->readStdin();
}

void QtPythonShell::readStdin()
{
  this->Implementation->Console.startReadingFromStdin();
//  this->Implementation->Interpreter->RawInput();
//  emit this->executing(false);
//  this->MultilineStatement = true;
//  this->promptForInput();
}

void QtPythonShell::printStdout(const QString& text)
{
  QTextCharFormat format = this->Implementation->Console.getFormat();
  format.setForeground(QColor(0, 150, 0));
  this->Implementation->Console.setFormat(format);
  
  this->Implementation->Console.printString(text);
  
  QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

void QtPythonShell::printMessage(const QString& text)
{
  QTextCharFormat format = this->Implementation->Console.getFormat();
  format.setForeground(QColor(0, 0, 150));
  this->Implementation->Console.setFormat(format);
  
  this->Implementation->Console.printString(text);
}

void QtPythonShell::printStderr(const QString& text)
{
  QTextCharFormat format = this->Implementation->Console.getFormat();
  format.setForeground(QColor(255, 0, 0));
  this->Implementation->Console.setFormat(format);
  
  this->Implementation->Console.printString(text);
  
  QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

void QtPythonShell::onExecuteCommand(const QString& Command)
{
  QString command = Command;
  command.replace(QRegExp("\\s*$"), "");
  this->internalExecuteCommand(command);
  this->promptForInput();
}

void QtPythonShell::promptForInput()
{
  this->Implementation->promptForInput();
}

void QtPythonShell::internalExecuteCommand(const QString& command)
{
  emit this->executing(true);  
  this->Implementation->executeCommand(command);
  emit this->executing(false);
}
