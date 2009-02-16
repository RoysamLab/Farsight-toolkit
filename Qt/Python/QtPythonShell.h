/*=========================================================================

   Program: ParaView
   Module:    $RCSfile: QtPythonShell.h,v $

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

#ifndef _QtPythonShell_h
#define _QtPythonShell_h

#include "QtPythonExport.h"
#include <QWidget>
/**
  Qt widget that provides an interactive "shell" interface to an embedded Python interpreter.
  You can put an instance of QtPythonShell in a dialog or a window, and the user will be able
  to enter Python commands and see their output, while the UI is still responsive.
  
  \sa pqConsoleWidget, pqPythonDialog
*/  
  
class vtkObject;
class QTPYTHON_EXPORT QtPythonShell :
  public QWidget
{
  Q_OBJECT
  friend class pqPythonShell;
  
public:
  QtPythonShell(QWidget* Parent);
  ~QtPythonShell();


  /// Initializes the interpretor. If an interpretor is already setup (by an
  /// earlier call to this method), it will be destroyed.
  virtual void initializeInterpretor(int argc, char* argv[]);
  virtual void initializeInterpretor();
  void setExecutable(char *c);

  /// Prints some text on the shell.
  void printMessage(const QString&);

signals:
  void executing(bool);

public slots:
  virtual void clear();
  virtual void executeScript(const QString& text);

private slots:
  virtual void printStderr(vtkObject* o, unsigned long l,
                           void* v, void* calldata);
  virtual void printStdout(vtkObject* o, unsigned long l,
                           void* v, void* calldata);
  virtual void readStdin(vtkObject* o, unsigned long l,
                         void* v, void* calldata);

  virtual void onExecuteCommand(const QString& Command);

private:
  QtPythonShell(const QtPythonShell&);
  QtPythonShell& operator=(const QtPythonShell&);

  void promptForInput();
  void internalExecuteCommand(const QString&);

  struct QtImplementation;
  QtImplementation* const Implementation;

  void printStderr(const QString&);
  void printStdout(const QString&);
  void readStdin();
  char *argv0;
};

#endif // !_QtPythonShell_h

