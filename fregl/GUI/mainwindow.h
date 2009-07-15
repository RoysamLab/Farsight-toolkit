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

#ifndef IMAGEWINDOW_H
#define IMAGEWINDOW_H

#include <QMainWindow>

#include "gnrcimagedialog.h"

class QAction;
class QActionGroup;
class QLabel;
class QMenu;
class QToolBar;
class QWorkspace;
class TransactionThread;


class MainWindow : public QMainWindow
{
  Q_OBJECT;
public:
  MainWindow();
  
protected:
  void closeEvent(QCloseEvent *event);
                                     
private slots:
  void    loadImage();
  void    saveImageAs();
  void    saveROIAs();
  void    deleteROI();
  void    defineROI();

private:
  void    createActions();
  void    createMenus();
  void    createStatusBar();
  
  QMenu   *fileMenu;
  QAction *openAction;
  QAction *saveAsAction;
  QAction *deleteAction;
  QAction *defineAction;
  QAction *exitAction;
  
  QAction *aboutAction;
  QLabel  *readyLabel;
  
  QString  fileFilters;
  
  QWorkspace   *workspace;
  QWidgetList   windows;
  QToolBar     *fileToolBar;
  QActionGroup *windowActionGroup;
  GNRCImageDialog *inDlg;

private: 
  GNRCImageDialog *createGNRCImage();

};

#endif
