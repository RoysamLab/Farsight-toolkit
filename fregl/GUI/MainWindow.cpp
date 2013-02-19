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

#include <QtGui>
#include <QChar>
#include "mainwindow.h"

//#define DEBUG_P(...)

////////////////////////////////////////////////////////////////////////////////
GNRCImageDialog *MainWindow::createGNRCImage( )
{
  GNRCImageDialog *idlg = new GNRCImageDialog;
  workspace->addWindow(idlg);
  workspace->addAction(idlg->windowMenuAction());
  return idlg;
}

MainWindow::MainWindow()
{
  workspace = new QWorkspace;
  setCentralWidget(workspace);
  
  /*
  connect(workspace, SIGNAL(windowActivated(QWidget *)),
          this, SLOT(updateMenus()));
  */

  createActions();
  createMenus();
  createStatusBar();
  statusBar()->showMessage(tr("Ready"),2000);
  
  fileFilters = tr("All Files (*.*)\n"
                   "TIF Files (*.tif *.tiff)\n"
                   "PNG Files (*.png)\n"
                   "PIC Files (*.pic)\n");
  
  setWindowTitle(tr("FARSIGHT MONTAGE NAVIGATOR"));

  inDlg = 0;
}


/*----------------------------------------------------------------------------*/
/* function: createActions()                                                  */
/*                                                                            */
/* This function is used to create Actions that are associated with various   */
/* Menu items. By Actions, we imply actual routines that will be called to    */
/* obtain some desired result. In order to create an Action, we need to do the*/
/* following:                                                                 */
/* 1.) Define a QAction (e.g., QAction *openAction)                           */
/* 2.) Label the QAction element (e.g., openAction = new QAction(QIcon(":src/ */
/*     images/open.png"), tr("&Open..."), this). The QIcon argumenet is       */
/*     optional.                                                              */
/* 3.) Add optional "setShortcut" and "setStatusTip".                         */
/* 4.) Finally, bind this item with a "connect" that essentially calls the    */
/*     module to implement the operation (e.g.,                               */
/*     connect(openAction, SIGNAL(triggered()), this, SLOT(loadImage())). In  */
/*     this example, "loadImage()" is the module that is being called. These  */
/*     modules should be defined as "private" operators in the main class.    */
/*     The actual routines performing the operations (e.g., an image          */
/*     thresholding operation) must be accessed from within the called module.*/
/* Finally, after all these action's, we bind them to a "QActionGroup".       */ 
/*----------------------------------------------------------------------------*/
void MainWindow::createActions()
{
  openAction = new QAction(tr("&Load Montage XML..."), this);
  openAction->setShortcut(tr("Ctrl+O"));
  openAction->setStatusTip(tr("Open an existing image file"));
  connect(openAction, SIGNAL(triggered()), this, SLOT(loadImage()));
  
  //saveAsAction = new QAction(tr("Save &As..."), this);
  //saveAsAction->setStatusTip(tr("Save the image under a new name"));
  //connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveImageAs()));
  
  saveAsAction = new QAction(tr("Save &ROI..."), this);
  saveAsAction->setStatusTip(tr("Save the selected region of interest as a new 3D image"));
  connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveROIAs()));
  
  exitAction = new QAction(tr("E&xit"), this);
  exitAction->setShortcut(tr("Ctrl+Q"));
  exitAction->setStatusTip(tr("Exit the application"));
  connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
  
  deleteAction = new QAction(tr("Delete ROI..."), this);
  deleteAction->setShortcut(tr("Ctrl+D"));
  deleteAction->setStatusTip(tr("Delete the selected ROI"));
  connect(deleteAction, SIGNAL(triggered()), this, SLOT(deleteROI()));

  defineAction = new QAction(tr("Define &ROI..."), this);
  defineAction->setStatusTip(tr("Define new ROI"));
  connect(defineAction, SIGNAL(triggered()), this, SLOT(defineROI()));
  
  windowActionGroup = new QActionGroup(this);
}

/*----------------------------------------------------------------------------*/
/* function: createMenus()                                                    */
/*                                                                            */
/* This function is used to create Menus that are associated with various     */
/* functions in the application. In order to add a menu, we need to do the    */
/* following:                                                                 */
/* 1.) Define a QMenu type (e.g., QMenu *fileMenu) and add it to menuBar()    */
/* 2.) Define QAction elements (e.g., QAction *openAction) associated with    */
/*     each QMenu                                                             */
/* 3.) Add a separator (menuBar()->addSeparator() after each menu item        */
/*----------------------------------------------------------------------------*/
void MainWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAction);
  fileMenu->addAction(saveAsAction);
  fileMenu->addAction(deleteAction);
  //fileMenu->addAction(defineAction);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAction);
}

////////////////////////////////////////////////////////////////////////////////
void MainWindow::loadImage()
{
  QString xfileFilters = tr("XML Files (*.xml)\n");
  QString filename = QFileDialog::getOpenFileName( this,"Open the montage XML",".",xfileFilters);
  
  if(filename != "")
    {  
      //Create new Image Viewer. If there is already one in the main
      //window, delete. Only one window is allowed, because I don't
      //know how to control multiple windows.
      if (inDlg) delete inDlg; 
      inDlg = createGNRCImage();
      inDlg->loadXML(filename);
      inDlg->show();
    }
}

////////////////////////////////////////////////////////////////////////////////

void MainWindow::saveImageAs()
{
 
}

////////////////////////////////////////////////////////////////////////////////

void MainWindow::saveROIAs()
{
  //QString in_filename = QFileDialog::getOpenFileName( this,"Choose the 3D montage directory for cropping",".",fileFilters);
  QString in_filename = QFileDialog::getExistingDirectory(this,"Choose the 3D montage directory from which cropping is performed");
  
  QString out_fileFilters = tr("TIF Files (*.tif *.tiff)\n");
  QString out_filename = QFileDialog::getSaveFileName( this,"Save the ROI as",".",out_fileFilters);
  if (out_filename != "") {  
    inDlg->saveROIAs( in_filename, out_filename );
  }
}

////////////////////////////////////////////////////////////////////////////////

void MainWindow::deleteROI()
{
  inDlg->deleteROI();
}

//////////////////////////////////////////////////////////////////////////////

void MainWindow::defineROI()
{
  // A function to be implemented. A pope-up window should allow input
  // of the position of the top-left corner and the size of the box. 
}
////////////////////////////////////////////////////////////////////////////////
void MainWindow::closeEvent(QCloseEvent *event)
{
  event->accept();
}

void MainWindow::createStatusBar()
{
  readyLabel = new QLabel(tr(" Ready"));
  statusBar()->addWidget(readyLabel, 1);
}

