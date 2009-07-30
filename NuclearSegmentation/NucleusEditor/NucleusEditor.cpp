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

#include "NucleusEditor.h"

//*******************************************************************************
// NucleusEditor
//
// This widget is a main window for editing nuclei.  
//********************************************************************************

NucleusEditor::NucleusEditor(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{
	createMenus();
	createStatusBar();

	setWindowTitle(tr("Nucleus Editor"));

	segResult = new ftk::NuclearSegmentation();
	segWin = new SegmentationWindow();
	currentModel = NULL;

	lastPath = ".";

	tblWin.clear();
	pltWin.clear();
	hisWin=NULL;
	

	//DEMO
	this->pythonProcess = new QProcess();
	this->settings = new QSettings("RPI", "Farsight");

	this->resize(500,500);
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );
}

//******************************************************************************
// Here we just show a message in the status bar when loading
//******************************************************************************
void NucleusEditor::createStatusBar()
{
    statusLabel = new QLabel(tr(" Ready"));
    statusBar()->addWidget(statusLabel, 1);
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
/* 3.) Add a separator (menuBar()->addSeparator() after each menu group       */
/*																			  */
/*In order to create an Action, we need to do the							  */
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
/*	   Finally, after all these action's, we bind them to a "QActionGroup".       */ 
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NucleusEditor::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));

	loadAction = new QAction(tr("Load Image..."), this);
	loadAction->setStatusTip(tr("Load an image into the 5D image browser"));
	connect(loadAction, SIGNAL(triggered()), this, SLOT(loadImage()));
	fileMenu->addAction(loadAction);

	fileMenu->addSeparator();

	segmentAction = new QAction(tr("Start Segmentation Wizard..."), this);
	segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation Wizard"));
	connect(segmentAction,SIGNAL(triggered()),this,SLOT(segmentImage()));
	fileMenu->addAction(segmentAction);

	fileMenu->addSeparator();

	xmlAction = new QAction(tr("Load Result..."), this);
	xmlAction->setStatusTip(tr("Open an XML result file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu->addAction(xmlAction);

	/*
	xmlAction = new QAction(tr("Load From a dat file..."), this);
	xmlAction->setStatusTip(tr("Open a .dat file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadDatFile()));
	fileMenu->addAction(datAction);
	*/

	saveAction = new QAction(tr("Save Result As..."), this);
	saveAction->setStatusTip(tr("Save Changes (Edits, etc)"));
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(saveResult()));
	fileMenu->addAction(saveAction);

	fileMenu->addSeparator();

    exitAction = new QAction(tr("Exit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
    fileMenu->addAction(exitAction);

	//VIEW MENU
	viewMenu = menuBar()->addMenu(tr("&View"));
	newScatterAction = new QAction(tr("New Scatter"), this);
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

	pythonAction = new QAction(tr("Open Python Window"), this);
	pythonAction->setStatusTip(tr("Start your favorite python interpreter"));
	pythonAction->setShortcut(tr("Ctrl+P"));
	connect(pythonAction, SIGNAL(triggered()), this, SLOT(OpenPythonWindow()));
	viewMenu->addAction(pythonAction);

	//EDITING MENU	
	editMenu = menuBar()->addMenu(tr("&Editing"));
	// There is nothing to edit initially. Disabled, merge,delete, and split
	// Enable them after loading results!
	//editMenu->setEnabled(false);

	mergeAction = new QAction(tr("Merge Cells"), this);
	mergeAction->setStatusTip(tr("Merge Cells"));
	mergeAction->setEnabled(false);
	connect(mergeAction, SIGNAL(triggered()), this, SLOT(mergeCells()));	
	editMenu->addAction(mergeAction);

	editMenu->addSeparator();

	deleteAction = new QAction(tr("Delete Cells"), this);
	deleteAction->setStatusTip(tr("Deletes the selected cells"));
	deleteAction->setEnabled(false);
	connect(deleteAction,SIGNAL(triggered()),this,SLOT(deleteCells()));	
	editMenu->addAction(deleteAction);

	editMenu->addSeparator();

	splitAction = new QAction(tr("Split Cells"), this);
	splitAction->setStatusTip(tr("Splits the selected cells"));
	splitAction->setEnabled(false);
	connect(splitAction,SIGNAL(triggered()),this,SLOT(splitCells()));
	editMenu->addAction(splitAction);	

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);
}

//****************************************************************************
// SLOT: about()
//   A brief message about Farsight is displayed
//****************************************************************************
void NucleusEditor::about()
{
	    QMessageBox::about(this, tr("About FARSIGHT"),
            tr("<h2>FARSIGHT</h2>"
			   "<h3>Renssalear Polytechnic Institute</h3>"
			   "<a><u>http://www.farsight-toolkit.org</a></u>"
               ));
}

ftk::Image::Pointer NucleusEditor::NewFTKImage(std::string filename)
{
	ftk::Image::Pointer img = ftk::Image::New();
	img->LoadFile(filename);
	return img;
}

//******************************************************************************
//Reimplement closeEvent to also close all other windows in the application
//******************************************************************************
void NucleusEditor::closeEvent(QCloseEvent *event)
{

	//Save changes:
	if(segResult)
	{
		if(segResult->editsNotSaved)
		{
			QString msg = tr("Recent Edits not saved do you want to save them before exiting?");
			QMessageBox::StandardButton button = QMessageBox::information ( 0, tr("Exit"), \
				msg, QMessageBox::Yes | QMessageBox::No , QMessageBox::NoButton );

			if(button == QMessageBox::Yes)
			{
				if(	!this->saveResult()	)
				{
					event->ignore();
					return;
				}
			}
		}
	}

	//First clear the model and its associated windows
	clearModel();

	
	//Then Close all other windows
	foreach (QWidget *widget, qApp->topLevelWidgets()) 
	{
		if (this != widget)
		{
			if(widget->isVisible())
				widget->close();
		}
    }



	//Then close myself
	event->accept();
	
} 

bool NucleusEditor::saveResult()
{
	if(segResult)
	{
		QString name = QFileDialog::getSaveFileName(this, tr("Save File As"), lastPath, tr("XML (*.xml)") );
		if(name == "")
		{
			return false;
		}
		else
		{
			segResult->WriteToXML( name.toStdString() );
		}

		if(segResult->editsNotSaved)
		{
			segResult->SaveLabel();
		}
		
	}
	return true;
}

//***************************************************************************
//  THIS FUNCTION CLEARS THE MODEL/SEGMENTATION BY CLOSING OPEN WINDOW 
//  GROUPS ASSOCIATED WITH THEM, THEN DELETING THE MODEL AND SEGMENTATION
//***************************************************************************
void NucleusEditor::clearModel(void)
{	
	for(unsigned int p=0; p<pltWin.size(); ++p)
		if ((pltWin.at(p))->isVisible())
			(pltWin.at(p))->close();
	
	
	for(unsigned int p=0; p<tblWin.size(); ++p)	
		if ((tblWin.at(p))->isVisible())
			(tblWin.at(p))->close();

	//Close the histogram
	if (hisWin !=NULL) 
	{
		hisWin->close();
		hisWin=NULL;
	}

	pltWin.clear();
	tblWin.clear();

	
	if(currentModel)
	{
		delete currentModel;
		currentModel = NULL;
	}

}

//*********************************************************************************
// This function initializes the model and selection model
//*********************************************************************************
void NucleusEditor::newModel(void)
{	
	if(currentModel)
		clearModel();

	if(segResult) {
		currentModel = new SegmentationModel(segResult);
	}
}

//******************************************************************************
// This function loads a segmentation result from XML
// The XML file should tell where to find the original image/data
//   It also contains all of the feature and segmentation information
//   As well as what type of segmentation has been performed
//******************************************************************************
void NucleusEditor::loadResult(void)
{
	QString filename  = QFileDialog::getOpenFileName(this,"Choose a Result",lastPath, 
			tr("XML Files (*.xml)\n"));

    if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();

	lastPath = path;

	//segResult = new ftk::NuclearSegmentation();
	if ( !segResult->RestoreFromXML(filename.toStdString()) )
	{
		std::cerr << segResult->GetErrorMessage() << std::endl;
		return;
	}

	//Now I have objects stored in memory - put the features into the model
	newModel();
	CreateNewTableWindow();
	CreateNewPlotWindow();
	
	hisWin = new HistoWindow(currentModel->GetSelectionModel());
	hisWin->show();

	segWin->SetModels(currentModel);
	segWin->SetChannelImage(segResult->getDataImage());
	segWin->SetLabelImage(segResult->getLabelImage());
	if( this->centralWidget() != this->segWin )
		this->setCentralWidget(segWin);

	this->update();

	// Enable the menu items for editing
	mergeAction->setEnabled(true);
	deleteAction->setEnabled(true);
	splitAction->setEnabled(true);

	//Set the status
	//editStatus=true;

}

void NucleusEditor::mergeCells(void)
{
	currentModel->mergeTrigger();
}

void NucleusEditor::deleteCells(void)
{
	currentModel->deleteTrigger();

}

void NucleusEditor::splitCells(void)
{
	currentModel->splitTrigger();

}

// Added by Aytekin Vargun 6/03/09
//
/*
void NucleusEditor::loadDatFile(void)
{
	QString filename  = QFileDialog::getOpenFileName(this,"Choose a dat file",lastPath, 
			tr("dat Files (*.dat)\n"));

    if(fileName == "")
		return;

	if(currentModel)
		clearModel();

	lastPath = QFileInfo(fileName).absolutePath();

	//ImageBrowser5D *browse = new ImageBrowser5D(fileName);
	yousef_nucleus_seg *browse = new yousef_nucleus_seg();
	browse->readFromIDLFormat(fileName);
	this->setCentralWidget(browse);

}
*/

void NucleusEditor::segmentImage()
{
	// Disable the menu items for editing
	mergeAction->setEnabled(false);
	deleteAction->setEnabled(false);
	splitAction->setEnabled(false);

	if(currentModel)
		clearModel();

	NuclearSegmentationWizard *wizard = new NuclearSegmentationWizard();
	wizard->setParent(NULL);
	this->setCentralWidget(wizard);
	//wizard->show();
}

void NucleusEditor::loadImage()
{
	segWin = new SegmentationWindow();
	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							    "All Files (*.*)"));

    if(fileName == "")
		return;

	if(currentModel)
		clearModel();

	lastPath = QFileInfo(fileName).absolutePath();

	ImageBrowser5D *browse = new ImageBrowser5D(fileName);
	this->setCentralWidget(browse);
	//browse->show();

	// Disable the menu items for editing
	mergeAction->setEnabled(false);
	deleteAction->setEnabled(false);
	splitAction->setEnabled(false);

}

//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void NucleusEditor::CreateNewPlotWindow(void)
{
	if(!currentModel)
		return;

	pltWin.push_back(new PlotWindow(currentModel->GetSelectionModel()));
	pltWin.back()->show();
}

//******************************************************************************
// Create a new table window
//******************************************************************************
void NucleusEditor::CreateNewTableWindow(void)
{
	if(!currentModel)
		return;

	tblWin.push_back(new TableWindow(currentModel->GetSelectionModel()));
	tblWin.back()->ResizeToOptimalSize();
	tblWin.back()->show();
}

//******************************************************************************
// Open a python interpreter.  Ask the user where one is located if this
// information hasn't been specified yet.
//******************************************************************************
void NucleusEditor::OpenPythonWindow()
{
/* this code is still in progress...
	if(this->settings->contains("farsight/installpath") == false)
	  {
		this->settings->setValue("farsight/installpath","");
	  }
	QString installPath = this->settings->value("farsight/installpath").toString();
	QString quote = QString("\"");
	while(!QDir(installPath).exists() || !QDir(installPath + "/python").exists()
        || !QDir(installPath+ "/bin").exists())
	  {
		//Need to browse for the farsight install folder
		installPath = QFileDialog::getExistingDirectory(this, tr("Select Farsight installation directory"));
		if(installPath == "")
		  {
			//user pressed cancel, bail out now.
			return;
		  }
		this->settings->setValue("farsight/installpath",installPath);
	  }

	this->pythonFiles = installPath + QString("/python");
	this->exeFiles = installPath + QString("/bin");

  this->PythonDialog = new QtPythonDialog(this, this->argv0);
  QObject::connect(this->PythonDialog, SIGNAL(interpreterInitialized()),
                   this, SLOT(initPythonInterpretor()));
  this->PythonDialog->initializeInterpretor();
  this->PythonDialog->show();
  this->PythonDialog->raise();
  this->PythonDialog->activateWindow();
*/
	if(this->settings->contains("python/window") == false)
    {
		if(this->BrowseForPythonExecutable() == false)
		{
		//user cancelled operation, abort
		return;
		}
    }

	if(this->settings->contains("farsight/installpath") == false)
	{
		this->settings->setValue("farsight/installpath","");
	}

	QString cmd = this->settings->value("python/window").toString();
	QString installPath = this->settings->value("farsight/installpath").toString();
	QString quote = QString("\"");

	while(!QDir(installPath).exists() || !QDir(installPath + "/python").exists() || !QDir(installPath+ "/bin").exists())
	{
		//Need to browse for the farsight install folder folder
		installPath = QFileDialog::getExistingDirectory(this, tr("Select directory containing Farsight python scripts"));
		if(installPath == "")
		{
			//user pressed cancel, bail out now.
			return;
		}
		this->settings->setValue("farsight/installpath",installPath);
	}

	QString pythonFiles = installPath + QString("/python");
	QString exeFiles = installPath + QString("/bin");
	QString saxonFile = installPath + QString("/bin/saxon9.jar");

	QString path1Cmd = QString("import sys;sys.path.append('") + pythonFiles + QString("');");
	QString path2Cmd = QString("import os;os.environ['PATH'] = os.environ['PATH'] + ';") + exeFiles + QString("';");
	QString path3Cmd = QString("os.environ['PATH'] = os.environ['PATH'] + ';") + pythonFiles + QString("';");
	QString path4Cmd = QString("os.environ['CLASSPATH'] = ';") + saxonFile + QString("';");
	QString importCmds = QString("from farsightutils import *;");
	QString printCmd = QString("print 'FARSIGHT ENVIRONMENT';");
	QString arg = QString(" -i -c ") + quote + path1Cmd + path2Cmd + path3Cmd + path4Cmd + importCmds + printCmd + quote;
	cmd.append(arg);
  
	this->pythonProcess->startDetached(cmd);
}

//******************************************************************************
// Ask the user to select a python interpreter to use with Farsight
//******************************************************************************
bool NucleusEditor::BrowseForPythonExecutable()
  {
  if(this->pythonProcess->state() != QProcess::NotRunning)
    {
    if(this->ConfirmClosePython() == false)
      {
      //we can't change farsight's associated python program command until the
      //user is ready to shut it down...
      return false;
      }
    this->pythonProcess->terminate();
    }
  QString pythonExecutable = QFileDialog::getOpenFileName(this, tr("Open Python"));
  if(pythonExecutable == "")
    {
    //user pressed cancel, bail out now.
    return false;
    }
  this->settings->setValue("python/window", pythonExecutable);
  this->currentPythonLabel->setText(pythonExecutable);
  return true;
  }

//******************************************************************************
// Ask the user for permission to shut down their Python process.
//******************************************************************************
bool NucleusEditor::ConfirmClosePython()
  {
   QMessageBox msgBox;
   msgBox.setText(tr("Python shutdown required to continue"));
   msgBox.setInformativeText(tr("OK to shutdown Farsight's python process?"));
   msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
   msgBox.setDefaultButton(QMessageBox::Cancel);
   if(msgBox.exec() == QMessageBox::Cancel)
     {
     return false;
     }
   return true;
  }

