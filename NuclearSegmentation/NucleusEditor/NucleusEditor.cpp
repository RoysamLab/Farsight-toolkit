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
	createSegmentToolBar();

	setWindowTitle(tr("FARSIGHT: Nuclear Segmentation Tool"));

	seg = NULL;
	segWin = new SegmentationWindow();
	currentModel = NULL;

	lastPath = ".";
	myImgName = "";

	tblWin.clear();
	pltWin.clear();
	hisWin=NULL;
	
	//DEMO
	this->pythonProcess = new QProcess();
	this->settings = new QSettings("RPI", "Farsight");

	loadThread = NULL;
	binaryThread = NULL;
	seedThread = NULL;
	clusterThread = NULL;
	finalizeThread = NULL;
	featuresThread = NULL;

	this->resize(500,500);
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );	
}

NucleusEditor::~NucleusEditor()
{
	if(seg)
		delete seg;

	if(currentModel)
		delete currentModel;

	delete segWin;
}

//******************************************************************************
// Here we just show a message in the status bar when loading
//******************************************************************************
void NucleusEditor::createStatusBar()
{
    statusLabel = new QLabel(tr(" Ready"));
    statusBar()->addWidget(statusLabel, 1);
}

//** ToolBar Setup
void NucleusEditor::createSegmentToolBar()
{
	segmentTool = new QToolBar();

	segmentAbort = new QAction(tr("ABORT"),this);
	connect(segmentAbort, SIGNAL(triggered()), this, SLOT(abortSegment()));
	segmentTool->addAction(segmentAbort);

	segmentContinue = new QAction(tr("-->"),this);
	segmentContinue->setToolTip(tr("Continue Segmentation"));
	segmentContinue->setEnabled(false);
	connect(segmentContinue, SIGNAL(triggered()),this, SLOT(segment()));
	segmentTool->addAction(segmentContinue);
	
    segmentProgress = new QProgressBar();
	segmentProgress->setRange(0,7);
	segmentProgress->setTextVisible(true);
	segmentTool->addWidget(segmentProgress); 

	segmentTaskLabel = new QLabel;
	segmentTaskLabel->setFixedWidth(80);
	segmentTool->addWidget(segmentTaskLabel);
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

	//segmentAction = new QAction(tr("Start Segmentation Wizard..."), this);
	//segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation Wizard"));
	segmentAction = new QAction(tr("Start Segmentation..."), this);
	segmentAction->setEnabled(false);
	segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation on this Image"));
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

	saveAction = new QAction(tr("Save Result"), this);
	saveAction->setEnabled(false);
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
	newScatterAction->setEnabled(false);
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
	connect(mergeAction, SIGNAL(triggered()), this, SLOT(mergeCells()));	
	editMenu->addAction(mergeAction);

	editMenu->addSeparator();

	deleteAction = new QAction(tr("Delete Cells"), this);
	deleteAction->setStatusTip(tr("Deletes the selected cells"));
	connect(deleteAction,SIGNAL(triggered()),this,SLOT(deleteCells()));	
	editMenu->addAction(deleteAction);

	editMenu->addSeparator();

	// Splitting has two modes and therefore has two submenu items:
	// Start Splitting
	// End Splitting

	//Main Splitting Menu
	splitMenu= editMenu->addMenu(tr("&Splitting"));

	//submenu1:
	splitStartAction = new QAction(tr("Start Splitting"), this);
	splitStartAction->setStatusTip(tr("Start Splitting Mode"));
	connect(splitStartAction,SIGNAL(triggered()),this,SLOT(startSplitting()));
	splitMenu->addAction(splitStartAction);

	//submenu2:
	splitEndAction = new QAction(tr("End Splitting"), this);
	splitEndAction->setStatusTip(tr("End Splitting Mode"));
	connect(splitEndAction,SIGNAL(triggered()),this,SLOT(endSplitting()));
	splitMenu->addAction(splitEndAction);

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);

	viewMenu->setEnabled(true);
	setEditsEnabled(false);

}

void NucleusEditor::setEditsEnabled(bool val)
{
	editMenu->setEnabled(val);
	mergeAction->setEnabled(val);
	deleteAction->setEnabled(val);
	splitMenu->setEnabled(val);
	splitStartAction->setEnabled(val);
	splitEndAction->setEnabled(val);
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
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
	//Stop any running threads:
	this->abortSegment();

	//Save changes:
	if(seg)
	{
		if(seg->editsNotSaved)
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
	if(seg)
	{
		//QString name = QFileDialog::getSaveFileName(this, tr("Save File As"), lastPath, tr("XML (*.xml)") );
		//if(name == "")
		//{
		//	return false;
		//}
		//else
		//{
		//	seg->WriteToXML( name.toStdString() );
		//}

		if(seg->editsNotSaved)
		{
			std::string name = seg->GetDataFilename();
			name.erase(name.find_first_of("."));
			name.append(".xml");
			seg->SaveLabel();
			seg->WriteToXML( name );
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

	if(seg) {
		currentModel = new SegmentationModel(seg);
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
	if(seg)
		delete seg;
	seg = new ftk::NuclearSegmentation();
	if ( !seg->RestoreFromXML(filename.toStdString()) )
	{
		std::cerr << seg->GetErrorMessage() << std::endl;
		return;
	}

	//Now I have objects stored in memory - put the features into the model
	newModel();
	CreateNewTableWindow();
	CreateNewPlotWindow();
	
	hisWin = new HistoWindow(currentModel->GetSelectionModel());
	hisWin->show();

	segWin->SetModels(currentModel);
	segWin->SetChannelImage(seg->getDataImage());
	segWin->SetLabelImage(seg->getLabelImage());
	if( this->centralWidget() != this->segWin )
		this->setCentralWidget(segWin);

	this->update();

	// Enable the menu items for editing
	setEditsEnabled(true);
	newScatterAction->setEnabled(true);
	saveAction->setEnabled(true);

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

void NucleusEditor::startSplitting(void)
{	
	currentModel->startSplitTrigger();
}

void NucleusEditor::endSplitting(void)
{
	currentModel->endSplitTrigger();
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

	if(myImgName == "")
	{
		return;
	}

	if(currentModel)
		clearModel();

	QString dataFile = lastPath + "/" + myImgName;
	QString paramFile = lastPath + "/" + "Seg_Params.ini";
	if(seg) delete seg;
	seg = new ftk::NuclearSegmentation();
	if( !seg->SetInputs( dataFile.toStdString(), paramFile.toStdString() ) )
	{
		delete seg;
		seg = NULL;
		return;
	}
	
	int numChannels = myImg->GetImageInfo()->numChannels;
	if(numChannels > 1)
	{
		//I need to select a channel to use for segmentation
		int useChannel = 0;
		for (int i=0; i<numChannels; ++i)
		{
			std::string name = myImg->GetImageInfo()->channelNames.at(i);
			if(name.find("Nuclei") != string::npos)
				break;
			useChannel++;
		}
		if(useChannel >= numChannels)
		{
			//didn't find the string, so should probably ask which one to use
			useChannel = 0;
		}
		seg->SetChannel(useChannel);
	}
	
	this->addToolBar(Qt::TopToolBarArea, segmentTool);
	segmentContinue->setEnabled(false);
	segmentTool->setVisible(true);
	
	segmentState = 0;
	segment();

	//NuclearSegmentationWizard *wizard = new NuclearSegmentationWizard();
	//wizard->setParent(NULL);
	//this->setCentralWidget(wizard);
	//wizard->show();
}

void NucleusEditor::abortSegment()
{
	if(loadThread)
	{
		loadThread->wait();
		delete loadThread;
		loadThread = NULL;
	}
	if(binaryThread)
	{
		binaryThread->wait();
		delete binaryThread;
		binaryThread = NULL;
	}
	if(seedThread)
	{
		seedThread->wait();
		delete seedThread;
		seedThread = NULL;
	}
	if(clusterThread)
	{
		clusterThread->wait();
		delete clusterThread;
		clusterThread = NULL;
	}
	if(finalizeThread)
	{
		finalizeThread->wait();
		delete finalizeThread;
		finalizeThread = NULL;
	}
	if(featuresThread)
	{
		featuresThread->wait();
		delete featuresThread;
		featuresThread = NULL;
	}

	this->removeToolBar(segmentTool);
	segmentState = -1;
	//if(currentModel)
	//	clearModel();
	//delete seg;
	//seg = NULL;

	segWin->SetLabelImage(NULL);
	QApplication::restoreOverrideCursor();
	fileMenu->setEnabled(true);
	viewMenu->setEnabled(true);
	this->setEditsEnabled(false);
}

void NucleusEditor::segment()
{
	switch(segmentState)
	{
	case 0:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		// Disable the menu items for editing
		fileMenu->setEnabled(false);
		editMenu->setEnabled(false);
		viewMenu->setEnabled(false);

		segmentTaskLabel->setText(tr(" Loading "));
		segmentProgress->setValue(0);
		loadThread = new Load(seg);
		connect(loadThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 1;
		loadThread->start();
		break;
	case 1:
		if(loadThread)
		{
			delete loadThread;
			loadThread = NULL;
		}
		segmentTaskLabel->setText(tr(" Binarizing "));
		segmentProgress->setValue(1);
		binaryThread = new Binarize(seg);
		connect(binaryThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 2;
		binaryThread->start();
		break;
	case 2:
		if(binaryThread)
		{
			delete binaryThread;
			binaryThread = NULL;
		}
		segmentProgress->setValue(2);
		segmentTaskLabel->setText(tr(" Seeds "));
		seedThread = new SeedDetect(seg);
		connect(seedThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 3;
		seedThread->start();
		break;
	case 3:
		if(seedThread)
		{
			delete seedThread;
			seedThread = NULL;
		}
		segmentProgress->setValue(3);
		segmentTaskLabel->setText(tr(" Clustering "));
		clusterThread = new Cluster(seg);
		connect(clusterThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 4;
		clusterThread->start();
		break;
	case 4:	//Final State we get to after successful clustering - for seed editing
		if(clusterThread)
		{
			delete clusterThread;
			clusterThread = NULL;
		}
		segmentProgress->setValue(4);
		segWin->SetLabelImage(seg->getLabelImage());
		segmentTaskLabel->setText(tr(" Editing Seeds "));

		QApplication::restoreOverrideCursor();
		fileMenu->setEnabled(true);
		this->setEditsEnabled(false);
		viewMenu->setEnabled(true);
		segmentContinue->setEnabled(true);
		segmentState = 5;
		break;
	case 5:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		// Disable the menu items for editing
		fileMenu->setEnabled(false);
		this->setEditsEnabled(false);
		viewMenu->setEnabled(false);
		segmentContinue->setEnabled(false);

		segmentTaskLabel->setText(tr(" Finalizing "));
		segmentProgress->setValue(5);
		finalizeThread = new Finalize(seg);
		connect(finalizeThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 6;
		finalizeThread->start();
		break;
	case 6:
		if(finalizeThread)
		{
			delete finalizeThread;
			finalizeThread = NULL;
		}
		segmentProgress->setValue(6);
		segWin->SetLabelImage(seg->getLabelImage());
		segmentTaskLabel->setText(tr(" Features "));

		featuresThread = new Features(seg);
		connect(featuresThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 7;
		featuresThread->start();
		break;
	case 7:
		if(featuresThread)
		{
			delete featuresThread;
			featuresThread = NULL;
		}
		segmentProgress->setValue(7);
		segmentTaskLabel->setText(tr(" DONE "));

		newModel();
		CreateNewTableWindow();
		CreateNewPlotWindow();
		hisWin = new HistoWindow(currentModel->GetSelectionModel());
		hisWin->show();
		segWin->SetModels(currentModel);
		//this->update();

		QApplication::restoreOverrideCursor();
		fileMenu->setEnabled(true);
		saveAction->setEnabled(true);
		newScatterAction->setEnabled(true);
		this->setEditsEnabled(true);
		viewMenu->setEnabled(true);
		segmentState = -1;

		//Now remove the toolbar:
		this->removeToolBar(segmentTool);

		break;
	}

}


void NucleusEditor::loadImage()
{
	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							    "All Files (*.*)"));

    if(fileName == "")
		return;

	if(currentModel)
		clearModel();

	lastPath = QFileInfo(fileName).absolutePath();
	myImgName = QFileInfo(fileName).fileName();

	//******************************************************

	// NEW BROWSER:
	//ImageBrowser5D *browse = new ImageBrowser5D(fileName);
	//this->setCentralWidget(browse);
	
	// OLD BROWSER:
	myImg = ftk::Image::New();
	myImg->LoadFile(fileName.toStdString());
	segWin->SetChannelImage(myImg);
	if( this->centralWidget() != this->segWin )
		this->setCentralWidget(segWin);
	this->update();

	//************************************************

	// Disable the menu items for editing
	this->setEditsEnabled(false);
	segmentAction->setEnabled(true);
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

Load::Load(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;
}

void Load::run()
{
	mySeg->LoadData();
}

Binarize::Binarize(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;	
}

void Binarize::run()
{
	mySeg->Binarize(false);
}

SeedDetect::SeedDetect(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;	
}

void SeedDetect::run()
{
	mySeg->DetectSeeds(false);
}

Cluster::Cluster(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;	
}

void Cluster::run()
{
	mySeg->RunClustering();
}


Finalize::Finalize(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;
}

void Finalize::run()
{
	mySeg->Finalize();
}

Features::Features(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;
}

void Features::run()
{
	mySeg->LabelsToObjects();
	mySeg->ReleaseSegMemory();
}