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
	segView = new SegmentationView();
	connect(segView, SIGNAL(mouseAt(int,int,int)), this, SLOT(setMouseStatus(int,int,int)));
	this->setCentralWidget(segView);

	createMenus();
	createStatusBar();
	createSegmentToolBar();

	setWindowTitle(tr("FARSIGHT: Nuclear Segmentation Tool"));

	lastPath = ".";
	myImgName = "";

	tblWin.clear();
	pltWin.clear();
	hisWin=NULL;
	pWizard=NULL;

	nucSeg = NULL;
	loadThread = NULL;
	binaryThread = NULL;
	seedThread = NULL;
	clusterThread = NULL;
	finalizeThread = NULL;
	featuresThread = NULL;
	table = NULL;
	selection = NULL;

	this->resize(500,500);
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );	
}

NucleusEditor::~NucleusEditor()
{
	if(nucSeg) delete nucSeg;
	if(selection) delete selection;
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

	segmentAbort = new QAction(QIcon(":/icons/abort.png"), tr("ABORT"), this);
	segmentAbort->setToolTip(tr("Abort Segmentation and save nothing"));
	segmentAbort->setEnabled(true);
	connect(segmentAbort, SIGNAL(triggered()), this, SLOT(abortSegment()));
	segmentTool->addAction(segmentAbort);

	segmentContinue = new QAction(QIcon(":/icons/go.png"), tr("GO"), this);
	segmentContinue->setToolTip(tr("Continue Segmentation"));
	segmentContinue->setEnabled(false);
	connect(segmentContinue, SIGNAL(triggered()),this, SLOT(segment()));
	segmentTool->addAction(segmentContinue);

	segmentStop = new QAction(QIcon(":/icons/end.png"), tr("STOP"), this);
	segmentStop->setToolTip(tr("Stop Segmentation (skip finalization), and use clustering as result"));
	segmentStop->setEnabled(false);
	connect(segmentStop, SIGNAL(triggered()), this, SLOT(stopSegment()));
	segmentTool->addAction(segmentStop);
	
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

	segmentAction = new QAction(tr("Start Segmentation..."), this);
	segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation on this Image"));
	connect(segmentAction,SIGNAL(triggered()),this,SLOT(segmentImage()));
	fileMenu->addAction(segmentAction);

	fileMenu->addSeparator();

	xmlAction = new QAction(tr("Load Result..."), this);
	xmlAction->setStatusTip(tr("Open an XML result file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu->addAction(xmlAction);

	saveAction = new QAction(tr("Save Result"), this);
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

	showBoundsAction = new QAction(tr("Show &Boundaries"), this);
	showBoundsAction->setCheckable(true);
	showBoundsAction->setChecked(false);
	showBoundsAction->setStatusTip(tr("Draw boundaries using a label image"));
	showBoundsAction->setShortcut(tr("Ctrl+B"));
	connect(showBoundsAction, SIGNAL(triggered()), this, SLOT(toggleBounds()));
	viewMenu->addAction(showBoundsAction);

	showIDsAction = new QAction(tr("Show &IDs"), this);
	showIDsAction->setCheckable(true);
	showIDsAction->setChecked(false);
	showIDsAction->setStatusTip(tr("Draw ID numbers at centroid locations"));
	showIDsAction->setShortcut(tr("Ctrl+I"));
	//connect(showIDsAction, SIGNAL(triggered()), this, SLOT(toggleIDs()));
	//viewMenu->addAction(showIDsAction);

	viewMenu->addSeparator();
	newScatterAction = new QAction(tr("New Scatter"), this);
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

	showHistoAction = new QAction(tr("Show Histogram"),this);
	showHistoAction->setStatusTip(tr("Show a Histogram"));
	connect(showHistoAction,SIGNAL(triggered()),this,SLOT(ShowHistogram()));
	//viewMenu->addAction(showHistoAction);

	imageIntensityAction = new QAction(tr("Adjust Image Intensity"), this);
	imageIntensityAction->setStatusTip(tr("Allows modification of image intensity"));
	imageIntensityAction->setShortcut(tr("Ctrl+G"));
	connect(imageIntensityAction, SIGNAL(triggered()), segView, SLOT(AdjustImageIntensity()));
	viewMenu->addAction(imageIntensityAction);

	//EDITING MENU	
	editMenu = menuBar()->addMenu(tr("&Editing"));

	clearSelectAction = new QAction(tr("Clear Selections"), this);
	clearSelectAction->setStatusTip(tr("Clear Current Object Selections"));
	clearSelectAction->setShortcut(tr("Ctrl+C"));
	connect(clearSelectAction, SIGNAL(triggered()), this, SLOT(clearSelections()));
	editMenu->addAction(clearSelectAction);

	visitAction = new QAction(tr("Mark as Visited"), this);
	visitAction->setStatusTip(tr("Mark Selected Objects as Visited"));
	visitAction->setShortcut(tr("Ctrl+V"));
	connect(visitAction, SIGNAL(triggered()), this, SLOT(markVisited()));
	editMenu->addAction(visitAction);

	editMenu->addSeparator();

	classAction = new QAction(tr("Change Class"), this);
	classAction->setStatusTip(tr("Modify the class designation for the selected objects"));
	classAction->setShortcut(tr("Ctrl+L"));
	connect(classAction, SIGNAL(triggered()), this, SLOT(changeClass()));
	editMenu->addAction(classAction);

	addAction = new QAction(tr("Add Cell"), this);
	addAction->setStatusTip(tr("Draw a Box to add a new cell"));
	addAction->setShortcut(tr("Ctrl+A"));
	connect(addAction, SIGNAL(triggered()), segView, SLOT(GetBox()));
	connect(segView, SIGNAL(boxDrawn(int,int,int,int,int)), this, SLOT(addCell(int,int,int,int,int)));
	editMenu->addAction(addAction);

	mergeAction = new QAction(tr("Merge Cells"), this);
	mergeAction->setStatusTip(tr("Merge Cells"));
	mergeAction->setShortcut(tr("Ctrl+M"));
	connect(mergeAction, SIGNAL(triggered()), this, SLOT(mergeCells()));	
	editMenu->addAction(mergeAction);

	deleteAction = new QAction(tr("Delete Cells"), this);
	deleteAction->setStatusTip(tr("Deletes the selected cells"));
	deleteAction->setShortcut(tr("Ctrl+D"));
	connect(deleteAction,SIGNAL(triggered()),this,SLOT(deleteCells()));	
	editMenu->addAction(deleteAction);

	splitZAction = new QAction(tr("Split Cell At Z"), this);
	splitZAction->setStatusTip(tr("Split selected cell along the current Z slice"));
	splitZAction->setShortcut(tr("Ctrl+T"));
	connect(splitZAction, SIGNAL(triggered()), this, SLOT(splitCellAlongZ()));
	editMenu->addAction(splitZAction);

	splitAction = new QAction(tr("Split Cell X-Y"), this);
	splitAction->setStatusTip(tr("Split a cell by choosing two seed points"));
	splitAction->setShortcut(tr("Ctrl+P"));
	connect(splitAction, SIGNAL(triggered()), segView, SLOT(Get2Points()));
	connect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));
	editMenu->addAction(splitAction);

	//editMenu->addSeparator();

	exclusionAction = new QAction(tr("Apply Exclusion Margin..."), this);
	exclusionAction->setStatusTip(tr("Set parameters for exclusion margin"));
	connect(exclusionAction, SIGNAL(triggered()), this, SLOT(applyExclusionMargin()));
	//editMenu->addAction(exclusionAction);

	//TOOL MENU
	toolMenu = menuBar()->addMenu(tr("Tools"));
	svmAction = new QAction(tr("Detect Outliers"), this);
	connect(svmAction, SIGNAL(triggered()), this, SLOT(startSVM()));
	toolMenu->addAction(svmAction);

	kplsAction = new QAction(tr("Classify"), this);
	connect(kplsAction, SIGNAL(triggered()), this, SLOT(startKPLS()));
	toolMenu->addAction(kplsAction);

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);

	setMenusForResult(false);
	segmentAction->setEnabled(false);
}

void NucleusEditor::setEditsEnabled(bool val)
{
	clearSelectAction->setEnabled(val);
	visitAction->setEnabled(val);
	addAction->setEnabled(val);
	mergeAction->setEnabled(val);
	deleteAction->setEnabled(val);
	splitZAction->setEnabled(val);
	splitAction->setEnabled(val);
	classAction->setEnabled(val);
	exclusionAction->setEnabled(val);
}

void NucleusEditor::setMenusForResult(bool val)
{
	saveAction->setEnabled(val);

	showBoundsAction->setEnabled(val);
	showIDsAction->setEnabled(val);
	newScatterAction->setEnabled(val);
	showHistoAction->setEnabled(val);
	imageIntensityAction->setEnabled(val);

	svmAction->setEnabled(val);
	kplsAction->setEnabled(val);

	this->setEditsEnabled(val);
}

void NucleusEditor::menusEnabled(bool val)
{
	fileMenu->setEnabled(val);
	viewMenu->setEnabled(val);
	editMenu->setEnabled(val);
	toolMenu->setEnabled(val);
}

//****************************************************************************
// SLOT: about()
//   A brief message about Farsight is displayed
//****************************************************************************
void NucleusEditor::about()
{
	    QMessageBox::about(this, tr("About FARSIGHT"),
            tr("<h2>FARSIGHT</h2>"
			   "<h3>Rensselear Polytechnic Institute</h3>"
			   "<a><u>http://www.farsight-toolkit.org</a></u>"
               ));
}

//******************************************************************************
// SLOT: changes the status bar to say the mouse coordinates
//******************************************************************************
void NucleusEditor::setMouseStatus(int x, int y, int z)
{
	(this->statusLabel)->setText(QString::number(x) + ", " + QString::number(y) + ", " + QString::number(z));
}

//******************************************************************************
//Reimplement closeEvent to also close all other windows in the application
//******************************************************************************
void NucleusEditor::closeEvent(QCloseEvent *event)
{
	//Save changes:
	if( !checkSaveSeg() )
	{
		event->ignore();
		return;
	}

	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

	//Stop any running threads:
	this->abortSegment();

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

//return true if ask to save and saved, or do not want to save, false if want to save and not saved
bool NucleusEditor::checkSaveSeg()
{
	if(nucSeg)
	{
		if(nucSeg->EditsNotSaved())
		{
			QString msg = tr("Recent Edits not saved do you want to save them before exiting?");
			QMessageBox::StandardButton button = QMessageBox::information ( 0, tr("Exit"), \
				msg, QMessageBox::Yes | QMessageBox::No , QMessageBox::NoButton );

			if(button == QMessageBox::Yes)
			{
				if(	!this->saveResult()	)
				{
					return false;
				}
			}
		}
	}
	return true;
}

bool NucleusEditor::saveResult()
{
	if(nucSeg)
	{
		if(nucSeg->EditsNotSaved())
		{
			std::string name = nucSeg->GetDataFilename();
			name.erase(name.find_first_of("."));
			name.append(".xml");
			nucSeg->SaveChanges(name);
		}
		
	}
	return true;
}

void NucleusEditor::clearSelections()
{
	if(selection)
	{
		selection->clear();
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
	if( !checkSaveSeg() )
		return;

	QString filename  = QFileDialog::getOpenFileName(this,"Choose a Result",lastPath, 
			tr("XML Files (*.xml)\n"));

    if(filename == "")
		return;

	abortSegment();

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();

	lastPath = path;

	//segResult = new ftk::NuclearSegmentation();
	if(nucSeg) delete nucSeg;
	nucSeg = new ftk::NuclearSegmentation();
	if ( !nucSeg->LoadAll(filename.toStdString()) )
	{
		std::cerr << nucSeg->GetErrorMessage() << std::endl;
		return;
	}

	table = nucSeg->GetFeatureTable();

	if(selection) delete selection;
	selection = new ObjectSelection();

	segView->SetChannelImage(nucSeg->GetDataImage());
	segView->SetLabelImage(nucSeg->GetLabelImage(), selection);

	CreateNewTableWindow();
	CreateNewPlotWindow();

	// Enable the menu items for editing
	setMenusForResult(true);
	segmentAction->setEnabled(true);
}

void NucleusEditor::loadImage()
{
	if( !checkSaveSeg() )
		return;

	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							    "All Files (*.*)"));

    if(fileName == "")
		return;

	abortSegment();

	lastPath = QFileInfo(fileName).absolutePath();
	myImgName = QFileInfo(fileName).fileName();

	//******************************************************
	// NEW BROWSER:
	//ImageBrowser5D *browse = new ImageBrowser5D(fileName);
	//this->setCentralWidget(browse);
	
	// OLD BROWSER:
	myImg = ftk::Image::New();
	myImg->LoadFile(fileName.toStdString());
	segView->SetChannelImage(myImg);

	// Disable the menu items for editing
	this->setMenusForResult(false);
	segmentAction->setEnabled(true);
	imageIntensityAction->setEnabled(true);
}

//**********************************************************************
// SLOT: start the pattern analysis widget:
//**********************************************************************
void NucleusEditor::startSVM()
{
	if(!table) return;

	if(pWizard)
	{
		delete pWizard;
	}
	pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_SVM, "", "outlier?", this);
	connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
	pWizard->show();
}

//**********************************************************************
// SLOT: start the pattern analysis widget:
//**********************************************************************
void NucleusEditor::startKPLS()
{
	if(!table) return;

	if(pWizard)
	{
		delete pWizard;
	}
	pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_KPLS, "class", "prediction", this);
	connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
	pWizard->show();
}

//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void NucleusEditor::CreateNewPlotWindow(void)
{
	if(!table) return;

	pltWin.push_back(new PlotWindow());
	pltWin.back()->setModels(table,selection);
	pltWin.back()->show();
}

//******************************************************************************
// Create a new table window
//******************************************************************************
void NucleusEditor::CreateNewTableWindow(void)
{
	if(!table) return;

	tblWin.push_back(new TableWindow());
	tblWin.back()->setModels(table,selection);
	tblWin.back()->show();

	//vtkQtTableView * view = vtkQtTableView::New();
	//view->AddRepresentationFromInput(table);
	//view->Update();
	//view->GetWidget()->show();
}

//*******************************************************************************
// Create new Histogram Window
//*******************************************************************************
void NucleusEditor::CreateNewHistoWindow(void)
{
	//if(this->hisWin)
	//	delete hisWin;
	//hisWin = new HistoWindow(currentModel->GetSelectionModel());
	//hisWin->show();
}

void NucleusEditor::ShowHistogram(void)
{
	if(this->hisWin)
		hisWin->show();
	else
		this->CreateNewHistoWindow();
}

void NucleusEditor::toggleBounds(void)
{
	if(!segView) return;

	if( showBoundsAction->isChecked() )
		segView->SetBoundsVisible(true);
	else
		segView->SetBoundsVisible(false);
}

void NucleusEditor::toggleIDs(void)
{
	//if(!segView) return;

	//if( showIDsAction->isChecked() )
	//	segView->SetIDsVisible(true);
	//else
	//	segView->SetIDsVisible(false);
}

void NucleusEditor::changeClass(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());

	//Get the new class number:
	bool ok;
	QString msg = tr("Change the class of all selected items to: \n");
	int newClass = QInputDialog::getInteger(NULL, tr("Change Class"), msg, -1, -1, 10, 1, &ok);
	
	//Change the class of these objects:
	if(ok)
	{
		nucSeg->SetClass(ids,newClass);
		this->updateViews();
	}
}

void NucleusEditor::markVisited(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	
	nucSeg->MarkAsVisited(ids,1);
	this->updateViews();
}

void NucleusEditor::addCell(int x1, int y1, int x2, int y2, int z)
{
	if(!nucSeg) return;
	int id = nucSeg->AddObject(x1, y1, z, x2, y2, z);
	if(id != 0)
	{
		this->updateViews();
		selection->select(id);
	}
}

void NucleusEditor::deleteCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	nucSeg->Delete(ids);
	selection->clear();
	this->updateViews();
}

void NucleusEditor::mergeCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	nucSeg->Merge(ids);
	selection->clear();
	this->updateViews();
}

void NucleusEditor::splitCell(int x1, int y1, int z1, int x2, int y2, int z2)
{
	if(!nucSeg) return;

	ftk::Object::Point P1;
	ftk::Object::Point P2;
	P1.t=0;
	P1.x=x1;
	P1.y=y1;
	P1.z=z1;
	P2.t=0;
	P2.x=x2;
	P2.y=y2;
	P2.z=z2;

	nucSeg->Split(P1, P2);
	selection->clear();
	this->updateViews();
}

void NucleusEditor::splitCellAlongZ(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	selection->clear();
	for ( set<long int>::iterator it=sels.begin(); it != sels.end(); it++ )
	{
		nucSeg->SplitAlongZ(*it,segView->GetCurrentZ());
	}
	this->updateViews();
}

void NucleusEditor::applyExclusionMargin(void)
{
	//Get the parameters to use for the brick rule:
	int xy = 0; 
	int z = 0;
	MarginDialog *dialog = new MarginDialog(this);
	if( dialog->exec() )	
	{
		xy = dialog->getMargin();
		z = dialog->getZ();
	}
	delete dialog;

	//Now apply the brick rule to my image!!!!
	//if(currentModel)
	//	currentModel->applyMargins(xy,z);

}

void NucleusEditor::segmentImage()
{
	if(myImgName == "")
	{
		return;
	}

	QString dataFile = lastPath + "/" + myImgName;

	QVector<QString> chs;
	int numChannels = myImg->GetImageInfo()->numChannels;
	for (int i=0; i<numChannels; ++i)
	{
		chs.push_back(QString::fromStdString(myImg->GetImageInfo()->channelNames.at(i)));
	}

	//Get the paramFile to use:
	QString paramFile = "";
	int useChannel;
	ParamsFileDialog *dialog = new ParamsFileDialog(lastPath,chs,this);
	if( dialog->exec() )	
	{
		paramFile = dialog->getFileName();
		useChannel = dialog->getChannelNumber();
	}
	delete dialog;

	if(nucSeg) delete nucSeg;
	nucSeg = new ftk::NuclearSegmentation();
	if( !nucSeg->SetInputs( dataFile.toStdString(), paramFile.toStdString() ) )
	{
		delete nucSeg;
		nucSeg = NULL;
		return;
	}
	
	nucSeg->SetChannel(useChannel);
	
	this->addToolBar(Qt::TopToolBarArea, segmentTool);
	segmentStop->setEnabled(false);
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

	quitNucSeg();

	QApplication::restoreOverrideCursor();
	menusEnabled(true);
	setMenusForResult(false);
}

//Call this slot when the table has been modified (new rows or columns) to update the views:
void NucleusEditor::updateViews()
{
	for(unsigned int p=0; p<pltWin.size(); ++p)
		pltWin.at(p)->update();
	
	for(unsigned int p=0; p<tblWin.size(); ++p)	
		tblWin.at(p)->update();

	if(segView)
		segView->update();
}

//***************************************************************************
// THIS FUNCTION CLEARS MEMORY BY CLOSING OPEN WINDOWS 
// THEN CLEARING THE SEGMENTATION
//***************************************************************************
void NucleusEditor::quitNucSeg(void)
{
	for(unsigned int p=0; p<pltWin.size(); ++p)
		if ((pltWin.at(p))->isVisible())
			(pltWin.at(p))->close();
	
	for(unsigned int p=0; p<tblWin.size(); ++p)	
		if ((tblWin.at(p))->isVisible())
			(tblWin.at(p))->close();

	pltWin.clear();
	tblWin.clear();

	//Close the histogram
	if (hisWin !=NULL) 
	{
		hisWin->close();
		hisWin=NULL;
	}

	//segView->SetChannelImage(NULL);
	segView->SetLabelImage(NULL);
	if(nucSeg)
	{
		delete nucSeg;
		nucSeg = NULL;
	}
}

//Go here when the "jump" button in clicked during nuclear segmentation:
void NucleusEditor::stopSegment(void)
{
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
	// Disable the menu items for editing
	menusEnabled(false);
	segmentStop->setEnabled(false);
	segmentContinue->setEnabled(false);

	segmentTaskLabel->setText(tr(" Skipping "));
	segmentProgress->setValue(5);
	segmentState = 6;
	this->segment();
}

void NucleusEditor::segment()
{
	switch(segmentState)
	{
	case 0:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		menusEnabled(false);

		segmentTaskLabel->setText(tr(" Loading "));
		segmentProgress->setValue(0);
		loadThread = new Load(nucSeg);
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
		binaryThread = new Binarize(nucSeg);
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
		seedThread = new SeedDetect(nucSeg);
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
		clusterThread = new Cluster(nucSeg);
		connect(clusterThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState = 4;
		clusterThread->start();
		break;
	case 4:	//Final State we get to after successful clustering - for checking parameters
		if(clusterThread)
		{
			delete clusterThread;
			clusterThread = NULL;
		}
		segmentProgress->setValue(4);
		segView->SetLabelImage(nucSeg->GetLabelImage());
		segmentTaskLabel->setText(tr(" Inspect "));

		QApplication::restoreOverrideCursor();
		segmentStop->setEnabled(true);
		segmentContinue->setEnabled(true);
		//I should pop-up a msg here:
		segmentState = 5;
		break;
	case 5:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		segmentStop->setEnabled(false);
		segmentContinue->setEnabled(false);

		segmentTaskLabel->setText(tr(" Finalizing "));
		segmentProgress->setValue(5);
		finalizeThread = new Finalize(nucSeg);
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
		if(selection) delete selection;
		selection = new ObjectSelection();
		segView->SetLabelImage(nucSeg->GetLabelImage(),selection);
		segmentTaskLabel->setText(tr(" Features "));

		featuresThread = new Features(nucSeg);
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
		
		table = nucSeg->GetFeatureTable();
		CreateNewTableWindow();
		CreateNewPlotWindow();

		QApplication::restoreOverrideCursor();
		menusEnabled(true);
		setMenusForResult(true);
		showBoundsAction->setChecked(true);
		showIDsAction->setChecked(true);

		segmentState = -1;

		//Now remove the toolbar:
		this->removeToolBar(segmentTool);

		break;
	}

}

//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
// A dialog for changing the exclusion margin used for this image:
//******************************************************************************************
MarginDialog::MarginDialog(QWidget *parent)
: QDialog(parent)
{
	QLabel * header = new QLabel(tr("Please set parameters for the Brick Rule:"));
	QLabel * mLabel = new QLabel(tr("XY Margin: "));
	marginSpin = new QSpinBox();
	marginSpin->setMinimum(0);
	marginSpin->setMaximum(100);
	QLabel * pixLabel = new QLabel(tr("pixels"));

	QLabel * zLabel = new QLabel(tr("Z Margin: "));
	zSpin = new QSpinBox();
	zSpin->setMinimum(0);
	zSpin->setMaximum(10);
	QLabel * slcLabel = new QLabel(tr("slices"));

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));

	QGridLayout * layout = new QGridLayout();
	layout->addWidget(header,0,0,1,3);
	layout->addWidget(mLabel,1,0,1,1);
	layout->addWidget(marginSpin,1,1,1,1);
	layout->addWidget(pixLabel,1,2,1,1);
	layout->addWidget(zLabel,2,0,1,1);
	layout->addWidget(zSpin,2,1,1,1);
	layout->addWidget(slcLabel,2,2,1,1);
	layout->addWidget(okButton,3,2,1,1);
	this->setLayout(layout);
	this->setWindowTitle(tr("Apply Exclusion Margin"));
}

int MarginDialog::getMargin()
{
	return marginSpin->value();
}

int MarginDialog::getZ()
{
	return zSpin->value();
}
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

//***********************************************************************************
//***********************************************************************************
// A dialog to get the paramaters file to use and specify the channel if image has
// more than one:
//***********************************************************************************
ParamsFileDialog::ParamsFileDialog(QString lastPth, QVector<QString> channels, QWidget *parent)
: QDialog(parent)
{
	this->lastPath = lastPth;

	channelLabel = new QLabel("Choose Channel: ");
	channelCombo = new QComboBox();
	for(int v = 0; v<channels.size(); ++v)
	{
			channelCombo->addItem(channels.at(v));
	}
	QHBoxLayout *chLayout = new QHBoxLayout;
	chLayout->addWidget(channelLabel);
	chLayout->addWidget(channelCombo);

	autoButton = new QRadioButton(tr("Automatic Parameter Selection"),this);
	autoButton->setChecked(true);

	fileButton = new QRadioButton(tr("Use Parameter File..."),this);

	fileCombo = new QComboBox();
	fileCombo->addItem(tr(""));
	fileCombo->addItem(tr("Browse..."));
	connect(fileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	QHBoxLayout *bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addLayout(chLayout);
	layout->addWidget(autoButton);
	layout->addWidget(fileButton);
	layout->addWidget(fileCombo);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Parameters"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}
QString ParamsFileDialog::getFileName()
{
	if(autoButton->isChecked())
	{
		return QString("");
	}
	else
	{
		return fileCombo->currentText();
	}
}

int ParamsFileDialog::getChannelNumber()
{
	return channelCombo->currentIndex();
}

void ParamsFileDialog::ParamBrowse(QString comboSelection)
{
	//First check to see if we have selected browse
	if (comboSelection != tr("Browse..."))
		return;

	QString newfilename  = QFileDialog::getOpenFileName(this,"Choose a Parameters File",lastPath, 
			tr("INI Files (*.ini)\n" 
			   "TXT Files (*.txt)\n" 
			   "All Files (*.*)\n"));

	if (newfilename == "")
	{
		fileCombo->setCurrentIndex(0);
		return;
	}

	if( newfilename == fileCombo->currentText() )
		return;

	fileCombo->setCurrentIndex(0);
	fileCombo->setItemText(0,newfilename);
}
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************


//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
// Threads for running the segmentation algorithm:
//***********************************************************************************
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
	mySeg->ComputeFeatures();
	mySeg->ReleaseSegMemory();
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************











