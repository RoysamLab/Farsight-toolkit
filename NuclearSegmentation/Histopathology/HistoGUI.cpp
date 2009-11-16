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

#include "HistoGUI.h"
//*******************************************************************************
// HistoGUI
//
// This widget is a main window for Histopathology Project.  
//********************************************************************************

HistoGUI::HistoGUI(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{
	segView = new LabelImageViewQT();
	connect(segView, SIGNAL(mouseAt(int,int,int)), this, SLOT(setMouseStatus(int,int,int)));
	this->setCentralWidget(segView);

	createMenus();
	createStatusBar();

	setWindowTitle(tr("FARSIGHT: Histopathology"));

	lastPath = ".";
	myImgName = "";

	tblWin.clear();
	pltWin.clear();
	hisWin=NULL;

	nucSeg = NULL;

	table = NULL;
	selection = NULL;

	this->resize(500,500);
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );	
}

HistoGUI::~HistoGUI()
{
	if(nucSeg) delete nucSeg;
	if(selection) delete selection;
}

//******************************************************************************
// Here we just show a message in the status bar when loading
//******************************************************************************
void HistoGUI::createStatusBar()
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
void HistoGUI::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));

	loadAction = new QAction(tr("Load Image..."), this);
	loadAction->setStatusTip(tr("Load an image into the 5D image browser"));
	connect(loadAction, SIGNAL(triggered()), this, SLOT(loadImage()));
	fileMenu->addAction(loadAction);

	fileMenu->addSeparator();

	//segmentAction = new QAction(tr("Start Segmentation..."), this);
	//segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation on this Image"));
	//connect(segmentAction,SIGNAL(triggered()),this,SLOT(segmentImage()));
	//fileMenu->addAction(segmentAction);

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
	viewMenu->addAction(showHistoAction);

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

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);

	setMenusForResult(false);
}

void HistoGUI::setEditsEnabled(bool val)
{
	clearSelectAction->setEnabled(val);
	visitAction->setEnabled(val);
	addAction->setEnabled(val);
	mergeAction->setEnabled(val);
	deleteAction->setEnabled(val);
	splitZAction->setEnabled(val);
	splitAction->setEnabled(val);
	classAction->setEnabled(val);
}

void HistoGUI::setMenusForResult(bool val)
{
	saveAction->setEnabled(val);

	showBoundsAction->setEnabled(val);
	showIDsAction->setEnabled(val);
	newScatterAction->setEnabled(val);
	showHistoAction->setEnabled(val);
	imageIntensityAction->setEnabled(val);

	this->setEditsEnabled(val);
}

void HistoGUI::menusEnabled(bool val)
{
	fileMenu->setEnabled(val);
	viewMenu->setEnabled(val);
	editMenu->setEnabled(val);
}

//****************************************************************************
// SLOT: about()
//   A brief message about Farsight is displayed
//****************************************************************************
void HistoGUI::about()
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
void HistoGUI::setMouseStatus(int x, int y, int z)
{
	(this->statusLabel)->setText(QString::number(x) + ", " + QString::number(y) + ", " + QString::number(z));
}

//******************************************************************************
//Reimplement closeEvent to also close all other windows in the application
//******************************************************************************
void HistoGUI::closeEvent(QCloseEvent *event)
{
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

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

bool HistoGUI::saveResult()
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

void HistoGUI::clearSelections()
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
void HistoGUI::loadResult(void)
{
	QString filename  = QFileDialog::getOpenFileName(this,"Choose a Result",lastPath, 
			tr("XML Files (*.xml)\n"));

    if(filename == "")
		return;

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
}

void HistoGUI::loadImage()
{
	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							    "All Files (*.*)"));

    if(fileName == "")
		return;

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
	imageIntensityAction->setEnabled(true);
}

//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void HistoGUI::CreateNewPlotWindow(void)
{
	if(!table) return;

	pltWin.push_back(new PlotWindow());
	pltWin.back()->setModels(table,selection);
	pltWin.back()->show();
}

//******************************************************************************
// Create a new table window
//******************************************************************************
void HistoGUI::CreateNewTableWindow(void)
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
void HistoGUI::CreateNewHistoWindow(void)
{
	if(!table) return;

	if(this->hisWin)
		delete hisWin;
	hisWin = new HistoWindow();
	hisWin->setModels(table,selection);
	hisWin->show();
}

void HistoGUI::ShowHistogram(void)
{
	if(this->hisWin)
		hisWin->show();
	else
		this->CreateNewHistoWindow();
}

void HistoGUI::toggleBounds(void)
{
	if(!segView) return;

	if( showBoundsAction->isChecked() )
		segView->SetBoundsVisible(true);
	else
		segView->SetBoundsVisible(false);
}

void HistoGUI::toggleIDs(void)
{
	//if(!segView) return;

	//if( showIDsAction->isChecked() )
	//	segView->SetIDsVisible(true);
	//else
	//	segView->SetIDsVisible(false);
}

void HistoGUI::changeClass(void)
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

void HistoGUI::markVisited(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	
	nucSeg->MarkAsVisited(ids,1);
	this->updateViews();
}

void HistoGUI::addCell(int x1, int y1, int x2, int y2, int z)
{
	if(!nucSeg) return;
	int id = nucSeg->AddObject(x1, y1, z, x2, y2, z);
	if(id != 0)
	{
		this->updateViews();
		selection->select(id);
	}
}

void HistoGUI::deleteCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	nucSeg->Delete(ids);
	selection->clear();
	this->updateViews();
}

void HistoGUI::mergeCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	nucSeg->Merge(ids);
	selection->clear();
	this->updateViews();
}

void HistoGUI::splitCell(int x1, int y1, int z1, int x2, int y2, int z2)
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

void HistoGUI::splitCellAlongZ(void)
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

//Call this slot when the table has been modified (new rows or columns) to update the views:
void HistoGUI::updateViews()
{
	for(unsigned int p=0; p<pltWin.size(); ++p)
		pltWin.at(p)->update();
	
	for(unsigned int p=0; p<tblWin.size(); ++p)	
		tblWin.at(p)->update();

	if(segView)
		segView->update();
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************











