/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


#include "NucleusEditor.h"
//*******************************************************************************
// NucleusEditor
//
// This widget is a main window for editing nuclei.  
//********************************************************************************

NucleusEditor::NucleusEditor(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{
	segView = new LabelImageViewQT();
	connect(segView, SIGNAL(mouseAt(int,int,int)), this, SLOT(setMouseStatus(int,int,int)));
	this->setCentralWidget(segView);

	createMenus();
	createStatusBar();
	createSegmentToolBar();

	setWindowTitle(tr("FARSIGHT: Nuclear Segmentation Tool"));

	lastPath = ".";

	tblWin.clear();
	pltWin.clear();
	hisWin=NULL;
	pWizard=NULL;

	myImg = NULL;
	labImg = NULL;

	nucSeg = NULL;
	nucChannel = 0;
	nucsegThread = NULL;
	featuresThread = NULL;
	table = NULL;
	selection = NULL;

	this->resize(800,800);
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
	segmentTool = this->addToolBar(tr("Segment"));
	segmentTool->setVisible(false);

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
	segmentProgress->setRange(0,6);
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

	xmlAction = new QAction(tr("Load Result..."), this);
	xmlAction->setStatusTip(tr("Open an XML result file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu->addAction(xmlAction);

	saveAction = new QAction(tr("Save Result"), this);
	saveAction->setStatusTip(tr("Save a segmentation result image"));
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(saveResult()));
	fileMenu->addAction(saveAction);

	saveTableAction = new QAction(tr("Save Table"), this);
	saveTableAction->setStatusTip(tr("Save the features table"));
	connect(saveTableAction, SIGNAL(triggered()), this, SLOT(saveTable()));
	fileMenu->addAction(saveTableAction);

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

	//TOOL MENU
	toolMenu = menuBar()->addMenu(tr("Tools"));

	segmentAction = new QAction(tr("Nuclear Segmentation"), this);
	segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation on this Image"));
	connect(segmentAction,SIGNAL(triggered()),this,SLOT(segmentImage()));
	toolMenu->addAction(segmentAction);

	cytoAction = new QAction(tr("Cytoplasm Segmentation"), this);
	cytoAction->setStatusTip(tr("Starts the Cytoplasm Segmentation"));
	connect(cytoAction, SIGNAL(triggered()), this, SLOT(cytoSeg()));
	toolMenu->addAction(cytoAction);

	toolMenu->addSeparator();

	svmAction = new QAction(tr("Detect Outliers"), this);
	connect(svmAction, SIGNAL(triggered()), this, SLOT(startSVM()));
	toolMenu->addAction(svmAction);

	kplsAction = new QAction(tr("Classify"), this);
	connect(kplsAction, SIGNAL(triggered()), this, SLOT(startKPLS()));
	toolMenu->addAction(kplsAction);

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

	editMenu->addSeparator();

	exclusionAction = new QAction(tr("Apply Exclusion Margin..."), this);
	exclusionAction->setStatusTip(tr("Set parameters for exclusion margin"));
	connect(exclusionAction, SIGNAL(triggered()), this, SLOT(applyExclusionMargin()));
	editMenu->addAction(exclusionAction);

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);

	setMenusForResult(false);
	segmentAction->setEnabled(false);
	cytoAction->setEnabled(false);
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
		if(nucSeg->EditsNotSaved)
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
	if(!labImg)
		return false;

	int ch = requestChannel(labImg);
	
	if(ch==-1)
		return false;

	QString filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("TIFF Image (*.tif)"));

	if(filename == "")
		return false;

	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();
	lastPath = path;

	QString fullBase = path + "/" + name;
	return labImg->SaveChannelAs(ch,fullBase.toStdString() ,"tif");
}

bool NucleusEditor::saveTable()
{
	if(!table)
		return false;

	QString filename = QFileDialog::getSaveFileName(this, tr("Save As..."),lastPath, tr("TEXT(*.txt)"));

	if(filename == "")
		return false;

	lastPath = QFileInfo(filename).absolutePath();

	return ftk::SaveTable(filename.toStdString(), table);
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
	if ( 0 )//!nucSeg->LoadAll(filename.toStdString()) )
	{
		std::cerr << nucSeg->GetErrorMessage() << std::endl;
		return;
	}

	//table = nucSeg->GetFeatureTable();

	if(selection) delete selection;
	selection = new ObjectSelection();

	segView->SetChannelImage(nucSeg->GetDataImage());
	segView->SetLabelImage(nucSeg->GetLabelImage(), selection);

	CreateNewTableWindow();
	CreateNewPlotWindow();

	// Enable the menu items for editing
	setMenusForResult(true);
	segmentAction->setEnabled(true);
	cytoAction->setEnabled(true);
}

void NucleusEditor::loadImage()
{
	if( !checkSaveSeg() )
		return;

	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
								"XML Image Definition (*.xml)\n"
							    "All Files (*.*)"));
    if(fileName == "")
		return;

	abortSegment();

	lastPath = QFileInfo(fileName).absolutePath();
	QString myExt = QFileInfo(fileName).suffix();
	if(QFileInfo(fileName).suffix() == "xml")
	{
		if(!this->loadXMLImage(fileName.toStdString())) return;
	}
	else
	{
		myImg = ftk::Image::New();
		if(!myImg->LoadFile(fileName.toStdString()))
		{
			myImg = NULL;
			return;
		}
	}
	segView->SetChannelImage(myImg);

	// Disable the menu items for editing
	this->setMenusForResult(false);
	segmentAction->setEnabled(true);
	cytoAction->setEnabled(false);
	imageIntensityAction->setEnabled(true);
}

bool NucleusEditor::loadXMLImage(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Image" ) != 0 )
		return false;

	std::vector<std::string> files;
	std::vector<std::string> chName;
	std::vector<unsigned char> color;

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "file" ) == 0 )
		{
			files.push_back( parentElement->GetText() );
			chName.push_back( parentElement->Attribute("chname") );
			color.push_back( atoi(parentElement->Attribute("r")) );
			color.push_back( atoi(parentElement->Attribute("g")) );
			color.push_back( atoi(parentElement->Attribute("b")) );
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();

	myImg = ftk::Image::New();
	if(!myImg->LoadFilesAsMultipleChannels(files,chName,color))	//Load for display
	{
		myImg = NULL;
		return false;
	}
	return true;
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
	if(!table) return;

	if(this->hisWin)
		delete hisWin;
	hisWin = new HistoWindow();
	hisWin->setModels(table,selection);
	hisWin->show();
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
		//nucSeg->SetClass(ids,newClass);
		this->updateViews();
	}
}

void NucleusEditor::markVisited(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	
	//nucSeg->MarkAsVisited(ids,1);
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

	selection->clear();
	nucSeg->Exclude(xy, z);
	this->updateViews();
}

int NucleusEditor::requestChannel(ftk::Image::Pointer img)
{
	QStringList chs;
	int numChannels = img->GetImageInfo()->numChannels;
	for (int i=0; i<numChannels; ++i)
	{
		chs << QString::fromStdString(img->GetImageInfo()->channelNames.at(i));
	}

	bool ok=false;
	QString choice = QInputDialog::getItem(this, tr("Choose Channel"), tr("Channel:"),chs,0,false,&ok);
	if( !ok || choice.isEmpty() )
		return -1;

	int ch=0;
	for(int i=0; i<numChannels; ++i)
	{
		if( choice == QString::fromStdString(img->GetImageInfo()->channelNames.at(i)) )
		{
			ch = i;
			break;
		}	
	}
	return ch;
}

void NucleusEditor::cytoSeg(void)
{
	cytChannel = requestChannel(myImg);

	ftk::CytoplasmSegmentation * cytoSeg = new ftk::CytoplasmSegmentation();
	cytoSeg->SetDataInput(myImg, "data_channel", cytChannel);
	cytoSeg->SetNucleiInput(labImg, "label_image");
	cytoSeg->Run();
	cytoSeg->SaveOutputImage();
	delete cytoSeg;

	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(myImg,labImg,cytChannel,1);
	iCalc->SetFeaturePrefix("cyto_");
	iCalc->Append(table);
	delete iCalc;

	this->updateViews();
}

void NucleusEditor::segmentImage()
{
	QVector<QString> chs;
	int numChannels = myImg->GetImageInfo()->numChannels;
	for (int i=0; i<numChannels; ++i)
	{
		chs.push_back(QString::fromStdString(myImg->GetImageInfo()->channelNames.at(i)));
	}

	//Get the paramFile to use:
	QString paramFile = "";
	ParamsFileDialog *dialog = new ParamsFileDialog(lastPath,chs,this);
	if( dialog->exec() )	
	{
		paramFile = dialog->getFileName();
		nucChannel = dialog->getChannelNumber();
	}
	delete dialog;

	if(nucSeg) delete nucSeg;
	nucSeg = new ftk::NuclearSegmentation();
	nucSeg->SetInput(this->myImg, "nuc_channel", nucChannel);
	nucSeg->SetParameters(paramFile.toStdString());

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
	if(nucsegThread)
	{
		nucsegThread->wait();
		delete nucsegThread;
		nucsegThread = NULL;
	}
	if(featuresThread)
	{
		featuresThread->wait();
		delete featuresThread;
		featuresThread = NULL;
	}
	segmentTool->setVisible(false);
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
	segmentProgress->setValue(segmentState);
	segmentState++;
	this->segment();
}

void NucleusEditor::segment()
{
	switch(segmentState)
	{
	case 0:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		menusEnabled(false);
		segmentTaskLabel->setText(tr(" Binarizing "));
		segmentProgress->setValue(segmentState);
		segmentState++;
		nucsegThread = new NucSegThread(nucSeg);
		connect(nucsegThread, SIGNAL(finished()), this, SLOT(segment()));
		nucsegThread->start();
		break;
	case 1:
		segmentProgress->setValue(segmentState);
		segmentTaskLabel->setText(tr(" Seeds "));
		segmentState++;
		nucsegThread->start();
		break;
	case 2:
		segmentProgress->setValue(segmentState);
		segmentTaskLabel->setText(tr(" Clustering "));
		segmentState++;
		nucsegThread->start();
		break;
	case 3:	//Final State we get to after successful clustering - for checking parameters
		segmentProgress->setValue(segmentState);
		segView->SetLabelImage(nucSeg->GetLabelImage());
		segmentTaskLabel->setText(tr(" Inspect "));
		QApplication::restoreOverrideCursor();
		segmentStop->setEnabled(true);
		segmentContinue->setEnabled(true);
		segmentState++;
		//I should pop-up a msg here:
		break;
	case 4:
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		segmentStop->setEnabled(false);
		segmentContinue->setEnabled(false);
		segmentTaskLabel->setText(tr(" Finalizing "));
		segmentProgress->setValue(segmentState);
		segmentState++;
		nucsegThread->start();
		break;
	case 5:
		if(nucsegThread)
		{
			delete nucsegThread;
			nucsegThread = NULL;
		}
		segmentProgress->setValue(segmentState);
		if(selection) delete selection;
		selection = new ObjectSelection();
		labImg = nucSeg->GetLabelImage();
		segView->SetLabelImage(labImg,selection);
		segmentTaskLabel->setText(tr(" Features "));

		featuresThread = new Features(myImg, labImg, nucChannel, "");
		connect(featuresThread, SIGNAL(finished()), this, SLOT(segment()));
		segmentState++;
		featuresThread->start();
		break;
	case 6:
		if(featuresThread)
		{
			table = featuresThread->GetTable();
			delete featuresThread;
			featuresThread = NULL;
		}
		segmentProgress->setValue(segmentState);
		segmentTaskLabel->setText(tr(" DONE "));
		
		CreateNewTableWindow();
		CreateNewPlotWindow();

		QApplication::restoreOverrideCursor();
		menusEnabled(true);
		cytoAction->setEnabled(true);
		setMenusForResult(true);
		showBoundsAction->setChecked(true);
		showIDsAction->setChecked(true);

		segmentState = -1;

		//Now remove the toolbar:
		segmentTool->setVisible(false);
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
NucSegThread::NucSegThread(ftk::NuclearSegmentation *seg)
: QThread()
{
	mySeg = seg;
	step = 0;
}

void NucSegThread::run()
{
	switch(step)
	{
	case 0:
		mySeg->Binarize(false);
		step++;
		break;
	case 1:
		mySeg->DetectSeeds(false);
		step++;
		break;
	case 2:
		mySeg->RunClustering(true);
		step++;
		break;
	case 3:
		mySeg->Finalize();
		mySeg->ReleaseSegMemory();
		step++;
		break;
	}
}

Features::Features(ftk::Image::Pointer dImg, ftk::Image::Pointer lImg, int chan, std::string pre)
: QThread()
{
	dataImg = dImg;
	labImg = lImg;
	channel = chan;
	prefix = pre;
	table = NULL;
}

void Features::run()
{
	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(dataImg,labImg,channel,0);
	iCalc->SetFeaturePrefix(prefix);
	table = iCalc->Compute();
	delete iCalc;
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************











