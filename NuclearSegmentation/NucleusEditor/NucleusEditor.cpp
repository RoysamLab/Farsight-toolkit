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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

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
	chSignalMapper = NULL;

	segView = new LabelImageViewQT(&colorItemsMap);
	connect(segView, SIGNAL(mouseAt(int,int,int, int,list<int>)), this, SLOT(setMouseStatus(int,int,int, int, list<int>)));
	connect(segView, SIGNAL(autoMerge()), this, SLOT(mergeCells()));
	selection = new ObjectSelection();
	this->setCentralWidget(segView);

	createMenus();
	createStatusBar();
	createProcessToolBar();

	setEditsEnabled(false);
	setCommonEnabled(true);
	modelsMenu->setEnabled(false);

	setWindowTitle(tr("FARSIGHT: Nucleus Editing Tool"));

	standardImageTypes = tr("Images (*.tif *.tiff *.pic *.png *.jpg *.lsm)\n"
							"XML Image Definition (*.xml)\n"
							"All Files (*.*)");

	tblWin.clear();
	pltWin.clear();
//	hisWin.clear();
	pWizard=NULL;

	myImg = NULL;
	labImg = NULL;

	nucSeg = NULL;
	pProc = NULL;
	processThread = NULL;
	table = NULL;
	NucAdjTable = NULL;
	CellAdjTable = NULL;
	
	kplsRun = 0;	//This flag gets set after kpls has run to make sure we show the colored centroids!!
	trainName = 0;
	predictName = 0;

	this->resize(800,800);

	this->readSettings();
  raise();
  activateWindow();
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );
}

NucleusEditor::~NucleusEditor()
{
	if(selection) delete selection;
	if(nucSeg) delete nucSeg;
	if(pProc) delete pProc;

}

//******************************************************************************
//Write/read settings functions
//******************************************************************************
void NucleusEditor::readSettings()
{
	QSettings settings;
	lastPath = settings.value("lastPath", ".").toString();

	colorItemsMap["Selected Objects"] = settings.value("colorForSelections", Qt::yellow).value<QColor>();
	colorItemsMap["Object Boundaries"] = settings.value("colorForBounds", Qt::cyan).value<QColor>();
	colorItemsMap["Object IDs"] = settings.value("colorForIDs", Qt::green).value<QColor>();
	colorItemsMap["ROI Boundary"] = settings.value("colorForROI", Qt::gray).value<QColor>();

}

void NucleusEditor::writeSettings()
{
	QSettings settings;
	settings.setValue("lastPath", lastPath);

	settings.setValue("colorForSelections", colorItemsMap["Selected Objects"]);
	settings.setValue("colorForBounds", colorItemsMap["Object Boundaries"]);
	settings.setValue("colorForIDs", colorItemsMap["Object IDs"]);
	settings.setValue("colorForROI", colorItemsMap["ROI Boundary"]);
}

//******************************************************************************

//******************************************************************************
// Here we just show a message in the status bar when loading
//******************************************************************************
void NucleusEditor::createStatusBar()
{
    statusLabel = new QLabel(tr(" Ready"));
    statusBar()->addWidget(statusLabel, 1);
}

//** ToolBar Setup
void NucleusEditor::createProcessToolBar()
{
	processToolbar = this->addToolBar(tr("Process"));
	processToolbar->setVisible(false);

	processAbort = new QAction(QIcon(":/icons/abort.png"), tr("ABORT"), this);
	processAbort->setToolTip(tr("Abort Processing"));
	processAbort->setEnabled(true);
	connect(processAbort, SIGNAL(triggered()), this, SLOT(abortProcess()));
	processToolbar->addAction(processAbort);

	processContinue = new QAction(QIcon(":/icons/go.png"), tr("GO"), this);
	processContinue->setToolTip(tr("Continue Processing"));
	processContinue->setEnabled(false);
	connect(processContinue, SIGNAL(triggered()),this, SLOT(continueProcess()));
	processToolbar->addAction(processContinue);

    processProgress = new QProgressBar();
	processProgress->setRange(0,0);
	processProgress->setTextVisible(false);
	processToolbar->addWidget(processProgress);

	processTaskLabel = new QLabel;
	processTaskLabel->setFixedWidth(250);
	processToolbar->addWidget(processTaskLabel);
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

	loadProjectAction = new QAction(tr("Load Project..."), this);
	loadProjectAction->setStatusTip(tr("Load ftk project"));
	loadProjectAction->setShortcut(tr("Ctrl+L"));
	connect(loadProjectAction, SIGNAL(triggered()), this, SLOT(loadProject()));
	fileMenu->addAction(loadProjectAction);

	loadImageAction = new QAction(tr("Load Image..."), this);
	loadImageAction->setStatusTip(tr("Load an image into the 5D image browser"));
	loadImageAction->setShortcut(tr("Ctrl+O"));
	connect(loadImageAction, SIGNAL(triggered()), this, SLOT(askLoadImage()));
	fileMenu->addAction(loadImageAction);

	loadLabelAction = new QAction(tr("Load Result..."), this);
	loadLabelAction->setStatusTip(tr("Load a result image into the image browser"));
	connect(loadLabelAction,SIGNAL(triggered()), this, SLOT(askLoadResult()));
	fileMenu->addAction(loadLabelAction);

	loadTableAction = new QAction(tr("Load Table..."), this);
	loadTableAction->setStatusTip(tr("Load data table from text file"));
	connect(loadTableAction, SIGNAL(triggered()), this, SLOT(askLoadTable()));
	fileMenu->addAction(loadTableAction);

	fileMenu->addSeparator();

	processProjectAction = new QAction(tr("Process Image..."), this);
	processProjectAction->setStatusTip(tr("Choose project definition file to use in processing the current image"));
	connect(processProjectAction, SIGNAL(triggered()), this, SLOT(processProject()));
	fileMenu->addAction(processProjectAction);

	fileMenu->addSeparator();

	saveProjectAction = new QAction(tr("Save Project.."), this);
	saveProjectAction->setStatusTip(tr("Save the active project files..."));
	saveProjectAction->setShortcut(tr("Ctrl+S"));
	connect(saveProjectAction, SIGNAL(triggered()), this, SLOT(saveProject()));
	fileMenu->addAction(saveProjectAction);

	saveImageAction = new QAction(tr("Save Image..."), this);
	saveImageAction->setStatusTip(tr("Save image currently displayed"));
	connect(saveImageAction, SIGNAL(triggered()), this, SLOT(askSaveImage()));
	fileMenu->addAction(saveImageAction);

	saveLabelAction = new QAction(tr("Save Result..."), this);
	saveLabelAction->setStatusTip(tr("Save a segmentation result image"));
	connect(saveLabelAction, SIGNAL(triggered()), this, SLOT(askSaveResult()));
	fileMenu->addAction(saveLabelAction);

	saveTableAction = new QAction(tr("Save Table..."), this);
	saveTableAction->setStatusTip(tr("Save the features table"));
	connect(saveTableAction, SIGNAL(triggered()), this, SLOT(askSaveTable()));
	fileMenu->addAction(saveTableAction);

	saveDisplayAction = new QAction(tr("Save Display Image..."), this);
	saveDisplayAction->setStatusTip(tr("Save displayed image to file"));
	connect(saveDisplayAction, SIGNAL(triggered()), this, SLOT(saveDisplayImageToFile()));
	fileMenu->addAction(saveDisplayAction);

	saveCompositeAction = new QAction(tr("Save Composite Image..."), this);
	saveCompositeAction->setStatusTip(tr("Save displayed composite image to file"));
	connect(saveCompositeAction, SIGNAL(triggered()), this, SLOT(saveCompositeImageToFile()));
	fileMenu->addAction(saveCompositeAction);

	fileMenu->addSeparator();

    exitAction = new QAction(tr("Exit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
    fileMenu->addAction(exitAction);

	//VIEW MENU
	viewMenu = menuBar()->addMenu(tr("&View"));

	setPreferencesAction = new QAction(tr("Preferences..."), this);
	connect(setPreferencesAction, SIGNAL(triggered()), this, SLOT(setPreferences()));
	viewMenu->addAction(setPreferencesAction);

	viewMenu->addSeparator();

	showCrosshairsAction = new QAction(tr("Show Selection Crosshairs"), this);
	showCrosshairsAction->setCheckable(true);
	showCrosshairsAction->setChecked( segView->GetCrosshairsVisible() );
	showCrosshairsAction->setStatusTip(tr("Show Crosshairs at selected object"));
	connect(showCrosshairsAction, SIGNAL(triggered()), this, SLOT(toggleCrosshairs()));
	viewMenu->addAction(showCrosshairsAction);

	showBoundsAction = new QAction(tr("Show &Boundaries"), this);
	showBoundsAction->setCheckable(true);
	showBoundsAction->setChecked( segView->GetBoundsVisible() );
	showBoundsAction->setStatusTip(tr("Draw boundaries using a label image"));
	showBoundsAction->setShortcut(tr("Ctrl+B"));
	connect(showBoundsAction, SIGNAL(triggered()), this, SLOT(toggleBounds()));
	viewMenu->addAction(showBoundsAction);

	showIDsAction = new QAction(tr("Show Object IDs"), this);
	showIDsAction->setCheckable(true);
	showIDsAction->setChecked( segView->GetIDsVisible() );
	showIDsAction->setStatusTip(tr("Show IDs of objects"));
	showIDsAction->setShortcut(tr("Ctrl+I"));
	connect(showIDsAction, SIGNAL(triggered()), this, SLOT(toggleIDs()));
	viewMenu->addAction(showIDsAction);

	showCentroidsAction = new QAction(tr("Show Object Centroids"), this);
	showCentroidsAction->setCheckable(true);
	showCentroidsAction->setChecked( segView->GetCentroidsVisible() );
	showCentroidsAction->setStatusTip(tr("Show Centroids of objects"));
	//showCentroidsAction->setShortcut(tr(""));
	connect(showCentroidsAction, SIGNAL(triggered()), this, SLOT(toggleCentroids()));
	viewMenu->addAction(showCentroidsAction);

	adjacencyMenu = viewMenu->addMenu(tr("Show Adjacency"));

	showNucAdjAction = new QAction(tr("Nuclear"), this);
	showNucAdjAction->setCheckable(true);
	showNucAdjAction->setChecked( segView->GetNucAdjVisible() );
	showNucAdjAction->setStatusTip(tr("Show adjacency of nuclei"));
	showNucAdjAction->setShortcut(tr("Shift+N"));
	connect(showNucAdjAction, SIGNAL(triggered()), this, SLOT(toggleNucAdjacency()));
	adjacencyMenu->addAction(showNucAdjAction);

	showCellAdjAction = new QAction(tr("Cellular"), this);
	showCellAdjAction->setCheckable(true);
	showCellAdjAction->setChecked( segView->GetCellAdjVisible() );
	showCellAdjAction->setStatusTip(tr("Show adjacency of cells"));
	showCellAdjAction->setShortcut(tr("Shift+C"));
	connect(showCellAdjAction, SIGNAL(triggered()), this, SLOT(toggleCellAdjacency()));
	adjacencyMenu->addAction(showCellAdjAction);

	zoomMenu = viewMenu->addMenu(tr("Zoom"));

	zoomInAction = new QAction(tr("Zoom In"), this);
	zoomInAction->setStatusTip(tr("Zoom In On The Displayed Image"));
	zoomInAction->setShortcut(tr("="));
	connect(zoomInAction, SIGNAL(triggered()), segView, SLOT(zoomIn()));
	zoomMenu->addAction(zoomInAction);

	zoomOutAction  = new QAction(tr("Zoom Out"), this);
	zoomOutAction->setStatusTip(tr("Zoom Out of The Displayed Image"));
	zoomOutAction->setShortcut(tr("-"));
	connect(zoomOutAction, SIGNAL(triggered()), segView, SLOT(zoomOut()));
	zoomMenu->addAction(zoomOutAction);

	displayChannelMenu = viewMenu->addMenu(tr("Display Channel"));
	connect(displayChannelMenu, SIGNAL(aboutToShow()), this, SLOT(DisplayChannelsMenu()));

	viewMenu->addSeparator();

	newTableAction = new QAction(tr("New Table"), this);
	newTableAction->setStatusTip(tr("Open a new table window"));
	connect(newTableAction, SIGNAL(triggered()), this, SLOT(CreateNewTableWindow()));
	viewMenu->addAction(newTableAction);

	newScatterAction = new QAction(tr("New Scatterplot"), this);
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

	newHistoAction = new QAction(tr("New Histogram"),this);
	newHistoAction->setStatusTip(tr("Open a new Histogram Window"));
	connect(newHistoAction,SIGNAL(triggered()),this,SLOT(CreateNewHistoWindow()));
	viewMenu->addAction(newHistoAction);

	ragMenu = viewMenu->addMenu(tr("New Region Adjacency Graph"));

	nucRagAction = new QAction(tr("Nuclear Adjacency Graph"), this);	
	connect(nucRagAction, SIGNAL(triggered()), this, SLOT(CreateNewNucRAG()));
	ragMenu->addAction(nucRagAction);

	cellRagAction = new QAction(tr("Cellular Adjacency Graph"), this);
	connect(cellRagAction, SIGNAL(triggered()), this, SLOT(CreateNewCellRAG()));
	ragMenu->addAction(cellRagAction);

	imageIntensityAction = new QAction(tr("Adjust Image Intensity"), this);
	imageIntensityAction->setStatusTip(tr("Allows modification of image intensity"));
	connect(imageIntensityAction, SIGNAL(triggered()), segView, SLOT(AdjustImageIntensity()));
	viewMenu->addAction(imageIntensityAction);

	//TOOL MENU
	toolMenu = menuBar()->addMenu(tr("Tools"));

	roiMenu = toolMenu->addMenu(tr("Region Of Interest"));

	drawROIAction = new QAction(tr("Draw ROI"), this);
	connect(drawROIAction, SIGNAL(triggered()), this, SLOT(startROI()));
	roiMenu->addAction(drawROIAction);
	
	clearROIAction = new QAction(tr("Clear ROI"), this);
	connect(clearROIAction, SIGNAL(triggered()), this, SLOT(clearROI()));
	roiMenu->addAction(clearROIAction);

	saveROIAction = new QAction(tr("Save ROI Mask..."), this);
	connect(saveROIAction, SIGNAL(triggered()), this, SLOT(saveROI()));
	roiMenu->addAction(saveROIAction);

	loadROIAction = new QAction(tr("Load ROI Mask..."), this);
	connect(loadROIAction, SIGNAL(triggered()), this, SLOT(loadROI()));
	roiMenu->addAction(loadROIAction);

	roiStatsAction = new QAction(tr("Compute ROI Statistics"), this);
	connect(roiStatsAction, SIGNAL(triggered()), this, SLOT(roiStatistics()));
	toolMenu->addAction(roiStatsAction);

	preprocessAction = new QAction(tr("Preprocess Image..."), this);
	connect(preprocessAction, SIGNAL(triggered()), this, SLOT(preprocessImage()));
	toolMenu->addAction(preprocessAction);

	segmentNucleiAction = new QAction(tr("Segment Nuclei..."), this);
	connect(segmentNucleiAction, SIGNAL(triggered()), this, SLOT(segmentNuclei()));
	toolMenu->addAction(segmentNucleiAction);

	editNucleiAction = new QAction(tr("Edit Nuclei"), this);
	connect(editNucleiAction, SIGNAL(triggered()), this, SLOT(startEditing()));
	toolMenu->addAction(editNucleiAction);

	svmAction = new QAction(tr("Detect Outliers"), this);
	connect(svmAction, SIGNAL(triggered()), this, SLOT(startSVM()));
	toolMenu->addAction(svmAction);

	databaseAction = new QAction(tr("Update Database"), this);
	connect(databaseAction, SIGNAL(triggered()), this, SLOT(updateDatabase()));
	toolMenu->addAction(databaseAction);

	classifyMenu = toolMenu->addMenu(tr("Classifier"));

	trainAction = new QAction(tr("Train"), this);
	connect(trainAction, SIGNAL(triggered()), this, SLOT(startTraining()));
	classifyMenu->addAction(trainAction);

	kplsAction = new QAction(tr("Classify"), this);
	connect(kplsAction, SIGNAL(triggered()), this, SLOT(startKPLS()));
	classifyMenu->addAction(kplsAction);

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
	//editMenu->addAction(visitAction);

	editMenu->addSeparator();

	classAction = new QAction(tr("Change Class"), this);
	classAction->setStatusTip(tr("Modify the class designation for the selected objects"));
	connect(classAction, SIGNAL(triggered()), this, SLOT(changeClass()));
	editMenu->addAction(classAction);

	addAction = new QAction(tr("Add Object"), this);
	addAction->setStatusTip(tr("Draw a Box to add a new object"));
	addAction->setShortcut(tr("Ctrl+A"));
	connect(addAction, SIGNAL(triggered()), segView, SLOT(GetBox()));
	connect(segView, SIGNAL(boxDrawn(int,int,int,int,int)), this, SLOT(addCell(int,int,int,int,int)));
	editMenu->addAction(addAction);

	mergeAction = new QAction(tr("Merge Objects"), this);
	mergeAction->setStatusTip(tr("Merge Objects"));
	mergeAction->setShortcut(tr("Ctrl+M"));
	connect(mergeAction, SIGNAL(triggered()), this, SLOT(mergeCells()));
	editMenu->addAction(mergeAction);

	deleteAction = new QAction(tr("Delete Objects"), this);
	deleteAction->setStatusTip(tr("Deletes the selected objects"));
	deleteAction->setShortcut(tr("Ctrl+D"));
	connect(deleteAction,SIGNAL(triggered()),this,SLOT(deleteCells()));
	editMenu->addAction(deleteAction);

	fillAction = new QAction(tr("Fill Objects"), this);
	fillAction->setStatusTip(tr("Fill holes in the selected objects"));
	fillAction->setShortcut(tr("Ctrl+F"));
	connect(fillAction,SIGNAL(triggered()),this,SLOT(fillCells()));	
	editMenu->addAction(fillAction);

	splitZAction = new QAction(tr("Split Objects At Z"), this);
	splitZAction->setStatusTip(tr("Split selected objects along the current Z slice"));
	splitZAction->setShortcut(tr("Ctrl+T"));
	connect(splitZAction, SIGNAL(triggered()), this, SLOT(splitCellAlongZ()));
	editMenu->addAction(splitZAction);

	splitAction = new QAction(tr("Split Objects X-Y"), this);
	splitAction->setStatusTip(tr("Split objects by choosing two seed points"));
	splitAction->setShortcut(tr("Ctrl+P"));
	splitAction->setCheckable(true);
	splitAction->setChecked(false);
	connect(splitAction, SIGNAL(triggered()), this, SLOT(splitCells()));
	connect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));
	editMenu->addAction(splitAction);

	editMenu->addSeparator();

	exclusionAction = new QAction(tr("Apply Exclusion Margin..."), this);
	exclusionAction->setStatusTip(tr("Set parameters for exclusion margin"));
	connect(exclusionAction, SIGNAL(triggered()), this, SLOT(applyExclusionMargin()));
	editMenu->addAction(exclusionAction);

	// PRE-PROCESSING:
	createPreprocessingMenu();

	// MODELS MENU
	modelsMenu = menuBar()->addMenu(tr("&Models"));
    
    createTrainingAction = new QAction(tr("Create Training Model..."), this);
    connect(createTrainingAction, SIGNAL(triggered()), this, SLOT(createTrainer()));
    modelsMenu->addAction(createTrainingAction);
    
    appendTrainingAction = new QAction(tr("Append Training Model..."), this);
    connect(appendTrainingAction, SIGNAL(triggered()), this, SLOT(appendTrainer()));
    modelsMenu->addAction(appendTrainingAction);

	//QUERIES MENU
	queriesMenu = menuBar()->addMenu(tr("Queries"));

	kNearestNeighborsAction = new QAction(tr("Query K Nearest Neighbors..."), this);
    connect(kNearestNeighborsAction, SIGNAL(triggered()), this, SLOT(queryKNearest()));
    queriesMenu->addAction(kNearestNeighborsAction);
    
    inRadiusNeighborsAction = new QAction(tr("Query Neighbors Within Radius..."), this);
    connect(inRadiusNeighborsAction, SIGNAL(triggered()), this, SLOT(queryInRadius()));
    queriesMenu->addAction(inRadiusNeighborsAction);

	queryViewsOffAction = new QAction(tr("Set Query Views Off"), this);
	queryViewsOffAction->setShortcut(tr("Shift+O"));
    connect(queryViewsOffAction, SIGNAL(triggered()), this, SLOT(queryViewsOff()));
    queriesMenu->addAction(queryViewsOffAction);

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);
}

void NucleusEditor::createPreprocessingMenu()
{
	//************************************************************************************************
	//************************************************************************************************
	// PRE-PROCESSING
	//************************************************************************************************
    PreprocessMenu = menuBar()->addMenu(tr("&Preprocessing"));

	cropAction = new QAction(tr("Crop Image"), this);
	cropAction->setStatusTip(tr("Draw and crop the image to rectangle"));
	connect(cropAction, SIGNAL(triggered()), this, SLOT(CropToRegion()));
	//PreprocessMenu->addAction(cropAction);

	blankAction = new QAction(tr("Mask Image"), this);
	blankAction->setStatusTip(tr("Draw a polygon region, and mask all pixels outside this region"));
	connect(blankAction, SIGNAL(triggered()), this, SLOT(startROI()));
	//PreprocessMenu->addAction(blankAction);

	//PreprocessMenu->addSeparator();

	invertPixAction = new QAction(tr("Invert Intensity Filter"), this);
	invertPixAction->setStatusTip(tr("Invert intensity value of each pixel"));
	connect(invertPixAction, SIGNAL(triggered()), this, SLOT(InvertIntensities()));
	PreprocessMenu->addAction(invertPixAction);

    MedianAction = new QAction(tr("Median Filter"), this);
   	MedianAction->setStatusTip(tr("Apply Median Filter "));
   	connect(MedianAction, SIGNAL(triggered()), this, SLOT(MedianFilter()));
   	PreprocessMenu->addAction(MedianAction);


   	AnisotropicAction = new QAction(tr("Gradient Anisotropic Diffusion Filter"), this);
	AnisotropicAction->setStatusTip(tr("Apply Gradient Anisotropic Diffusion Filtering "));
	connect(AnisotropicAction, SIGNAL(triggered()), this, SLOT(AnisotropicDiffusion()));
	PreprocessMenu->addAction(AnisotropicAction);

	CurvAnisotropicAction = new QAction(tr("Curvature Anisotropic Diffusion Filter"), this);
	CurvAnisotropicAction->setStatusTip(tr("Apply Curvature Anisotropic Diffusion Filtering "));
	connect(CurvAnisotropicAction, SIGNAL(triggered()), this, SLOT(CurvAnisotropicDiffusion()));
	PreprocessMenu->addAction(CurvAnisotropicAction);

	GSErodeAction = new QAction(tr("Grayscale Erosion Filter"), this);
	GSErodeAction->setStatusTip(tr("Apply Grayscale Erosion"));
	connect(GSErodeAction, SIGNAL(triggered()), this, SLOT(GrayscaleErode()));
	PreprocessMenu->addAction(GSErodeAction);

	GSDilateAction = new QAction(tr("Grayscale Dilation Filter"), this);
	GSDilateAction->setStatusTip(tr("Apply Grayscale Dilation"));
	connect(GSDilateAction, SIGNAL(triggered()), this, SLOT(GrayscaleDilate()));
	PreprocessMenu->addAction(GSDilateAction);

	GSOpenAction = new QAction(tr("Grayscale Open Filter"), this);
	GSOpenAction->setStatusTip(tr("Apply Grayscale Opening"));
	connect(GSOpenAction, SIGNAL(triggered()), this, SLOT(GrayscaleOpen()));
	PreprocessMenu->addAction(GSOpenAction);

	GSCloseAction = new QAction(tr("Grayscale Close Filter"), this);
	GSCloseAction->setStatusTip(tr("Apply Grayscale Closing"));
	connect(GSCloseAction, SIGNAL(triggered()), this, SLOT(GrayscaleClose()));
	PreprocessMenu->addAction(GSCloseAction);


	SigmoidAction = new QAction(tr("Sigmoid Filter"), this);
	SigmoidAction->setStatusTip(tr("Apply Sigmoid Filter "));
	connect(SigmoidAction, SIGNAL(triggered()), this, SLOT(SigmoidFilter()));
	PreprocessMenu->addAction(SigmoidAction);

	//ResampleAction = new QAction(tr("Resample Image Filter"), this);
	//ResampleAction->setStatusTip(tr("Resample the Image"));
	//connect(ResampleAction, SIGNAL(triggered()), this, SLOT(Resample()));
	//PreprocessMenu->addAction(ResampleAction);

	//************************************************************************************************
	//************************************************************************************************
}

void NucleusEditor::setEditsEnabled(bool val)
{
	editMenu->setEnabled(val);
	clearSelectAction->setEnabled(val);
	visitAction->setEnabled(val);
	addAction->setEnabled(val);
	mergeAction->setEnabled(val);
	deleteAction->setEnabled(val);
	fillAction->setEnabled(val);
	splitZAction->setEnabled(val);
	splitAction->setEnabled(val);
	classAction->setEnabled(val);
	exclusionAction->setEnabled(val);
	if( val )
		editNucleiAction->setEnabled(false);
	else
		editNucleiAction->setEnabled(true);
	if( val )
		segmentNucleiAction->setEnabled(false);
	else
		segmentNucleiAction->setEnabled(true);
}

void NucleusEditor::setCommonEnabled(bool val){
	toolMenu->setEnabled(val);
	viewMenu->setEnabled(val);
}

void NucleusEditor::setPreprocessingEnabled(bool val)
{
	PreprocessMenu->setEnabled(val);
	//AnisotropicAction->setEnabled(val);
	//CurvAnisotropicAction->setEnabled(val);
	//SigmoidAction->setEnabled(val);
	MedianAction->setEnabled(val);
	//GSErodeAction->setEnabled(val);
	//GSDilateAction->setEnabled(val);
	//GSOpenAction->setEnabled(val);
	//GSCloseAction->setEnabled(val);
}

void NucleusEditor::menusEnabled(bool val)
{
	fileMenu->setEnabled(val);
	viewMenu->setEnabled(val);
	editMenu->setEnabled(val);
	toolMenu->setEnabled(val);
	PreprocessMenu->setEnabled(val);
}

//****************************************************************************
// SLOT: about()
//   A brief message about Farsight is displayed
//****************************************************************************
void NucleusEditor::about()
{
	QString version = QString::number(CPACK_PACKAGE_VERSION_MAJOR);
	version += ".";
	version += QString::number(CPACK_PACKAGE_VERSION_MINOR);
	version += ".";
	version += QString::number(CPACK_PACKAGE_VERSION_PATCH);

	QString text = tr("<h2>FARSIGHT ") + version + tr("</h2>");
	text += tr("<h3>Rensselear Polytechnic Institute</h3>");
	text += tr("<a><u>http://www.farsight-toolkit.org</a></u>");

	QMessageBox::about(this, tr("About FARSIGHT"), text);
}

//******************************************************************************
// SLOT: changes the status bar to say the mouse coordinates
//******************************************************************************
void NucleusEditor::setMouseStatus(int x, int y, int z, int t, list<int> v)
{
	QString statusMsg("X: " + QString::number(x) + ", Y: " + QString::number(y) + ", Z: " + QString::number(z) + ", T: " + QString::number(t));
  
  int i = 1;
  for(list<int>::iterator it = v.begin(); it != v.end(); it++) {
   statusMsg.append(", Value " + QString::number(i)  + ": " + QString::number(*it));
   ++i;
  }

  (this->statusLabel)->setText(statusMsg);
}

//Pop up a message box that asks if you want to save changes 
bool NucleusEditor::askSaveChanges(QString text)
{
	QMessageBox::StandardButtons buttons = QMessageBox::Yes | QMessageBox::No;
	QMessageBox::StandardButton defaultButton = QMessageBox::Yes;
	QMessageBox::StandardButton returnButton = QMessageBox::question(this, tr("Save Changes"), text, buttons, defaultButton);

	if(returnButton == QMessageBox::Yes)
		return true;
	else
		return false;
}
//******************************************************************************
//Reimplement closeEvent to also close all other windows in the application
//******************************************************************************
void NucleusEditor::closeEvent(QCloseEvent *event)
{
	this->abortProcess();

	if(!projectFiles.inputSaved || !projectFiles.outputSaved || !projectFiles.definitionSaved || !projectFiles.tableSaved || !projectFiles.adjTablesSaved)
	{
		if( askSaveChanges(tr("Save changes to the project?")) )
			this->saveProject();
		else
		{
			if(!projectFiles.inputSaved && askSaveChanges(tr("Save changes to the input image?")) )
				this->askSaveImage();
			if(!projectFiles.outputSaved && askSaveChanges(tr("Save changes to the result/label image?")) )
				this->askSaveResult();
			if(!projectFiles.tableSaved && askSaveChanges(tr("Save changes to the table?")) )
				this->askSaveTable();
			if(!projectFiles.adjTablesSaved && askSaveChanges(tr("Save changes to the adjacency tables?")) )
				this->askSaveAdjTables();
		}
	}

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
	this->writeSettings();
	event->accept();
}

//************************************************************************************************
//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
// Saving SLOTS:
//*********************************************************************************
//*********************************************************************************
bool NucleusEditor::saveProject()
{

	if(projectFiles.path == "")
	{
		projectFiles.path = lastPath.toStdString();
	}

	//Make up defaults if needed:
	if(projectFiles.input.size() != 0)
	{
		QString fname = QString::fromStdString(projectFiles.input);
		QString bname = QFileInfo(fname).baseName();

		if(projectFiles.output == "")
			projectFiles.output = bname.toStdString() + "_label.xml";

		if(projectFiles.definition == "")
			projectFiles.definition = bname.toStdString() + "_def.xml";

		if(projectFiles.table == "")
			projectFiles.table = bname.toStdString() + "_table.txt";

		if(projectFiles.adjTables == "")
			projectFiles.adjTables = bname.toStdString() + "_adjTables.txt";

		if(projectFiles.log == "")
			createDefaultLogName();
	}

	ProjectFilenamesDialog *dialog = new ProjectFilenamesDialog(&projectFiles, this);
	if(dialog->exec() == QDialog::Rejected)
		return false;

	if(projectFiles.name != "")
	{
		projectFiles.Write();
		lastPath = QString::fromStdString(projectFiles.path);
	}
	else
		return false;

	if(projectFiles.input != "" && !projectFiles.inputSaved)
	{
		this->saveImage();
	}
	if(projectFiles.output != "" && !projectFiles.outputSaved)
	{
		this->saveResult();
	}
	//if(projectFiles.log != "" && !projectFiles.logSaved)
	//{
	//}
	if(projectFiles.definition != "" && !projectFiles.definitionSaved)
	{
		if( projectDefinition.Write(projectFiles.GetFullDef()) )
			projectFiles.definitionSaved = true;
	}
	if(projectFiles.table != "" && !projectFiles.tableSaved)
	{
		this->saveTable();
	}
	if(projectFiles.adjTables != "" && !projectFiles.adjTablesSaved)
	{
		this->saveAdjTables();
	}

	return true;
}

bool NucleusEditor::askSaveImage()
{
	if(!myImg)
		return false;

	QString filename;

	if(myImg->GetImageInfo()->numChannels == 1)
		filename = QFileDialog::getSaveFileName(this, tr("Save Image As..."),lastPath, standardImageTypes);
	else
		filename = QFileDialog::getSaveFileName(this, tr("Save Image As..."),lastPath, tr("XML Image Definition(*.xml)"));

	if(filename == "")
	{
		return false;
	}
	else
	{
		lastPath = QFileInfo(filename).absolutePath() + QDir::separator();
		projectFiles.path = lastPath.toStdString();
		projectFiles.input = QFileInfo(filename).fileName().toStdString();
	}

	return this->saveImage();
}

bool NucleusEditor::saveImage()
{
	if(!myImg)
		return false;

	if(projectFiles.input == "")
		return false;

	QString fullname = QString::fromStdString( projectFiles.GetFullInput() );
	QString fullbase = QFileInfo(fullname).completeBaseName();
	QString ext = QFileInfo(fullname).suffix();

	bool ok;
	if(ext == "xml")
		ok = ftk::SaveXMLImage(fullname.toStdString(), myImg);
	else
		ok = myImg->SaveChannelAs(0, fullbase.toStdString(), ext.toStdString());

	projectFiles.inputSaved = ok;
	return ok;
}

bool NucleusEditor::askSaveResult()
{
	if(!labImg)
		return false;

	QString filename;

	if(labImg->GetImageInfo()->numChannels == 1)
		filename = QFileDialog::getSaveFileName(this, tr("Save Result As..."),lastPath, tr("TIFF Image (*.tif)"));
	else
		filename = QFileDialog::getSaveFileName(this, tr("Save Result As..."),lastPath, tr("XML Image Definition(*.xml)"));

	if(filename == "")
	{
		return false;
	}
	else
	{
		lastPath = QFileInfo(filename).absolutePath() + QDir::separator();
		projectFiles.path = lastPath.toStdString();
		projectFiles.output = QFileInfo(filename).fileName().toStdString();
	}

	return this->saveResult();
}

bool NucleusEditor::saveResult()
{
	if(!labImg)
		return false;

	if(projectFiles.output == "")
		return false;

	QString fullname = QString::fromStdString( projectFiles.GetFullOutput() );
	QString fullbase = QFileInfo(fullname).completeBaseName();
	QString ext = QFileInfo(fullname).suffix();

	bool ok;
	if(ext == "xml")
		ok = ftk::SaveXMLImage(fullname.toStdString(), labImg);
	else
		ok = labImg->SaveChannelAs(0, fullbase.toStdString(), ext.toStdString());

	projectFiles.outputSaved = ok;
	return ok;
}

bool NucleusEditor::askSaveTable()
{
	if(!table)
		return false;

	QString filename = QFileDialog::getSaveFileName(this, tr("Save Table As..."),lastPath, tr("TEXT(*.txt)"));
	if(filename == "")
		return false;

	lastPath = QFileInfo(filename).absolutePath() + QDir::separator();
	projectFiles.path = lastPath.toStdString();
	projectFiles.table = QFileInfo(filename).fileName().toStdString();

	return this->saveTable();
}

bool NucleusEditor::saveTable()
{
	if(!table)
		return false;

	if(projectFiles.table == "")
		return false;

	bool ok = ftk::SaveTable( projectFiles.GetFullTable(), table);
	projectFiles.tableSaved = ok;
	return ok;
}

bool NucleusEditor::askSaveAdjTables()
{
	if(!NucAdjTable)
		return false;

	QString filename = QFileDialog::getSaveFileName(this, tr("Save Adj_Tables As..."),lastPath, tr("TEXT(*.txt)"));
	if(filename == "")
		return false;

	lastPath = QFileInfo(filename).absolutePath() + QDir::separator();
	projectFiles.path = lastPath.toStdString();
	projectFiles.adjTables = QFileInfo(filename).fileName().toStdString();

	return this->saveAdjTables();
}

bool NucleusEditor::saveAdjTables()
{
	if(!NucAdjTable)
		return false;

	if(projectFiles.adjTables == "")
		return false;

	vtkSmartPointer<vtkTable> adjTables = vtkSmartPointer<vtkTable>::New();
	adjTables->Initialize();
	for(int c=0; c<(int)NucAdjTable->GetNumberOfColumns(); ++c)
	{
	    vtkSmartPointer<vtkStringArray> column = vtkSmartPointer<vtkStringArray>::New();
		column->SetName( NucAdjTable->GetColumnName(c) );
		adjTables->AddColumn(column);
	}
	for(int row1=0; row1<(int)NucAdjTable->GetNumberOfRows(); ++row1)
		adjTables->InsertNextRow(NucAdjTable->GetRow(row1));
		
	if(CellAdjTable)
	{
		vtkSmartPointer<vtkVariantArray> blankRow = vtkSmartPointer<vtkVariantArray>::New();
		blankRow->InsertNextValue(0);
		blankRow->InsertNextValue(0);
		adjTables->InsertNextRow(blankRow);
		for(int row2=0; row2<(int)CellAdjTable->GetNumberOfRows(); ++row2)
			adjTables->InsertNextRow(CellAdjTable->GetRow(row2));		
	}

	bool ok = ftk::SaveTable( projectFiles.GetFullAdjTables(), adjTables);
	projectFiles.adjTablesSaved = ok;
	return ok;
}

void NucleusEditor::createDefaultLogName(void)
{
	QString filename = QString::fromStdString( projectFiles.GetFullInput() );
	QString name = QFileInfo(filename).baseName() + "_log.txt";
	projectFiles.log = name.toStdString();
}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************
// LOADERS:
//************************************************************************************************
//************************************************************************************************
void NucleusEditor::loadProject()
{
	QString filename = QFileDialog::getOpenFileName(this, "Open project...", lastPath, tr("XML Project File(*.xml)"));
	if(filename == "")
		return;

	QString path = QFileInfo(filename).absolutePath();
	//QString name = QFileInfo(filename).baseName();
	lastPath = path;

	projectFiles.Read(filename.toStdString());

	ProjectFilenamesDialog *dialog = new ProjectFilenamesDialog(&projectFiles, this);
	if(dialog->exec() == QDialog::Rejected)
		return;

	if(projectFiles.input != "")
	{
		this->loadImage( QString::fromStdString(projectFiles.GetFullInput()) );
	}
	if(projectFiles.output != "")
	{
		this->loadResult( QString::fromStdString(projectFiles.GetFullOutput()) );
	}
	if(projectFiles.log == "")	//Not opposite boolean here
	{
		this->createDefaultLogName();		//If log file not given, use the default name.
	}
	if(projectFiles.definition != "")
	{
		projectDefinition.Load( projectFiles.GetFullDef() );
		projectFiles.definitionSaved = true;
	}
	if(projectFiles.table != "")
	{
		this->loadTable( QString::fromStdString(projectFiles.GetFullTable()) );
	}
	if(projectFiles.adjTables != "")
	{
		this->loadAdjTables( QString::fromStdString(projectFiles.GetFullAdjTables()) );
	}

	this->startEditing();
}

void NucleusEditor::askLoadTable()
{
	if(!projectFiles.tableSaved && askSaveChanges(tr("Save changes to the current table?")) )
	{
		this->askSaveTable();
	}

	QString fileName  = QFileDialog::getOpenFileName(this, "Select table file to open", lastPath,
								tr("TXT Files (*.txt)"));

    if(fileName == "")
		return;

	this->loadTable(fileName);
}

void NucleusEditor::loadTable(QString fileName)
{
	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();

	table = ftk::LoadTable(fileName.toStdString());
	if(!table) return;

	projectFiles.path = lastPath.toStdString();
	projectFiles.table = QFileInfo(fileName).fileName().toStdString();
	projectFiles.tableSaved = true;

	selection->clear();

	this->closeViews();
	this->CreateNewTableWindow();

	//Get prediction colums, if any, for center map coloring
	prediction_names.clear();
	prediction_names = ftk::GetColumsWithString( "prediction" , table );

	if( !prediction_names.empty() ){
		kplsRun = 1;
		this->updateViews();
	}
	else
		kplsRun = 0;
}

void NucleusEditor::loadAdjTables(QString fileName)
{
	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();

	vtkSmartPointer<vtkTable> adjTables = vtkSmartPointer<vtkTable>::New();
	adjTables = ftk::LoadTable(fileName.toStdString());
	if(!adjTables) return;

	projectFiles.path = lastPath.toStdString();
	projectFiles.adjTables = QFileInfo(fileName).fileName().toStdString();
	projectFiles.adjTablesSaved = true;

	int row = 0;
	int flag = 0; 
	vtkSmartPointer<vtkTable> NuclearTable = vtkSmartPointer<vtkTable>::New();
	NuclearTable->Initialize();
	for(int c=0; c<(int)adjTables->GetNumberOfColumns(); ++c)
	{
	    vtkSmartPointer<vtkStringArray> column = vtkSmartPointer<vtkStringArray>::New();
		column->SetName( adjTables->GetColumnName(c) );
		NuclearTable->AddColumn(column);
	}
	while((adjTables->GetValue(row,0).ToInt() != 0) && (adjTables->GetValue(row,1).ToInt() != 0))
	{
		NuclearTable->InsertNextRow(adjTables->GetRow(row));
		adjTables->RemoveRow(row);
		if((int)adjTables->GetNumberOfRows() == 0)
		{
			flag = 1;
			break;
		}
	}
	NucAdjTable = NuclearTable;
	segView->SetNucAdjTable(NucAdjTable);
	
	if(flag == 0)
	{
		adjTables->RemoveRow(row);

		vtkSmartPointer<vtkTable> CellularTable = vtkSmartPointer<vtkTable>::New();
		CellularTable->Initialize();
		for(int c=0; c<(int)adjTables->GetNumberOfColumns(); ++c)
		{
		    vtkSmartPointer<vtkStringArray> column = vtkSmartPointer<vtkStringArray>::New();
			column->SetName( adjTables->GetColumnName(c) );
			CellularTable->AddColumn(column);
		}
		while((int)adjTables->GetNumberOfRows() != 0)
		{
			CellularTable->InsertNextRow(adjTables->GetRow(row));
			adjTables->RemoveRow(row);
		}
		CellAdjTable = CellularTable;
		segView->SetCellAdjTable(CellAdjTable);
	}
}

void NucleusEditor::askLoadResult(void)
{
	if(!projectFiles.outputSaved && askSaveChanges(tr("Save changes to the result/label image?")) )
	{
		this->askSaveResult();
	}

	QString fileName  = QFileDialog::getOpenFileName(this, "Open File", lastPath, standardImageTypes);
	if(fileName == "") return;

	this->loadResult(fileName);
}

void NucleusEditor::loadResult(QString fileName)
{
	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();
	QString name = QFileInfo(fileName).fileName();
	QString myExt = QFileInfo(fileName).suffix();
	if(myExt == "xml")
	{
		labImg = ftk::LoadXMLImage(fileName.toStdString());
	}
	else
	{
		labImg = ftk::Image::New();
		if(!labImg->LoadFile(fileName.toStdString()))
			labImg = NULL;
	}
	selection->clear();
	segView->SetLabelImage(labImg, selection);
	this->updateNucSeg();

	projectFiles.path = lastPath.toStdString();
	projectFiles.output = name.toStdString();
	projectFiles.outputSaved = true;
}

void NucleusEditor::askLoadImage()
{
	if(!projectFiles.inputSaved && askSaveChanges(tr("Save changes to the input image?")) )
	{
		this->askSaveImage();
	}

	//Get the filename of the new image
	QString fileName = QFileDialog::getOpenFileName(this, "Open Image", lastPath, standardImageTypes);
	//If no filename do nothing
    if(fileName == "")
		return;

	projectFiles.ClearAll();
	this->loadImage(fileName);
}

void NucleusEditor::loadImage(QString fileName)
{
	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();
	QString name = QFileInfo(fileName).fileName();
	QString myExt = QFileInfo(fileName).suffix();
	if(myExt == "xml")
	{
		myImg = ftk::LoadXMLImage(fileName.toStdString());
	}
	else
	{
		myImg = ftk::Image::New();
		if(!myImg->LoadFile(fileName.toStdString()))
			myImg = NULL;
	}
	
	//this->updateNucSeg();
	if(nucSeg)
	{
		segView->SetCenterMapPointer( 0 );
		segView->SetBoundingBoxMapPointer( 0 );
		delete nucSeg;
		nucSeg = NULL;
	}

	segView->SetChannelImage(myImg);
	projectFiles.path = lastPath.toStdString();
	projectFiles.input = name.toStdString();
	projectFiles.inputSaved = true;
}

void NucleusEditor::saveDisplayImageToFile(void)
{
	if(!myImg && !labImg)
		return;

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Display Image"), lastPath, standardImageTypes);
	if(fileName.size() == 0)
		return;

	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();

	segView->SaveDisplayImageToFile(fileName);
}

void NucleusEditor::saveCompositeImageToFile(void)
{
	if(!myImg)
		return;

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Composite Image"), lastPath, standardImageTypes);
	if(fileName.size() == 0)
		return;

	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();

	segView->SaveCompositeImageToFile(fileName);
}

void NucleusEditor::loadROI(void)
{
	QString fileName  = QFileDialog::getOpenFileName(this, tr("Load ROI Mask Image"), lastPath, standardImageTypes);
	if(fileName == "") return;

	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();
	
	segView->GetROIMaskImage()->load(fileName);
	segView->SetROIVisible(true);

	updateROIinTable();
}

void NucleusEditor::saveROI(void)
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save ROI Mask Image As..."),lastPath, standardImageTypes);
	if(filename == "") return;

	lastPath = QFileInfo(filename).absolutePath() + QDir::separator();

	segView->GetROIMaskImage()->save(filename);
}
//************************************************************************************************
//************************************************************************************************
//************************************************************************************************
//************************************************************************************************
//**********************************************************************
// SLOT: start the nuclear associations tool:
//**********************************************************************
/*
void NucleusEditor::startAssociations()
{
	QString fileName = QFileDialog::getOpenFileName(
                             this, "Select file to open", lastPath,
                             tr("XML Association Definition (*.xml)\n"
							    "All Files (*.*)"));
    if(fileName == "")
		return;

	lastPath = QFileInfo(fileName).absolutePath();

	ftk::AssociativeFeatureCalculator * assocCal = new ftk::AssociativeFeatureCalculator();
	assocCal->SetInputFile(fileName.toStdString());

	if(!table)
	{
		table = assocCal->Compute();
		CreateNewTableWindow();
	}
	else
		assocCal->Append(table);

	this->updateViews();

	delete assocCal;
}
*/

void NucleusEditor::startROI(void)
{
	segView->GetROI();
	connect(segView, SIGNAL(roiDrawn()), this, SLOT(endROI()));
}

void NucleusEditor::endROI()
{
	segView->ClearGets();
	disconnect(segView, SIGNAL(roiDrawn()), this, SLOT(endROI()));

	updateROIinTable();
}

void NucleusEditor::updateROIinTable()
{
	if(!table) return;
	if(!nucSeg) return;

	std::map<int, ftk::Object::Point> * cmap = nucSeg->GetCenterMapPointer();
	if(!cmap) return;

	QImage * t_img = segView->GetROIMaskImage();

	const char * columnForROI = "roi";

	//If need to create a new column do so now:
	vtkAbstractArray * output = table->GetColumnByName(columnForROI);
	if(output == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForROI );
		column->SetNumberOfValues( table->GetNumberOfRows() );
		table->AddColumn(column);
	}

	for(int row = 0; (int)row < table->GetNumberOfRows(); ++row)  
	{
		int id = table->GetValue(row,0).ToInt();
		ftk::Object::Point center = (*cmap)[id];
		int val = t_img->pixelIndex( center.x, center.y );
		table->SetValueByName( row, columnForROI, vtkVariant(val) );
	}
	projectFiles.tableSaved = false;
	updateViews();
}

void NucleusEditor::clearROI(void)
{
	const char * columnForROI = "roi";
	
	segView->GetROIMaskImage()->fill(Qt::white);
	segView->SetROIVisible(false);

	if(table)
	{
		table->RemoveColumnByName(columnForROI);
	}
	projectFiles.tableSaved = false;
	updateViews();
}

void NucleusEditor::roiStatistics(void)
{
	QImage * t_img = segView->GetROIMaskImage();
	if(t_img == NULL)
		return;

	typedef itk::Image<unsigned char, 3> ImageType;
	typedef itk::LabelStatisticsImageFilter< ImageType , ImageType > LabelStatisticsType;

	//Change QImage into ITK IMAGE:
	ImageType::Pointer roiImg = ImageType::New();
	ImageType::PointType origin;
   	origin[0] = 0; 
	origin[1] = 0;
	origin[2] = 0;
    roiImg->SetOrigin( origin );
	ImageType::IndexType start = {{ 0,0,0 }};    
	ImageType::SizeType size = {{ t_img->width(), t_img->height(), myImg->GetImageInfo()->numZSlices }};
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    roiImg->SetRegions( region );
	roiImg->Allocate();
	roiImg->FillBuffer(0);
	for(int i=0; i<t_img->width(); ++i)
	{
		for(int j=0; j<t_img->height(); ++j)
		{
			ImageType::IndexType ind = {{ i,j, segView->GetCurrentZ() }};
			roiImg->SetPixel(ind, t_img->pixelIndex( i, j ) );
		}
	}

	//iterate through each channel and compute the statistics
	for(int c=0; c < myImg->GetImageInfo()->numChannels; ++c)
	{
		ImageType::Pointer chImg = myImg->GetItkPtr<unsigned char>(0, c);

		LabelStatisticsType::Pointer labelStatisticsFilter = LabelStatisticsType::New();	//ITK label statistics filter
		labelStatisticsFilter->SetInput( chImg );
		labelStatisticsFilter->SetLabelInput( roiImg );
		labelStatisticsFilter->UseHistogramsOff();

		std::cout << "Calculating Measures for channel # " << c << std::endl;
		try
		{
			labelStatisticsFilter->Update();
		}
		catch (itk::ExceptionObject & e) 
		{
			std::cerr << "Exception in ITK Label Statistics Filter: " << e << std::endl;
			continue;
		}
		
		int label = 1;
		//std::cout << "  sum: " << (float)labelStatisticsFilter->GetSum( label ) << std::endl;
		std::cout << "  min: " << (float)labelStatisticsFilter->GetMinimum( label ) << std::endl;
		std::cout << "  max: " << (float)labelStatisticsFilter->GetMaximum( label ) << std::endl;
		std::cout << "  mean: " << (float)labelStatisticsFilter->GetMean( label ) << std::endl;
		std::cout << "  sigma: " << (float)labelStatisticsFilter->GetSigma( label ) << std::endl;
		std::cout << "  variance: " << (float)labelStatisticsFilter->GetVariance( label ) << std::endl;

	}
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
	pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
	pWizard->show();
}

void NucleusEditor::updateDatabase()
{
	std::vector<std::string> col_names;
	sqlite3 *dbConn;

	if(table){
		dbConn = ftk::sqliteOpenConnection();
		if( dbConn ){
			for (int col = 1; col< table->GetNumberOfColumns(); ++col){
				std::string temp3=table->GetColumnName(col);
				col_names.push_back(temp3);
			}
			std::string table_name = "IMAGE_TEST";
			ftk::checkForUpdate( dbConn, col_names );

			std::string image_name;
			image_name = lastPath.toStdString() + "Nucleus_Editor_Image";
			char *im_nm_cstr = new char [image_name.size()+1];
			strcpy (im_nm_cstr, image_name.c_str());
			char *path_nm_cstr = new char [lastPath.toStdString().size()+1];
			strcpy (path_nm_cstr, lastPath.toStdString().c_str());
			std::vector< double > table_array;
			for (int row = 0; row< table->GetNumberOfRows(); ++row){
				for (int col = 0; col< table->GetNumberOfColumns(); ++col){
					table_array.push_back(table->GetValue(row,col).ToDouble());
				}
			}
			int sql_db_img_id = ftk::GenericInsert( dbConn, im_nm_cstr, table_name.c_str(), path_nm_cstr, table_array,table->GetNumberOfColumns(), table->GetNumberOfRows(), col_names );
			std::cout << "The image ID on the database is: " << sql_db_img_id << std::endl;
		}
	}
}

//**********************************************************************
// SLOT: start the training dialog:
//**********************************************************************
void NucleusEditor::startTraining()
{
	if(!table) return;

	TrainingDialog *d = new TrainingDialog(table, "train", this);
	connect(d, SIGNAL(changedTable()), this, SLOT(updateViews()));
	d->show();
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
	training_names.clear();
	for( int i=0; i<table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = table->GetColumnName(i);
		if( current_column.find("train") != std::string::npos ){
			training_names.push_back( current_column );
			std::string::iterator it;
			it=current_column.begin();
			current_column.erase ( current_column.begin(), current_column.begin()+6 );
			class_names.push_back( current_column );
		}
	}
	if( training_names.empty() ) return;
	trainName = 0;
	if( training_names.size() > 1 ){
		QVector<QString> qtraining_names;
		for (int i=0; i<(int)training_names.size(); ++i){
			qtraining_names << QString::fromStdString(class_names.at(i));
		}
		PredictionDialog *pred_dial = new PredictionDialog(qtraining_names,this);
		pred_dial->show();
		if( pred_dial->exec() ){
			trainName = pred_dial->getTrainNumber();
		}
		delete pred_dial;
	}
	p_name.clear();
	p_name = "prediction_" + class_names.at(trainName);
	std::vector<std::string>::iterator str_it;
	for( str_it = prediction_names.begin(); str_it != prediction_names.end(); ++str_it )
		if( (*str_it).find( p_name.c_str() ) != std::string::npos )
			break;
	if( str_it == prediction_names.end() )
		prediction_names.push_back( p_name );
	pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_KPLS, training_names.at(trainName).c_str(), p_name.c_str(), this);
	connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
	connect(pWizard, SIGNAL(enableModels()), this, SLOT(EnableModels()));
	pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
	pWizard->show();
	kplsRun = 1;
    
}

//**********************************************************************
// SLOT: start the create trainer widget:
//**********************************************************************
void NucleusEditor::createTrainer()
{
	if(!table) return;

	if(pWizard)
	{
		delete pWizard;
	}

	pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_SEGMODEL,"","", this);
	pWizard->setWindowTitle(tr("Select Training Features"));
	pWizard->show();
	
}

//**********************************************************************
// SLOT: start the append trainer widget:
//**********************************************************************
void NucleusEditor::appendTrainer()
{

	//open pattern analysis wizard	
	if(!table) return;
	if(pWizard)
	{
		delete pWizard;
	}

	QString fileName  = QFileDialog::getOpenFileName(this, "Select training model to open", lastPath,
									tr("TXT Files (*.txt)"));
	if(fileName == "")
		return;
	lastPath = QFileInfo(fileName).absolutePath();

	this->loadModelFromFile(fileName.toStdString());

	pWizard = new PatternAnalysisWizard( table, model_table, fileName, PatternAnalysisWizard::_APPENDMODEL,"","", this);
	pWizard->show();
	
}

void NucleusEditor::loadModelFromFile( std::string file_name ){
	model_table = ftk::LoadTable(file_name);
	if(!model_table) return;
	//Append the data to the current table
	//this->GetTrainingNames( model_table );
	//this->Append();
	return;
}

//**********************************************************************
// SLOT: Query K Nearest Neighbors
//**********************************************************************
void NucleusEditor::queryKNearest()
{
	if(!table) return;
	if(!segView) return;
	
	std::vector<unsigned int> IDs;
	unsigned int k;
	unsigned short Class_dest, Class_src = 0;
	
	QVector<QString> classes;
	int max_class = 0;
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		if(table->GetValueByName(row,"prediction_default1").ToInt() > max_class)
			max_class = table->GetValueByName(row,"prediction_default1").ToInt();
	}
	for(int i=0; i<max_class; ++i)
		classes.push_back(QString::number(i+1));
	

	QueryDialog *dialog = new QueryDialog(1,classes,this);
	if( dialog->exec() )
	{
		IDs = dialog->parseIDs();
		k = dialog->parseK();
		if(IDs.at(0) == 0)
			Class_src = dialog->getSourceClass();
		Class_dest = dialog->getDestClass();
    }
	delete dialog;

	std::map<int, ftk::Object::Point> *	centerMap;
	centerMap = segView->GetCenterMapPointer();
	std::map<int, ftk::Object::Point>::iterator it;
	std::map< unsigned int, std::vector<float> > centroidMap;
	for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	{
		unsigned int id = (unsigned int)(*it).first;
		std::vector<float> c;
		c.push_back((float)(*it).second.x);
		c.push_back((float)(*it).second.y);
		c.push_back((float)(*it).second.z);
		centroidMap[id] = c;				
	}

	kNearestObjects* KNObj = new kNearestObjects(centroidMap);
	KNObj->setFeatureTable(table);
	std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	if(IDs.at(0) == 0)
		kNeighborIDs = KNObj->k_nearest_neighbors_All(k, Class_dest, Class_src);
	else
		kNeighborIDs = KNObj->k_nearest_neighbors_IDs(IDs, k, Class_dest);

	std::string Filename = "K_Nearest_Neighbors.txt";
	ofstream outFile; 
	outFile.open(Filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return;
	}
	//Write out the average distance:
	for(int i=0; i < (int)kNeighborIDs.size(); ++i)
	{
		outFile <<  kNeighborIDs.at(i).at(0).first << "\t";
		outFile << average(kNeighborIDs.at(i)) << "\n";
	}
	outFile.close();
	vtkSmartPointer<vtkTable> kNeighborTable = KNObj->vectorsToGraphTable(kNeighborIDs);
	segView->SetKNeighborTable(kNeighborTable);
	segView->SetKNeighborsVisibleOn();		
}

double NucleusEditor::average(std::vector< std::pair<unsigned int, double> > ID)
{
	double dist = 0;
	for(int i=1; i<(int)ID.size(); ++i)
	{
		dist += ID.at(i).second;
	}
	double average = dist/(int)(ID.size()-1);
	return average;
}


//**********************************************************************
// SLOT: Query Neighbors Within Radius
//**********************************************************************
void NucleusEditor::queryInRadius()
{
	if(!table) return;
	if(!segView) return;
	
	std::vector<unsigned int> IDs;
	double radius;
	unsigned short Class_dest, Class_src = 0;

	QVector<QString> classes;
	int max_class = 0;
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		if(table->GetValueByName(row,"prediction_default1").ToInt() > max_class)
			max_class = table->GetValueByName(row,"prediction_default1").ToInt();
	}
	for(int i=0; i<max_class; ++i)
		classes.push_back(QString::number(i+1));

	QueryDialog *dialog = new QueryDialog(2,classes,this);
	if( dialog->exec() )
	{
		IDs = dialog->parseIDs();
		radius = dialog->parseRad();
		if(IDs.at(0) == 0)
			Class_src = dialog->getSourceClass();
		Class_dest = dialog->getDestClass();
    }
	delete dialog;

	std::map<int, ftk::Object::Point> *	centerMap;
	centerMap = segView->GetCenterMapPointer();
	std::map<int, ftk::Object::Point>::iterator it;
	std::map< unsigned int, std::vector<float> > centroidMap;
	for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	{
		unsigned int id = (unsigned int)(*it).first;
		std::vector<float> c;
		c.push_back((float)(*it).second.x);
		c.push_back((float)(*it).second.y);
		c.push_back((float)(*it).second.z);
		centroidMap[id] = c;
	}

	kNearestObjects* KNObj = new kNearestObjects(centroidMap);
	KNObj->setFeatureTable(table);
	std::vector<std::vector< std::pair<unsigned int, double> > > radNeighborIDs;
	if(IDs.at(0) == 0)
		radNeighborIDs = KNObj->neighborsWithinRadius_All(radius, Class_dest, Class_src);
	else
		radNeighborIDs = KNObj->neighborsWithinRadius_IDs(IDs, radius, Class_dest);
	
	std::string Filename = "Neighbors_Within_Radius.txt";
	ofstream outFile; 
	outFile.open(Filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return;
	}
	//Write out the average distance:
	for(int i=0; i < (int)radNeighborIDs.size(); ++i)
	{
		outFile <<  radNeighborIDs.at(i).at(0).first << "\t";
		outFile << average(radNeighborIDs.at(i)) << "\n";
	}
	outFile.close();

	vtkSmartPointer<vtkTable> radNeighborTable = KNObj->vectorsToGraphTable(radNeighborIDs);
	
	segView->SetRadNeighborTable(radNeighborTable);
	segView->SetRadNeighborsVisibleOn();
	
}

//**********************************************************************
// SLOT: Turn Off Query Views
//**********************************************************************
void NucleusEditor::queryViewsOff()
{
	if(!segView) return;
	segView->SetQueryViewsOff();	

}


//*********************************************************************************************************
//Connect the closing signal from views to this slot to remove it from my lists of open views:
//*********************************************************************************************************
void NucleusEditor::viewClosing(QWidget * view)
{
	std::vector<TableWindow *>::iterator table_it;
	for ( table_it = tblWin.begin(); table_it < tblWin.end(); table_it++ )
	{
		if( *table_it == view )
		{
			tblWin.erase(table_it);
			return;
		}
	}

	std::vector<PlotWindow *>::iterator plot_it;
	for ( plot_it = pltWin.begin(); plot_it < pltWin.end(); plot_it++ )
	{
		if( *plot_it == view )
		{
			pltWin.erase(plot_it);
			return;
		}
	}

	std::vector<HistoWindow *>::iterator hist_it;
	for ( hist_it = hisWin.begin(); hist_it < hisWin.end(); hist_it++ )
	{
		if( *hist_it == view )
		{
			hisWin.erase(hist_it);
			return;
		}
	}
}

void NucleusEditor::closeViews()
{
	for(int p=0; p<(int)tblWin.size(); ++p)
		tblWin.at(p)->close();

	for(int p=0; p<(int)pltWin.size(); ++p)
		pltWin.at(p)->close();

	for(int p=0; p<(int)hisWin.size(); ++p)
		hisWin.at(p)->close();
}

//Call this slot when the table has been modified (new rows or columns) to update the views:
void NucleusEditor::updateViews(void)
{
	//Show colored seeds after kPLS has run
	if( kplsRun )
	{
		segView->SetClassMap(table, prediction_names);
		showCentroidsAction->setChecked(true);
		segView->SetCentroidsVisible(true);
		kplsRun = 0;
	}

	segView->update();

	for(int p=0; p<(int)tblWin.size(); ++p)
		tblWin.at(p)->update();

	for(int p=0; p<(int)pltWin.size(); ++p)
		pltWin.at(p)->update();

	for(int p=0; p<(int)hisWin.size(); ++p)
		hisWin.at(p)->update();


}
//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void NucleusEditor::CreateNewPlotWindow(void)
{
	if(!table) return;

	pltWin.push_back(new PlotWindow());
	connect(pltWin.back(), SIGNAL(closing(QWidget *)), this, SLOT(viewClosing(QWidget *)));
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
	connect(tblWin.back(), SIGNAL(closing(QWidget *)), this, SLOT(viewClosing(QWidget *)));
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

	hisWin.push_back(new HistoWindow());
	connect(hisWin.back(), SIGNAL(closing(QWidget *)), this, SLOT(viewClosing(QWidget *)));
	hisWin.back()->setModels(table,selection);
	hisWin.back()->show();
}

//******************************************************************************
// Create Region Adjacency Graphs
//******************************************************************************
void NucleusEditor::CreateNewNucRAG(void)
{
	FTKgraph* NucRAG = new FTKgraph();	
	NucRAG->DisplayGraph(NucAdjTable);
}

void NucleusEditor::CreateNewCellRAG(void)
{
	FTKgraph* CellRAG = new FTKgraph();
	CellRAG->DisplayGraph(CellAdjTable);
}

//******************************************************************************

void NucleusEditor::clearSelections()
{
	if(selection)
	{
		selection->clear();
	}
}

void NucleusEditor::setPreferences()
{
	PreferencesDialog dialog(this->colorItemsMap);
	if(dialog.exec())
	{
		this->colorItemsMap = dialog.GetColorItemsMap();
		segView->update();
	}
}

void NucleusEditor::toggleCrosshairs(void)
{
	if(!segView) return;

	if( showCrosshairsAction->isChecked() )
		segView->SetCrosshairsVisible(true);
	else
		segView->SetCrosshairsVisible(false);
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
	if(!segView) return;

	if( showIDsAction->isChecked() )
		segView->SetIDsVisible(true);
	else
		segView->SetIDsVisible(false);
}

void NucleusEditor::toggleCentroids(void)
{
	if(!segView) return;

	if( prediction_names.size() )
		segView->SetClassMap(table, prediction_names);

	if( showCentroidsAction->isChecked() )
		segView->SetCentroidsVisible(true);
	else
		segView->SetCentroidsVisible(false);
}

void NucleusEditor::toggleNucAdjacency(void)
{
	if(!segView) return;

	if( showNucAdjAction->isChecked() )
		segView->SetNucAdjVisible(true);
	else
		segView->SetNucAdjVisible(false);
}

void NucleusEditor::toggleCellAdjacency(void)
{
	if(!segView) return;

	if( showCellAdjAction->isChecked() )
		segView->SetCellAdjVisible(true);
	else
		segView->SetCellAdjVisible(false);
}

void NucleusEditor::DisplayChannelsMenu()
{
	if( !myImg )
		return;

	std::vector<std::string> channel_names = myImg->GetChannelNames();
	std::vector<bool> channel_status = segView->GetChannelFlags();

	//remove all existing actions;
	for(int i=0; i<(int)displayChannelAction.size(); ++i)
	{
		delete displayChannelAction.at(i);
	}
	displayChannelMenu->clear();
	displayChannelAction.clear();

	if(chSignalMapper)
		delete chSignalMapper;
	chSignalMapper = new QSignalMapper(this);
	for(int i=0; i<(int)channel_names.size(); ++i)
	{
		QAction * action = new QAction( tr(channel_names.at(i).c_str()), this );
		action->setCheckable(true);
		action->setChecked( channel_status.at(i) );
		action->setStatusTip( tr("Turn on/off this channel") );
		action->setShortcut( QString::number(i) );
		connect(action, SIGNAL(triggered()), chSignalMapper, SLOT(map()));
		chSignalMapper->setMapping( action, i );
		displayChannelMenu->addAction(action);
	}
	connect(chSignalMapper, SIGNAL(mapped(int)), this, SLOT(toggleChannel(int)));
}

void NucleusEditor::toggleChannel( int chNum )
{
	std::vector<bool> ch_stats = segView->GetChannelFlags();
	ch_stats[chNum] = !ch_stats[chNum];
	segView->SetChannelFlags( ch_stats );
}
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
// EDITING SLOTS:
//******************************************************************************************
//******************************************************************************************
void NucleusEditor::updateNucSeg(bool ask)
{
	if(nucSeg)
	{
		segView->SetCenterMapPointer( 0 );
		segView->SetBoundingBoxMapPointer( 0 );
		delete nucSeg;
		nucSeg = NULL;
	}

	if(!myImg || !labImg)
		return;

	int nucChannel = projectDefinition.FindInputChannel("NUCLEAR");
	if(nucChannel == -1 && ask)
	{
		nucChannel=requestChannel(myImg);
		projectDefinition.MakeDefaultNucleusSegmentation(nucChannel);
	}
	else
	{
		nucChannel=0;
		if( !projectDefinition.inputs.size() ) //ISAAC: THIS IS A HACK TO GET THE HISTOPATHOGLOGY PROJECT WORKING WILL FIX SOON SORRY -KEDAR
			projectDefinition.MakeDefaultNucleusSegmentation(nucChannel);
	}

	nucSeg = new ftk::NuclearSegmentation();
	nucSeg->SetInput(myImg,"nuc_img", nucChannel);
	nucSeg->SetLabelImage(labImg,"lab_img");
	nucSeg->ComputeAllGeometries();

	
	segView->SetCenterMapPointer( nucSeg->GetCenterMapPointer() );
	segView->SetBoundingBoxMapPointer( nucSeg->GetBoundingBoxMapPointer() );
}

void NucleusEditor::startEditing(void)
{
	std::string log_entry = "NUCLEAR_SEGMENTATION , ";
	log_entry += ftk::NumToString(nucSeg->GetNumberOfObjects()) + " , ";
	log_entry += ftk::TimeStamp();
	ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);

	projectFiles.nucSegValidated = false;
	setEditsEnabled(true);
	setCommonEnabled(true);
}

void NucleusEditor::stopEditing(void)
{
	setEditsEnabled(false);
	setCommonEnabled(true);
	
	if(splitAction->isChecked())
	{
		segView->ClearGets();
		splitAction->setChecked(false);
		disconnect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));
	}

	std::string log_entry = "NUCLEAR_SEGMENTATION_VALIDATED , ";
	log_entry += ftk::TimeStamp();
	ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	projectFiles.nucSegValidated = true;
}

void NucleusEditor::changeClass(void)
{

	//if(!nucSeg) return;

	//std::set<long int> sels = selection->getSelections();
	//std::vector<int> ids(sels.begin(), sels.end());

	//Get the new class number:
	//bool ok;
	//QString msg = tr("Change the class of all selected items to: \n");
	//int newClass = QInputDialog::getInteger(NULL, tr("Change Class"), msg, -1, -1, 10, 1, &ok);

	//Change the class of these objects:
	//if(ok)
	//{
		//nucSeg->SetClass(ids,newClass);
		//this->updateViews();
	//}
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
	int id = nucSeg->AddObject(x1, y1, z, x2, y2, z, table);
	if(id != 0)
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		this->updateViews();
		selection->select(id);

		std::string log_entry = "ADD , ";
		log_entry += ftk::NumToString(id) + " , ";
		log_entry += ftk::TimeStamp();
		ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	}
}

void NucleusEditor::deleteCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());

	if(ids.size() == 0)
		return;

	if(nucSeg->Delete(ids, table))
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		selection->clear();
		this->updateViews();
		for(int j=0; j<(int)ids.size(); ++j)
		{
			int ID = ids.at(j);
			for(int row=0; row<(int)NucAdjTable->GetNumberOfRows(); ++row)
			{
				if((NucAdjTable->GetValue(row,0).ToInt() == ID) || (NucAdjTable->GetValue(row,1).ToInt() == ID))
				{
					NucAdjTable->RemoveRow(row);
					--row;
				}
			}
		}
		segView->SetNucAdjTable(NucAdjTable);

		std::string log_entry = "DELETE , ";
		for(int i=0; i<(int)ids.size(); ++i)
			log_entry += " " + ftk::NumToString(ids.at(i));
		log_entry += " , ";
		log_entry += ftk::TimeStamp();
		ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	}
}

void NucleusEditor::mergeCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	//int newObj = nucSeg->Merge(ids, table);
	std::vector< std::vector<int> > new_grps = nucSeg->GroupMerge(ids, table, NucAdjTable);
	if(new_grps.size() != 0)
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		selection->clear();
		this->updateViews();
		segView->SetNucAdjTable(NucAdjTable);

		for(int i=0; i<(int)new_grps.size(); ++i)
		{
			std::string log_entry = "MERGE , ";
			for(int j=0; j<(int)new_grps.at(i).size()-1; ++j)
			{
				log_entry += " " + ftk::NumToString(new_grps.at(i).at(j));
			}
			log_entry += " , ";
			log_entry += ftk::NumToString(new_grps.at(i).at(new_grps.at(i).size()-1));
			log_entry += " , ";
			log_entry += ftk::TimeStamp();
			ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
		}
	}
}


void NucleusEditor::splitCells(void)
{
	if(splitAction->isChecked())
	{
		segView->Get2Points();
		segView->DoubleClicksOff();
		connect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));
	}
	else
	{
		segView->ClearGets();
		segView->DoubleClicksOn();
		disconnect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));	
	}
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

	std::vector<int> ret = nucSeg->Split(P1, P2, table, NucAdjTable);
	if(ret.size() != 0)
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		selection->clear();
		this->updateViews();

		std::string log_entry = "SPLIT , ";
		log_entry += ftk::NumToString(ret.at(0)) + " , ";
		for(int i=1; i<(int)ret.size(); ++i)
			log_entry += ftk::NumToString(ret.at(i)) + " ";
		log_entry += ", ";
		log_entry += ftk::TimeStamp();
		ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	}

	if(splitAction->isChecked())
	{
		segView->Get2Points();
	}
}

void NucleusEditor::splitCellAlongZ(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	if(sels.size() == 0)
		return;
	if(labImg->GetImageInfo()->numZSlices == 1)
		return;

	selection->clear();
	for ( set<long int>::iterator it=sels.begin(); it != sels.end(); it++ )
	{
		std::vector<int> ret = nucSeg->SplitAlongZ(*it,segView->GetCurrentZ(), table);
		if(ret.size() != 0)
		{
			projectFiles.outputSaved = false;
			projectFiles.tableSaved = false;
			projectFiles.adjTablesSaved = false;

			std::string log_entry = "SPLIT , ";
			log_entry += ftk::NumToString(ret.at(0)) + " , ";
			for(int i=1; i<(int)ret.size(); ++i)
				log_entry += ftk::NumToString(ret.at(i)) + " ";
			log_entry += ", ";
			log_entry += ftk::TimeStamp();
			ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
		}
	}
	this->updateViews();
}

void NucleusEditor::fillCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	if(sels.size() == 0)
		return;
	std::vector<int> ids(sels.begin(), sels.end());

	nucSeg->FillObjects(ids);
	std::string log_entry = "FILL , ";
	for(int i=0; i<(int)ids.size(); ++i)
		log_entry += ftk::NumToString(ids.at(i)) + " ";
	log_entry += ", ";
	log_entry += ftk::TimeStamp();
	ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);

	selection->clear();
	this->updateViews();
}

void NucleusEditor::applyExclusionMargin(void)
{
	//Get the parameters to use for the exclusion margin:
	int l = 0;
	int r = 0;
	int t = 0;
	int b = 0;
	int z1 = 0;
	int z2 = 0;
	ExclusionDialog *dialog = new ExclusionDialog(segView->GetDisplayImage(), this);
	if( !dialog->exec() )
	{
		delete dialog;
		return;
	}
	
	l = dialog->getMargin(0);
	r = dialog->getMargin(1);
	t = dialog->getMargin(2);
	b = dialog->getMargin(3);
	z1 = dialog->getMargin(4);
	z2 = dialog->getMargin(5);

	delete dialog;

	selection->clear();
	if(nucSeg->Exclude(l, r, t, b, z1, z2, table))
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		this->updateViews();

		std::string log_entry = "EXCLUSION_MARGIN, ";
		log_entry += ftk::NumToString(l) + " " + ftk::NumToString(r) + " ";
		log_entry += ftk::NumToString(t) + " " + ftk::NumToString(b) + " ";
		log_entry += ftk::NumToString(z1)+ " " + ftk::NumToString(z2)+ " , ";
		log_entry += ftk::TimeStamp();
		ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	}
}

int NucleusEditor::requestChannel(ftk::Image::Pointer img)
{
	QStringList chs;
	int numChannels = img->GetImageInfo()->numChannels;
	for (int i=0; i<numChannels; ++i)
	{
		chs << QString::number(i) + ". " + QString::fromStdString(img->GetImageInfo()->channelNames.at(i));
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

void NucleusEditor::CreateDefaultAssociationRules()
{
	//filename of the label image to use:
	QString fname = QString::fromStdString(projectFiles.output);
	QFileInfo inf(QDir(lastPath), QFileInfo(fname).baseName() + "_nuc.tif");
	std::string label_name = inf.absoluteFilePath().toStdString();
	//Create association rules:
	//1. Pass the constructor the filename of the label image to use and the number of associations to compute:
	ftk::NuclearAssociationRules * objAssoc = new ftk::NuclearAssociationRules(label_name, 1);
	std::vector<std::string> targFileNames = myImg->GetFilenames();
	std::string targFileName = targFileNames.at(projectDefinition.FindInputChannel("CYTOPLASM"));
	//2. Create each association rule: name, filename of target image, outside distance, inside distance, whole object, type
	objAssoc->AddAssociation("nuc_CK", targFileName, 0, 0, true, false, false, 0, 0, 3, "");
	//filename of the output xml and save:
	QFileInfo inf2(QDir(lastPath), QFileInfo(fname).baseName() + "_assoc.xml");
	//3. Write to file:
	objAssoc->WriteRulesToXML(inf2.absoluteFilePath().toStdString());
	//4. put the output xml in the project definition:
	//typedef struct { std::string name; std::string value; } StringParameter;
	//ftk::ProjectDefinition::StringParameter tpm;	//ISAAC: THIS IS A HACK TO GET THE HISTOPATHOGLOGY PROJECT WORKING WILL FIX SOON SORRY -KEDAR
	//tpm.value=inf2.absoluteFilePath().toStdString();	//ISAAC: THIS IS A HACK TO GET THE HISTOPATHOGLOGY PROJECT WORKING WILL FIX SOON SORRY -KEDAR
	//projectDefinition.associationRules.push_back( tpm );	//ISAAC: THIS IS A HACK TO GET THE HISTOPATHOGLOGY PROJECT WORKING WILL FIX SOON SORRY -KEDAR

	//Now create other rules:

	//AND create other XML association definition files for other label images:

}

//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//*****************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
// CREATE A NUCLEAR SEGMENTATION PROJECT AND THEN GET IT STARTED FOR THE CURRENT IMAGE:
void NucleusEditor::segmentNuclei()
{
	if(!myImg) return;

	//Get Channels in current Image:
	QVector<QString> chs = getChannelStrings();

	//Get the paramFile and channel to use:
	QString paramFile = "";
	int nucChannel = 0;
	ParamsFileDialog *dialog = new ParamsFileDialog(lastPath,chs,this);
	if( dialog->exec() )
	{
		//paramFile = dialog->getFileName();
		nucChannel = dialog->getChannelNumber();
		projectDefinition.MakeDefaultNucleusSegmentation(nucChannel);
		projectFiles.definitionSaved = false;
		projectFiles.nucSegValidated = false;
        
		startProcess();
		;
		
	}
	delete dialog;

}

QVector<QString> NucleusEditor::getChannelStrings(void)
{
	QVector<QString> chs;
	int numChannels = myImg->GetImageInfo()->numChannels;
	for (int i=0; i<numChannels; ++i)
	{
		chs << QString::number(i) + ". " + QString::fromStdString(myImg->GetImageInfo()->channelNames.at(i));
	}
	return chs;
}

//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
//******************************************************************************************
// TAKE CARE OF PROCESSING PROJECTS
//******************************************************************************************
//******************************************************************************************
void NucleusEditor::processProject(void)
{
	QString projectName = QFileDialog::getOpenFileName(
                             this, "Select Definition File", lastPath,
                             tr("XML Project Definition (*.xml)\n"
							    "All Files (*.*)"));
	if(projectName == "")  return;
	lastPath = QFileInfo(projectName).absolutePath() + QDir::separator();

	if(!projectFiles.definitionSaved)
		this->saveProject();

	//Load up the definition
	if( !projectDefinition.Load(projectName.toStdString()) ) return;

	projectFiles.path = lastPath.toStdString();
	projectFiles.definition = QFileInfo(projectName).fileName().toStdString();
	projectFiles.definitionSaved = true;
	projectFiles.nucSegValidated = false;

	startProcess();
}

void NucleusEditor::startProcess()
{
	//Assumes Image is already loaded up:
	if(!myImg) return;

	if(projectFiles.log == "")
		this->createDefaultLogName();

	//Set up a new processor:
	pProc = new ftk::ProjectProcessor();
	pProc->SetPath( lastPath.toStdString() );
	pProc->SetInputImage(myImg);
	if(labImg)
		pProc->SetOutputImage(labImg);
	if(table)
		pProc->SetTable(table);
	pProc->SetDefinition(&projectDefinition);
	pProc->Initialize();

	//Set up the process toolbar:
	processContinue->setEnabled(false);
	processToolbar->setVisible(true);
	processTaskLabel->setText(tr("  Processing Project: ") + tr(projectDefinition.name.c_str()) );
	abortProcessFlag = false;
	continueProcessFlag = false;
	menusEnabled(false);

	//start a new thread for the process:
	processThread = new ProcessThread(pProc);
	connect(processThread, SIGNAL(finished()), this, SLOT(process()));
	processThread->start();
}

//This slot keeps the processing going and updates the processToolbar
void NucleusEditor::process()
{
	//Check to see if abort has been clicked:
	if(abortProcessFlag)
	{
		abortProcessFlag = false;
		deleteProcess();
		QApplication::restoreOverrideCursor();
		return;
	}
	//Check to see if continue process has been clicked:
	if(continueProcessFlag)
	{
		processProgress->setRange(0,0);
		processTaskLabel->setText(tr("  Processing Project: ") + tr(projectDefinition.name.c_str()) );
		if(pProc->ReadyToEdit())	//Means I was editing
		{
			stopEditing();
			//this->saveProject();	//Will create default filenames and save the files
		}
		continueProcessFlag = false;
	}

	//Means I just finished nuclear segmentation and editing is not done:
	if(pProc->ReadyToEdit() && projectFiles.nucSegValidated == false)
	{
		selection->clear();
		labImg = pProc->GetOutputImage();
		segView->SetLabelImage(labImg,selection);
		this->updateNucSeg();
		table = pProc->GetTable();
		NucAdjTable = pProc->GetNucAdjTable();
		segView->SetNucAdjTable(NucAdjTable);
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		projectFiles.definitionSaved = false;

		processAbort->setEnabled(false);

		if(pProc->DoneProcessing())
		{
			NucAdjTable = pProc->GetNucAdjTable();
			segView->SetNucAdjTable(NucAdjTable);
			CellAdjTable = pProc->GetCellAdjTable();
			segView->SetCellAdjTable(CellAdjTable);
			deleteProcess();
		}
		else
		{
			processContinue->setEnabled(true);
			processTaskLabel->setText(tr("  You May Edit Nuclear Segmentation"));
			processProgress->setRange(0,2);
			processProgress->setValue(1);
		}

		this->closeViews();
		CreateNewTableWindow();
		CreateNewPlotWindow();
		this->startEditing();
	}
	else if(pProc->DoneProcessing()) //All done, so get results and clean up:
	{
		selection->clear();
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		projectFiles.definitionSaved = false;
		NucAdjTable = pProc->GetNucAdjTable();
		segView->SetNucAdjTable(NucAdjTable);
		CellAdjTable = pProc->GetCellAdjTable();
		segView->SetCellAdjTable(CellAdjTable);

		deleteProcess();
		this->updateNucSeg();
		this->updateViews();
	}
	else		//Not done, so continue:
	{
		//I need some type of input for the processor to continue:
		if( pProc->NeedInput() )
		{
			switch( pProc->NeedInput() )
			{
			case 3:			//Need something for association rules
				this->CreateDefaultAssociationRules();
				break;
			}
		}
		processThread->start();
	}
}

//Call this slot when the abort process button is clicked
void NucleusEditor::abortProcess()
{
	if(processThread)
	{
		QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
		processTaskLabel->setText(tr("  Aborting Process"));
		abortProcessFlag = true;
	}
}

//Call this slot when the continue process button is clicked
void NucleusEditor::continueProcess()
{
	if(processThread)
	{
		processAbort->setEnabled(true);
		processContinue->setEnabled(false);
		continueProcessFlag = true;
		process();
	}
}

void NucleusEditor::deleteProcess()
{
	if(processThread)
	{
		delete processThread;
		processThread = NULL;
	}
	if(pProc)
	{
		delete pProc;
		pProc = NULL;
	}
	processToolbar->setVisible(false);
	menusEnabled(true);
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
// Thread for running the segmentation algorithm:
//***********************************************************************************
ProcessThread::ProcessThread(ftk::ProjectProcessor *proc)
: QThread()
{
	myProc = proc;
}

void ProcessThread::run()
{
	if(!myProc->DoneProcessing())
		myProc->ProcessNext();
}

//***************************************************************************
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

	//autoButton = new QRadioButton(tr("Automatic Parameter Selection"),this);
	//autoButton->setChecked(true);

	//fileButton = new QRadioButton(tr("Use Parameter File..."),this);

	//fileCombo = new QComboBox();
	//fileCombo->addItem(tr(""));
	//fileCombo->addItem(tr("Browse..."));
	//connect(fileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	QHBoxLayout *bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addLayout(chLayout);
	//layout->addWidget(autoButton);
	//layout->addWidget(fileButton);
	//layout->addWidget(fileCombo);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Channel"));

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

//***************************************************************************
//***********************************************************************************
//***********************************************************************************
// A dialog to get the paramaters file to use and specify the channel if image has
// more than one:
//***********************************************************************************
QueryDialog::QueryDialog(int QueryType, QVector<QString> classes, QWidget *parent)
: QDialog(parent)
{
	idLabel = new QLabel("Choose Object IDs ('0' for all IDs) :");
	IDs = new QLineEdit();
	IDs->setMinimumWidth(100);
	IDs->setFocusPolicy(Qt::StrongFocus);
	idLayout = new QHBoxLayout;
	idLayout->addWidget(idLabel);
	idLayout->addWidget(IDs);

	if(QueryType == 1)
	{
		kLabel = new QLabel("Choose Number Of Nearest Neighbors K: ");
		K = new QLineEdit();
		K->setMinimumWidth(50);
		K->setFocusPolicy(Qt::StrongFocus);
		kLayout = new QHBoxLayout;
		kLayout->addWidget(kLabel);
		kLayout->addWidget(K);
	}
	else if(QueryType == 2)
	{
		radLabel = new QLabel("Choose Radius: ");
		RAD = new QLineEdit();
		RAD->setMinimumWidth(50);
		RAD->setFocusPolicy(Qt::StrongFocus);
		radLayout = new QHBoxLayout;
		radLayout->addWidget(radLabel);
		radLayout->addWidget(RAD);
	}	

	classLabel1 = new QLabel("Choose Source Class: ");
	classCombo1 = new QComboBox();
	classCombo1->addItem("All");
	for(int v = 0; v<classes.size(); ++v)
	{
		classCombo1->addItem(classes.at(v));
	}
	QHBoxLayout *classLayout1 = new QHBoxLayout;
	classLayout1->addWidget(classLabel1);
	classLayout1->addWidget(classCombo1);

	classLabel2 = new QLabel("Choose Destination Class: ");
	classCombo2 = new QComboBox();
	classCombo2->addItem("All");
	for(int v = 0; v<classes.size(); ++v)
	{
		classCombo2->addItem(classes.at(v));
	}
	QHBoxLayout *classLayout2 = new QHBoxLayout;
	classLayout2->addWidget(classLabel2);
	classLayout2->addWidget(classCombo2);

	//autoButton = new QRadioButton(tr("Automatic Parameter Selection"),this);
	//autoButton->setChecked(true);

	//fileButton = new QRadioButton(tr("Use Parameter File..."),this);

	//fileCombo = new QComboBox();
	//fileCombo->addItem(tr(""));
	//fileCombo->addItem(tr("Browse..."));
	//connect(fileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	layout = new QVBoxLayout;
	layout->addLayout(idLayout);
	if(QueryType == 1)
		layout->addLayout(kLayout);
	if(QueryType == 2)
		layout->addLayout(radLayout);
	layout->addLayout(classLayout1);
	layout->addLayout(classLayout2);
	//layout->addWidget(autoButton);
	//layout->addWidget(quitButton);
	//layout->addWidget(okButton);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	if(QueryType == 1)
		this->setWindowTitle(tr("K Nearest Neighbors"));
	else if(QueryType == 2)
		this->setWindowTitle(tr("Neighbors Within Radius"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

std::vector<unsigned int> QueryDialog::parseIDs(void)
{
	std::vector<unsigned int> ids;
	QString input = IDs->displayText();
	QStringList values = input.split(",");
	for(int i=0; i<values.size(); ++i)
	{
		QString str = values.at(i);
		unsigned int v = str.toUInt();
		ids.push_back( v );
	}
	return ids;
}

unsigned int QueryDialog::parseK(void)
{
	unsigned int k;
	QString input = K->displayText();
	k = input.toUInt();
	return k;
}

double QueryDialog::parseRad(void)
{
	double rad;
	QString input = RAD->displayText();
	rad = input.toUInt();
	return rad;
}

unsigned short QueryDialog::getSourceClass()
{
	return classCombo1->currentIndex();
}

unsigned short QueryDialog::getDestClass()
{
	return classCombo2->currentIndex();
}



//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

//**
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
// PRE-PROCESSING FUNCTIONS:
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
void NucleusEditor::preprocessImage(void)
{
	if(!myImg)
		return;

	int nucChannel = 0;
	if(myImg->GetImageInfo()->numChannels > 1)
		nucChannel = this->requestChannel(myImg);

	PreprocessDialog * dialog = new PreprocessDialog(lastPath,this);
	dialog->SetImage( myImg->GetItkPtr<unsigned char>(0,nucChannel) );

	if(dialog->exec())
	{
		const ftk::Image::Info * info = myImg->GetImageInfo();
		int bpChunk1 = info->numZSlices * info->numRows * info->numColumns * info->bytesPerPix;

		ftk::Preprocess::ImageType3D::Pointer img = dialog->GetImage();
		ftk::Preprocess::ImageType3D::SizeType size = img->GetLargestPossibleRegion().GetSize();
		int bpChunk2 = size[2] * size[1] * size[0];

		if( bpChunk1 == bpChunk2 )	//Can't handle resizing in the GUI quite yet.
		{
		memcpy( myImg->GetDataPtr(0,nucChannel), dialog->GetImage()->GetBufferPointer(), bpChunk2 );
		segView->update();
		projectFiles.inputSaved = false;
		}
	}

	delete dialog;
}
void NucleusEditor::CropToRegion(void)
{

}

void NucleusEditor::InvertIntensities(void)
{
	this->preprocess("Invert");
}

void NucleusEditor::MedianFilter()
{
	this->preprocess("Median");
}

void NucleusEditor::AnisotropicDiffusion()
{
	this->preprocess("GAnisotropicDiffusion");
}

void NucleusEditor::CurvAnisotropicDiffusion()
{
	this->preprocess("CAnisotropicDiffusion");
}

void NucleusEditor::GrayscaleErode()
{
	this->preprocess("GSErode");
}

void NucleusEditor::GrayscaleDilate()
{
	this->preprocess("GSDilate");
}

void NucleusEditor::GrayscaleOpen()
{
	this->preprocess("GSOpen");
}

void NucleusEditor::GrayscaleClose()
{
	this->preprocess("GSClose");
}

void NucleusEditor::SigmoidFilter()
{
	this->preprocess("Sigmoid");
}

void NucleusEditor::preprocess(QString id)
{
	if(!myImg)
		return;

 	ftkPreprocessDialog *dialog = new ftkPreprocessDialog(this->getChannelStrings(), id.toStdString(), this->myImg, this);
  	if(dialog->exec())
 	{
		segView->update();
		projectFiles.inputSaved = false;
 	}
}


//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
// KPLS Prediction Dialog:
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//*******************************************************************************************************************************
//A dialog if there are more than one training columns to choose from
PredictionDialog::PredictionDialog(QVector<QString> training_fields, QWidget *parent)
: QDialog(parent){
	channelLabel = new QLabel("Choose Training Set: ");
	channelCombo = new QComboBox();
	
	for(int v = 0; v<training_fields.size(); ++v){
		channelCombo->addItem(training_fields.at(v));
	}
	
	QGridLayout *layout = new QGridLayout;
	this->setLayout(layout);
	this->setWindowTitle(tr("Select Classifier"));

	layout->addWidget(channelLabel,0,0);
	layout->addWidget(channelCombo,0,1);

	cancelButton = new QPushButton(tr("Cancel"),this);
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
	layout->addWidget(cancelButton,10,1);
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	layout->addWidget(okButton,10,2);
}

int PredictionDialog::getTrainNumber()
{
	return channelCombo->currentIndex();
}
