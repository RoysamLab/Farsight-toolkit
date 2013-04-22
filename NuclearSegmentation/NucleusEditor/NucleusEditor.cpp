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
#include "../exe/SomaExtraction.h"
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
	AL = new ALforNucEd();
	connect(segView, SIGNAL(mouseAt(int,int,int, int,std::list<int>)), this, SLOT(setMouseStatus(int,int,int, int, std::list<int>)));
	connect(segView, SIGNAL(autoMerge()), this, SLOT(mergeCells()));
	connect(segView, SIGNAL(emitTimeChanged()), this, SLOT(update5DTable()));
	connect(AL, SIGNAL(Classification_Done()), this, SLOT(ExtractClassificationResult()));
	selection = new ObjectSelection();
	connect(selection, SIGNAL(MultiChanged()), this, SLOT(updateMultiLabels()));
	this->setCentralWidget(segView);

  #ifdef USE_QT_TESTING
  this->TestInputFile = "";
  this->TestBaselineImageFileName = "";
  this->Tester = new GUITester(this);
  #endif

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
	renWin.clear();
	//	hisWin.clear();
	pWizard=NULL;

	myImg = NULL;
	labImg = NULL;

	nucSeg = NULL;
	pProc = NULL;
	diffusion_map = NULL;
	processThread = NULL;
	table = NULL;
	NucAdjTable = NULL;
	CellAdjTable = NULL;
	this->HeatmapWin = NULL;

#ifdef USE_TRACKING
	mfcellTracker = NULL;
#endif

	kplsRun = 0;	//This flag gets set after kpls has run to make sure we show the colored centroids!!
	trainName = 0;
	predictName = 0;
	activeRun = 0;
	saveSettingsOnExit = true;

	this->resize(800,800);

	this->readSettings();
	raise();
	activateWindow();
	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );

  //parse command line arguments (if any)
  QStringList args = QCoreApplication::arguments();
  if(args.size() > 1)
    {
    this->parseArguments(args);
    }
}

NucleusEditor::~NucleusEditor()
{
	this->closeViews();
	if(selection) delete selection;
	if(nucSeg) delete nucSeg;
	if(pProc) delete pProc;
	delete AL;
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

void NucleusEditor::clearSettings()
{
	int reply = QMessageBox::question(this, tr("Clear settings"),
		"This action will delete all custom Nucleus Editor settings, reverting them back\nto their default values.  Are you sure you'd like to do this?",
		QMessageBox::Yes |QMessageBox::No, QMessageBox::No);

	if (reply == QMessageBox::Yes)
	{
		QSettings settings;
		settings.clear();
		segView->clearSettings();
		saveSettingsOnExit = false;
	}
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
	processToolbar->setObjectName("processToolbar");
	processToolbar->setVisible(false);

	processAbort = new QAction(QIcon(":/icons/abort.png"), tr("ABORT"), this);
	processAbort->setObjectName("processAbort");
	processAbort->setToolTip(tr("Abort Processing"));
	processAbort->setEnabled(true);
	connect(processAbort, SIGNAL(triggered()), this, SLOT(abortProcess()));
	processToolbar->addAction(processAbort);

	processContinue = new QAction(QIcon(":/icons/go.png"), tr("GO"), this);
	processAbort->setObjectName("processAbort");
	processContinue->setToolTip(tr("Continue Processing"));
	processContinue->setEnabled(false);
	connect(processContinue, SIGNAL(triggered()),this, SLOT(continueProcess()));
	processToolbar->addAction(processContinue);

	processProgress = new QProgressBar();
	processProgress->setObjectName("processProgress");
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
/*	   Finally, after all these action's, we bind them to a "QActionGroup".   */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NucleusEditor::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->setObjectName("fileMenu");

	loadProjectAction = new QAction(tr("Load Project..."), this);
	loadProjectAction->setObjectName("loadProjectAction");
	loadProjectAction->setStatusTip(tr("Load ftk project"));
	loadProjectAction->setShortcut(tr("Ctrl+L"));
	connect(loadProjectAction, SIGNAL(triggered()), this, SLOT(loadProject()));
	fileMenu->addAction(loadProjectAction);

	loadImageAction = new QAction(tr("Load Image..."), this);
	loadImageAction->setObjectName("loadImageAction");
	loadImageAction->setStatusTip(tr("Load an image into the 5D image browser"));
	loadImageAction->setShortcut(tr("Ctrl+O"));
	connect(loadImageAction, SIGNAL(triggered()), this, SLOT(askLoadImage()));
	fileMenu->addAction(loadImageAction);

	//Amin
	load5DImageAction = new QAction(tr("Load 5D Images"), this);
	load5DImageAction->setObjectName("load5DImageAction");
	connect(load5DImageAction, SIGNAL(triggered()), this, SLOT(askload5DImage()));
	fileMenu->addAction(load5DImageAction);

	load5DLabelImageAction = new QAction(tr("Load 5D Labels"), this);
	load5DLabelImageAction->setObjectName("load5DLabelImageAction");
	connect(load5DLabelImageAction, SIGNAL(triggered()), this, SLOT(askload5DLabelImage()));
	fileMenu->addAction(load5DLabelImageAction);


	loadLabelAction = new QAction(tr("Load Result..."), this);
	loadLabelAction->setObjectName("loadLabelAction");
	loadLabelAction->setStatusTip(tr("Load a result image into the image browser"));
	connect(loadLabelAction,SIGNAL(triggered()), this, SLOT(askLoadResult()));
	fileMenu->addAction(loadLabelAction);

	loadTableAction = new QAction(tr("Load Table..."), this);
	loadTableAction->setObjectName("loadTableAction");
	loadTableAction->setStatusTip(tr("Load data table from text file"));
	connect(loadTableAction, SIGNAL(triggered()), this, SLOT(askLoadTable()));
	fileMenu->addAction(loadTableAction);

	fileMenu->addSeparator();

	processProjectAction = new QAction(tr("Process Image..."), this);
	processProjectAction->setObjectName("processProjectAction");
	processProjectAction->setStatusTip(tr("Choose project definition file to use in processing the current image"));
	connect(processProjectAction, SIGNAL(triggered()), this, SLOT(processProject()));
	fileMenu->addAction(processProjectAction);

	fileMenu->addSeparator();

	saveProjectAction = new QAction(tr("Save Project..."), this);
	saveProjectAction->setObjectName("saveProjectAction");
	saveProjectAction->setStatusTip(tr("Save the active project files..."));
	saveProjectAction->setShortcut(tr("Ctrl+S"));
	connect(saveProjectAction, SIGNAL(triggered()), this, SLOT(saveProject()));
	fileMenu->addAction(saveProjectAction);

	saveSomaImageAction = new QAction(tr("Save Soma Image..."), this);
	saveSomaImageAction->setObjectName("saveSomaImageAction");
	connect(saveSomaImageAction, SIGNAL(triggered()), this, SLOT(saveSomaImage()));
	fileMenu->addAction(saveSomaImageAction);
    
    saveNeuronImageAction = new QAction(tr("Save Neuron Image..."), this);
	saveNeuronImageAction->setObjectName("saveNeuronImageAction");
	connect(saveNeuronImageAction, SIGNAL(triggered()), this, SLOT(saveNeuronImage()));
	fileMenu->addAction(saveNeuronImageAction);

	saveImageAction = new QAction(tr("Save Image..."), this);
	saveImageAction->setObjectName("saveImageAction");
	saveImageAction->setStatusTip(tr("Save image currently displayed"));
	connect(saveImageAction, SIGNAL(triggered()), this, SLOT(askSaveImage()));
	fileMenu->addAction(saveImageAction);

	saveLabelAction = new QAction(tr("Save Result..."), this);
	saveLabelAction->setObjectName("saveLabelAction");
	saveLabelAction->setStatusTip(tr("Save a segmentation result image"));
	connect(saveLabelAction, SIGNAL(triggered()), this, SLOT(askSaveResult()));
	fileMenu->addAction(saveLabelAction);

	saveTableAction = new QAction(tr("Save Table..."), this);
	saveTableAction->setObjectName("saveTableAction");
	saveTableAction->setStatusTip(tr("Save the features table"));
	connect(saveTableAction, SIGNAL(triggered()), this, SLOT(askSaveTable()));
	fileMenu->addAction(saveTableAction);

	saveDisplayAction = new QAction(tr("Save Display Image..."), this);
	saveDisplayAction->setObjectName("saveDisplayAction");
	saveDisplayAction->setStatusTip(tr("Save displayed image to file"));
	connect(saveDisplayAction, SIGNAL(triggered()), this, SLOT(saveDisplayImageToFile()));
	fileMenu->addAction(saveDisplayAction);

	saveCompositeAction = new QAction(tr("Save Composite Image..."), this);
	saveCompositeAction->setObjectName("saveCompositeAction");
	saveCompositeAction->setStatusTip(tr("Save displayed composite image to file"));
	connect(saveCompositeAction, SIGNAL(triggered()), this, SLOT(saveCompositeImageToFile()));
	fileMenu->addAction(saveCompositeAction);

	fileMenu->addSeparator();

	exitAction = new QAction(tr("Exit"), this);
	exitAction->setObjectName("exitAction");
	exitAction->setShortcut(tr("Ctrl+Q"));
	exitAction->setStatusTip(tr("Exit the application"));
	connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
	fileMenu->addAction(exitAction);

	//VIEW MENU
	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->setObjectName("viewMenu");

	setPreferencesAction = new QAction(tr("Preferences..."), this);
	setPreferencesAction->setObjectName("setPreferencesAction");
	connect(setPreferencesAction, SIGNAL(triggered()), this, SLOT(setPreferences()));
	viewMenu->addAction(setPreferencesAction);

	viewMenu->addSeparator();

	showCrosshairsAction = new QAction(tr("Show Selection Crosshairs"), this);
	showCrosshairsAction->setObjectName("showCrosshairsAction");
	showCrosshairsAction->setCheckable(true);
	showCrosshairsAction->setChecked( segView->GetCrosshairsVisible() );
	showCrosshairsAction->setStatusTip(tr("Show Crosshairs at selected object"));
	connect(showCrosshairsAction, SIGNAL(triggered()), this, SLOT(toggleCrosshairs()));
	viewMenu->addAction(showCrosshairsAction);

	showBoundsAction = new QAction(tr("Show &Boundaries"), this);
	showBoundsAction->setObjectName("showBoundsAction");
	showBoundsAction->setCheckable(true);
	showBoundsAction->setChecked( segView->GetBoundsVisible() );
	showBoundsAction->setStatusTip(tr("Draw boundaries using a label image"));
	showBoundsAction->setShortcut(tr("Ctrl+B"));
	connect(showBoundsAction, SIGNAL(triggered()), this, SLOT(toggleBounds()));
	viewMenu->addAction(showBoundsAction);

	showIDsAction = new QAction(tr("Show Object IDs"), this);
	showIDsAction->setObjectName("showIDsAction");
	showIDsAction->setCheckable(true);
	showIDsAction->setChecked( segView->GetIDsVisible() );
	showIDsAction->setStatusTip(tr("Show IDs of objects"));
	showIDsAction->setShortcut(tr("Ctrl+I"));
	connect(showIDsAction, SIGNAL(triggered()), this, SLOT(toggleIDs()));
	viewMenu->addAction(showIDsAction);

	showCentroidsAction = new QAction(tr("Show Object Centroids"), this);
	showCentroidsAction->setObjectName("showCentroidsAction");
	showCentroidsAction->setCheckable(true);
	showCentroidsAction->setChecked( segView->GetCentroidsVisible() );
	showCentroidsAction->setStatusTip(tr("Show Centroids of objects"));
	//showCentroidsAction->setShortcut(tr(""));
	connect(showCentroidsAction, SIGNAL(triggered()), this, SLOT(toggleCentroids()));
	viewMenu->addAction(showCentroidsAction);

	adjacencyMenu = viewMenu->addMenu(tr("Show Adjacency"));
	adjacencyMenu->setObjectName("adjacencyMenu");

	showNucAdjAction = new QAction(tr("Nuclear"), this);
	showNucAdjAction->setObjectName("showNucAdjAction");
	showNucAdjAction->setCheckable(true);
	showNucAdjAction->setChecked( segView->GetNucAdjVisible() );
	showNucAdjAction->setStatusTip(tr("Show adjacency of nuclei"));
	showNucAdjAction->setShortcut(tr("Shift+N"));
	connect(showNucAdjAction, SIGNAL(triggered()), this, SLOT(toggleNucAdjacency()));
	adjacencyMenu->addAction(showNucAdjAction);

	showCellAdjAction = new QAction(tr("Cellular"), this);
	showCellAdjAction->setObjectName("showCellAdjAction");
	showCellAdjAction->setCheckable(true);
	showCellAdjAction->setChecked( segView->GetCellAdjVisible() );
	showCellAdjAction->setStatusTip(tr("Show adjacency of cells"));
	showCellAdjAction->setShortcut(tr("Shift+C"));
	connect(showCellAdjAction, SIGNAL(triggered()), this, SLOT(toggleCellAdjacency()));
	adjacencyMenu->addAction(showCellAdjAction);

	zoomMenu = viewMenu->addMenu(tr("Zoom"));
	zoomMenu->setObjectName("zoomMenu");

	zoomInAction = new QAction(tr("Zoom In"), this);
	zoomInAction->setObjectName("zoomInAction");
	zoomInAction->setStatusTip(tr("Zoom In On The Displayed Image"));
	zoomInAction->setShortcut(tr("="));
	connect(zoomInAction, SIGNAL(triggered()), segView, SLOT(zoomIn()));
	zoomMenu->addAction(zoomInAction);

	zoomOutAction  = new QAction(tr("Zoom Out"), this);
	zoomOutAction->setObjectName("zoomOutAction");
	zoomOutAction->setStatusTip(tr("Zoom Out of The Displayed Image"));
	zoomOutAction->setShortcut(tr("-"));
	connect(zoomOutAction, SIGNAL(triggered()), segView, SLOT(zoomOut()));
	zoomMenu->addAction(zoomOutAction);

	displayChannelMenu = viewMenu->addMenu(tr("Display Channel"));
	displayChannelMenu->setObjectName("displayChannelMenu");
	connect(displayChannelMenu, SIGNAL(aboutToShow()), this, SLOT(DisplayChannelsMenu()));

	viewMenu->addSeparator();

	newTableAction = new QAction(tr("New Table"), this);
	newTableAction->setObjectName("newTableAction");
	newTableAction->setStatusTip(tr("Open a new table window"));
	connect(newTableAction, SIGNAL(triggered()), this, SLOT(CreateNewTableWindow()));
	viewMenu->addAction(newTableAction);

	newScatterAction = new QAction(tr("New Scatterplot"), this);
	newScatterAction->setObjectName("newScatterAction");
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

	newHistoAction = new QAction(tr("New Histogram"),this);
	newHistoAction->setObjectName("newHistoAction");
	newHistoAction->setStatusTip(tr("Open a new Histogram Window"));
	connect(newHistoAction,SIGNAL(triggered()),this,SLOT(CreateNewHistoWindow()));
	viewMenu->addAction(newHistoAction);

	newRenderAction = new QAction(tr("New Render Window"),this);
	newRenderAction->setObjectName("newRenderAction");
	newRenderAction->setStatusTip(tr("Open a new Render Window"));
	connect(newRenderAction,SIGNAL(triggered()),this,SLOT(CreateNewRenderWindow()));
	viewMenu->addAction(newRenderAction);

	//ragMenu = viewMenu->addMenu(tr("New Region Adjacency Graph"));

	//nucRagAction = new QAction(tr("Nuclear Adjacency Graph"), this);	
	//connect(nucRagAction, SIGNAL(triggered()), this, SLOT(CreateNewNucRAG()));
	//ragMenu->addAction(nucRagAction);

	//cellRagAction = new QAction(tr("Cellular Adjacency Graph"), this);
	//connect(cellRagAction, SIGNAL(triggered()), this, SLOT(CreateNewCellRAG()));
	//ragMenu->addAction(cellRagAction);

	imageIntensityAction = new QAction(tr("Adjust Image Intensity"), this);
	imageIntensityAction->setObjectName("imageIntensityAction");
	imageIntensityAction->setStatusTip(tr("Allows modification of image intensity"));
	connect(imageIntensityAction, SIGNAL(triggered()), segView, SLOT(AdjustImageIntensity()));
	viewMenu->addAction(imageIntensityAction);

	//TOOL MENU
	toolMenu = menuBar()->addMenu(tr("Tools"));
	toolMenu->setObjectName("toolMenu");

	getCentroidAction = new QAction(tr("Get Object Centroids"), this);
	getCentroidAction->setObjectName("getCentroidAction");
	connect(getCentroidAction, SIGNAL(triggered()), this, SLOT(getCentroids()));
	toolMenu->addAction(getCentroidAction);

	roiMenu = toolMenu->addMenu(tr("Region Of Interest"));
	roiMenu->setObjectName("roiMenu");

	drawROIAction = new QAction(tr("Draw ROI"), this);
	drawROIAction->setObjectName("drawROIAction");
	connect(drawROIAction, SIGNAL(triggered()), this, SLOT(startROI()));
	roiMenu->addAction(drawROIAction);

	clearROIAction = new QAction(tr("Clear ROI"), this);
	clearROIAction->setObjectName("clearROIAction");
	connect(clearROIAction, SIGNAL(triggered()), this, SLOT(clearROI()));
	roiMenu->addAction(clearROIAction);

	saveROIAction = new QAction(tr("Save ROI Mask..."), this);
	saveROIAction->setObjectName("saveROIAction");
	connect(saveROIAction, SIGNAL(triggered()), this, SLOT(saveROI()));
	roiMenu->addAction(saveROIAction);

	loadROIAction = new QAction(tr("Load ROI Mask..."), this);
	loadROIAction->setObjectName("loadROIAction");
	connect(loadROIAction, SIGNAL(triggered()), this, SLOT(loadROI()));
	roiMenu->addAction(loadROIAction);

	roiStatsAction = new QAction(tr("Compute ROI Statistics"), this);
	roiStatsAction->setObjectName("roiStatsAction");
	connect(roiStatsAction, SIGNAL(triggered()), this, SLOT(roiStatistics()));
	toolMenu->addAction(roiStatsAction);

	preprocessAction = new QAction(tr("Preprocess Image..."), this);
	preprocessAction->setObjectName("preprocessAction");
	connect(preprocessAction, SIGNAL(triggered()), this, SLOT(preprocessImage()));
	toolMenu->addAction(preprocessAction);

	unmixChannelsAction = new QAction(tr("Unmix Channels..."), this);
	unmixChannelsAction->setObjectName("unmixChannelsAction");
	connect(unmixChannelsAction, SIGNAL(triggered()), this, SLOT(unmixChannels()));
	toolMenu->addAction(unmixChannelsAction);

	segmentNucleiAction = new QAction(tr("Segment Nuclei..."), this);
	segmentNucleiAction->setObjectName("segmentNucleiAction");
	connect(segmentNucleiAction, SIGNAL(triggered()), this, SLOT(segmentNuclei()));
	toolMenu->addAction(segmentNucleiAction);

	editNucleiAction = new QAction(tr("Edit Nuclei"), this);
	editNucleiAction->setObjectName("editNucleiAction");
	connect(editNucleiAction, SIGNAL(triggered()), this, SLOT(startEditing()));
	toolMenu->addAction(editNucleiAction);

	svmAction = new QAction(tr("Detect Outliers"), this);
	svmAction->setObjectName("svmAction");
	connect(svmAction, SIGNAL(triggered()), this, SLOT(startSVM()));
	toolMenu->addAction(svmAction);

	databaseAction = new QAction(tr("Update Database"), this);
	databaseAction->setObjectName("databaseAction");
	connect(databaseAction, SIGNAL(triggered()), this, SLOT(updateDatabase()));
	toolMenu->addAction(databaseAction);

	activeContourAction = new QAction(tr("Active Contour"), this);
	activeContourAction->setObjectName("activeContourAction");
	connect(activeContourAction, SIGNAL(triggered()), this, SLOT(segmentByActiveContour()));
	toolMenu->addAction(activeContourAction);

	activeMenu   = toolMenu->addMenu(tr("Active Learning"));
	activeMenu->setObjectName("activeMenu");

	activeAction = new QAction(tr("Start Active Learning"), this);
	activeAction->setObjectName("activeAction");
	connect(activeAction, SIGNAL(triggered()), this, SLOT(startActiveLearningwithFeat()));
	activeMenu->addAction(activeAction);

	saveActiveResultsAction = new QAction(tr("Save Active Learning Results (Only for Multiple Images)"), this);
	saveActiveResultsAction->setObjectName("saveActiveResultsAction");
	connect(saveActiveResultsAction, SIGNAL(triggered()), this, SLOT(SaveActiveLearningResults()));
	activeMenu->addAction(saveActiveResultsAction);

	saveActiveLearningModelAction = new QAction(tr("Save Active Learning Model"), this);
	saveActiveLearningModelAction->setObjectName("saveActiveLearningModelAction");
	connect(saveActiveLearningModelAction, SIGNAL(triggered()), this, SLOT(SaveActiveLearningModel()));
	activeMenu->addAction(saveActiveLearningModelAction);

	classifyFromActiveLearningModelAction = new QAction(tr("Load Active Learning Model"), this);
	classifyFromActiveLearningModelAction->setObjectName("classifyFromActiveLearningModelAction");
	connect(classifyFromActiveLearningModelAction, SIGNAL(triggered()), this, SLOT(classifyFromActiveLearningModel()));
	activeMenu->addAction(classifyFromActiveLearningModelAction);

	// Will be enabled and disabled separately
	//showGalleryAction->setEnabled(false);
	saveActiveResultsAction->setEnabled(false);

	classifyMenu = toolMenu->addMenu(tr("Classifier"));
	classifyMenu->setObjectName("classifyMenu");

	trainAction = new QAction(tr("Train"), this);
	trainAction->setObjectName("trainAction");
	connect(trainAction, SIGNAL(triggered()), this, SLOT(startTraining()));
	classifyMenu->addAction(trainAction);

	kplsAction = new QAction(tr("Classify"), this);
	kplsAction->setObjectName("kplsAction");
	connect(kplsAction, SIGNAL(triggered()), this, SLOT(startKPLS()));
	classifyMenu->addAction(kplsAction);

	validationMenu = toolMenu->addMenu(tr("Validation"));
	validationMenu->setObjectName("validationMenu");
	activeValidation = new QAction(tr("Start Active Validation"), this);
	activeValidation->setObjectName("activeValidation");
	connect(activeValidation, SIGNAL(triggered()), this, SLOT(startActiveValidation()));
	validationMenu->addAction(activeValidation);

	runClusAction = new QAction(tr("Clustering Heatmap"), this);
	runClusAction->setObjectName("runClusAction");
	connect(runClusAction, SIGNAL(triggered()), this, SLOT(runClus()));
	toolMenu->addAction(runClusAction);


	//EDITING MENU
	editMenu = menuBar()->addMenu(tr("&Editing"));
	editMenu->setObjectName("editMenu");

	clearSelectAction = new QAction(tr("Clear Selections"), this);
	clearSelectAction->setObjectName("clearSelectAction");
	clearSelectAction->setStatusTip(tr("Clear Current Object Selections"));
	clearSelectAction->setShortcut(tr("Ctrl+C"));
	connect(clearSelectAction, SIGNAL(triggered()), this, SLOT(clearSelections()));
	editMenu->addAction(clearSelectAction);

	visitAction = new QAction(tr("Mark as Visited"), this);
	visitAction->setObjectName("visitAction");
	visitAction->setStatusTip(tr("Mark Selected Objects as Visited"));
	visitAction->setShortcut(tr("Ctrl+V"));
	connect(visitAction, SIGNAL(triggered()), this, SLOT(markVisited()));
	//editMenu->addAction(visitAction);

	editMenu->addSeparator();

	classAction = new QAction(tr("Change Class"), this);
	classAction->setObjectName("classAction");
	classAction->setStatusTip(tr("Modify the class designation for the selected objects"));
	connect(classAction, SIGNAL(triggered()), this, SLOT(changeClass()));
	editMenu->addAction(classAction);

	addAction = new QAction(tr("Add Object"), this);
	addAction->setObjectName("addAction");
	addAction->setStatusTip(tr("Draw a Box to add a new object"));
	addAction->setShortcut(tr("Ctrl+A"));
	connect(addAction, SIGNAL(triggered()), segView, SLOT(GetBox()));
	connect(segView, SIGNAL(boxDrawn(int,int,int,int,int)), this, SLOT(addCell(int,int,int,int,int)));
	editMenu->addAction(addAction);

	mergeAction = new QAction(tr("Merge Objects"), this);
	mergeAction->setObjectName("mergeAction");
	mergeAction->setStatusTip(tr("Merge Objects"));
	mergeAction->setShortcut(tr("Ctrl+M"));
	connect(mergeAction, SIGNAL(triggered()), this, SLOT(mergeCells()));
	editMenu->addAction(mergeAction);

	deleteAction = new QAction(tr("Delete Objects"), this);
	deleteAction->setObjectName("deleteAction");
	deleteAction->setStatusTip(tr("Deletes the selected objects"));
	deleteAction->setShortcut(tr("Ctrl+D"));
	connect(deleteAction,SIGNAL(triggered()),this,SLOT(deleteCells()));
	editMenu->addAction(deleteAction);

	fillAction = new QAction(tr("Fill Objects"), this);
	fillAction->setObjectName("fillAction");
	fillAction->setStatusTip(tr("Fill holes in the selected objects"));
	fillAction->setShortcut(tr("Ctrl+F"));
	connect(fillAction,SIGNAL(triggered()),this,SLOT(fillCells()));	
	editMenu->addAction(fillAction);

	splitZAction = new QAction(tr("Split Objects At Z"), this);
	splitZAction->setObjectName("splitZAction");
	splitZAction->setStatusTip(tr("Split selected objects along the current Z slice"));
	splitZAction->setShortcut(tr("Ctrl+T"));
	connect(splitZAction, SIGNAL(triggered()), this, SLOT(splitCellAlongZ()));
	editMenu->addAction(splitZAction);


	splitAction = new QAction(tr("Split Objects X-Y"), this);
	splitAction->setObjectName("splitAction");
	splitAction->setStatusTip(tr("Split objects by choosing two seed points"));
	splitAction->setShortcut(tr("Ctrl+P"));
	splitAction->setCheckable(true);
	splitAction->setChecked(false);
	connect(splitAction, SIGNAL(triggered()), this, SLOT(splitCells()));
	connect(segView, SIGNAL(pointsClicked(int,int,int,int,int,int)), this, SLOT(splitCell(int,int,int,int,int,int)));
	editMenu->addAction(splitAction);

	batchSplitAction = new QAction(tr("Batch Split Objects "), this);
	batchSplitAction->setObjectName("batchSplitAction");
	batchSplitAction->setStatusTip(tr("Batch Split objects by selecting the objects"));
	batchSplitAction->setShortcut(tr("Ctrl+Shift+P"));
	connect(batchSplitAction, SIGNAL(triggered()), this, SLOT(batchSplitCells()));
	editMenu->addAction(batchSplitAction);

	editMenu->addSeparator();

	exclusionAction = new QAction(tr("Apply Exclusion Margin..."), this);
	exclusionAction->setObjectName("exclusionAction");
	exclusionAction->setStatusTip(tr("Set parameters for exclusion margin"));
	connect(exclusionAction, SIGNAL(triggered()), this, SLOT(applyExclusionMargin()));
	editMenu->addAction(exclusionAction);


	// MODELS MENU
	modelsMenu = menuBar()->addMenu(tr("&Models"));
	modelsMenu->setObjectName("modelsMenu");

	createTrainingAction = new QAction(tr("Create Training Model..."), this);
	createTrainingAction->setObjectName("createTrainingAction");
	connect(createTrainingAction, SIGNAL(triggered()), this, SLOT(createTrainer()));
	modelsMenu->addAction(createTrainingAction);

	appendTrainingAction = new QAction(tr("Append Training Model..."), this);
	appendTrainingAction->setObjectName("appendTrainingAction");
	connect(appendTrainingAction, SIGNAL(triggered()), this, SLOT(appendTrainer()));
	modelsMenu->addAction(appendTrainingAction);

	//QUERIES MENU
	queriesMenu = menuBar()->addMenu(tr("Queries"));
	queriesMenu->setObjectName("queriesMenu");

	kNearestNeighborsAction = new QAction(tr("Query K Nearest Neighbors..."), this);
	kNearestNeighborsAction->setObjectName("kNearestNeighborsAction");
	connect(kNearestNeighborsAction, SIGNAL(triggered()), this, SLOT(queryKNearest()));
	queriesMenu->addAction(kNearestNeighborsAction);

	inRadiusNeighborsAction = new QAction(tr("Query Neighbors Within Radius..."), this);
	inRadiusNeighborsAction->setObjectName("inRadiusNeighborsAction");
	connect(inRadiusNeighborsAction, SIGNAL(triggered()), this, SLOT(queryInRadius()));
	queriesMenu->addAction(inRadiusNeighborsAction);

	computeDiffusionMapAction = new QAction(tr("Compute Diffusion Map..."), this);
	computeDiffusionMapAction->setObjectName("computeDiffusionMapAction");
	connect(computeDiffusionMapAction, SIGNAL(triggered()), this, SLOT(computeDistanceMap()));
	queriesMenu->addAction(computeDiffusionMapAction);

	kNearestDiffusionAction = new QAction(tr("Query K Nearest Diffusion Neighbors..."), this);
	kNearestDiffusionAction->setObjectName("kNearestDiffusionAction");
	connect(kNearestDiffusionAction, SIGNAL(triggered()), this, SLOT(queryKDiffusion()));
	queriesMenu->addAction(kNearestDiffusionAction);

	queryViewsOffAction = new QAction(tr("Set Query Views Off"), this);
	queryViewsOffAction->setObjectName("queryViewsOffAction");
	queryViewsOffAction->setShortcut(tr("Shift+O"));
	connect(queryViewsOffAction, SIGNAL(triggered()), this, SLOT(queryViewsOff()));
	queriesMenu->addAction(queryViewsOffAction);

#ifdef USE_TRACKING
	//5-D MENU
	fiveDMenu = menuBar()->addMenu(tr("5D"));
	fiveDMenu->setObjectName("fiveDMenu");

	trackingAction = new QAction(tr("Tracking..."),this);
	trackingAction->setObjectName("trackingAction");
	connect(trackingAction, SIGNAL(triggered()), this, SLOT(startTracking()));
	fiveDMenu->addAction(trackingAction);

	kymoViewAction = new QAction(tr("Kymograph..."),this);
	kymoViewAction->setObjectName("kymoViewAction");
	connect(kymoViewAction, SIGNAL(triggered()), this, SLOT(displayKymoGraph()));
	fiveDMenu->addAction(kymoViewAction);
#endif

  #ifdef USE_QT_TESTING
  //Testing menu & actions
	testingMenu = menuBar()->addMenu(tr("Testing"));
  testingMenu->setObjectName("testingMenu");

  this->recordAction = new QAction("Record Test", this);
  this->recordAction->setObjectName("recordAction");
  this->recordAction->setStatusTip("Record a test to a .xml file");
  connect(this->recordAction, SIGNAL(triggered()), this->Tester, SLOT(record()));
  testingMenu->addAction(this->recordAction);

  this->playAction = new QAction("Play Test", this);
  this->playAction->setObjectName("playAction");
  this->playAction->setStatusTip("Run a previously recorded test");
  connect(this->playAction, SIGNAL(triggered()), this->Tester, SLOT(play()));
  testingMenu->addAction(this->playAction);

  this->clearSettingsAction = new QAction("Clear Settings", this);
  this->clearSettingsAction->setObjectName("clearSettingsAction");
  this->clearSettingsAction->setStatusTip("Revert settings back to default values");
  connect(this->clearSettingsAction, SIGNAL(triggered()), this, SLOT(clearSettings()));
  testingMenu->addAction(this->clearSettingsAction);
  #endif

	//HELP MENU
	helpMenu = menuBar()->addMenu(tr("Help"));
	helpMenu->setObjectName("helpMenu");
	aboutAction = new QAction(tr("About"),this);
	aboutAction->setObjectName("aboutAction");
	aboutAction->setStatusTip(tr("About the application"));
	connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));
	helpMenu->addAction(aboutAction);
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
	QString version = QString::number(CPACK_PACKAGE_VERSION_MAJOR);
	version += ".";
	version += QString::number(CPACK_PACKAGE_VERSION_MINOR);
	version += ".";
	version += QString::number(CPACK_PACKAGE_VERSION_PATCH);

	QString text = tr("<h2>FARSIGHT ") + version + tr("</h2>");
	text += tr("<h3>Rensselear Polytechnic Institute</h3>");
	text += tr("<a><u>http://www.farsight-toolkit.org</a></u>");

	QMessageBox::about(this, tr("About FARSIGHT"), text);
	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "About FARSIGHT" << std::endl;
	}
	#endif
}

//******************************************************************************
// SLOT: changes the status bar to say the mouse coordinates
//******************************************************************************
void NucleusEditor::setMouseStatus(int x, int y, int z, int t, std::list<int> v)
{
	QString statusMsg("X: " + QString::number(x) + ", Y: " + QString::number(y) + ", Z: " + QString::number(z) + ", T: " + QString::number(t));

	int i = 1;
	for(std::list<int>::iterator it = v.begin(); it != v.end(); it++) {
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
	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "closing NucleusEditor" << std::endl;
	}
  #endif
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
bool NucleusEditor::saveSomaImage()
{
	if(!labImg){
		std::cout<<"No label image\n";
		return false;
	}

	// Soma Montage
	typedef itk::Image<unsigned short, 3> rawImageType_uint;
	typedef itk::ImageFileWriter< rawImageType_uint > uintImageWriter;

	rawImageType_uint::Pointer imageLabelMontage = labImg->GetItkPtr<rawImageType_uint::PixelType>(0,0);
	rawImageType_uint::RegionType ImageMontageRegion = imageLabelMontage->GetLargestPossibleRegion();
	rawImageType_uint::Pointer imageSomaMontage = rawImageType_uint::New();
	rawImageType_uint::SizeType im_size = imageLabelMontage->GetBufferedRegion().GetSize();
	rawImageType_uint::RegionType region;
	region.SetSize( ImageMontageRegion.GetSize() );
	region.SetIndex( ImageMontageRegion.GetIndex() );
	imageSomaMontage->SetRegions( region );
	imageSomaMontage->Allocate();
	imageSomaMontage->FillBuffer(0);
	try
	{
		imageSomaMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	rawImageType_uint::PixelType * imageSomaArray = imageSomaMontage->GetBufferPointer();
	rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();

	itk::SizeValueType sizeXY = im_size[1] * im_size[0];
	itk::SizeValueType sizeX = im_size[0];
	std::map<unsigned int, unsigned int> classMap;
	std::string prediction_column;
	for( unsigned i=0; i<table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = table->GetColumnName(i);
		if( current_column.find("prediction") != std::string::npos ){
			prediction_column = current_column;
			break;
		}
	}

	if(prediction_column.empty()){
		std::cout<<"No prediction column\n";
		return false;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Load Label Image"), lastPath, standardImageTypes);

	std::cout<<"Writing the labels for class 1 in "<<prediction_column<<" to soma file\n"; 

	for(itk::SizeValueType row=0; row<table->GetNumberOfRows(); ++row)
	{
		classMap[table->GetValue(row,0).ToUnsignedInt()] = table->GetValueByName( row, prediction_column.c_str() ).ToInt();
	}

	#pragma omp parallel for
	for(itk::IndexValueType i=0; i<im_size[2]; ++i)
	{
		for(itk::IndexValueType j=0; j<im_size[1]; ++j)
		{
			for(itk::IndexValueType k=0; k<im_size[0]; ++k)
			{
				itk::IndexValueType offset = (i*sizeXY)+(j*sizeX)+k;
				if( imageLabelArray[offset] != 0 )
					if(classMap[imageLabelArray[offset]] == 1)
						imageSomaArray[offset] = imageLabelArray[offset];
					else
						imageSomaArray[offset] = 0;
			}
		}
	}	

	std::string nameSomaMontage = fileName.toStdString();
	uintImageWriter::Pointer writer = uintImageWriter::New();
	writer->SetFileName(nameSomaMontage);
	writer->SetInput(imageSomaMontage);
	writer->Update();

	/*
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		if(table->GetValueByName(row, "prediction_active_neu").ToInt() != 1)
		{
			table->RemoveRow(row);
			--row;
		}
	}*/


	//if(!labImg || !table) return false;

	//std::vector< vtkSmartPointer< vtkTable > > featureTable;
	//if(myImg->GetImageInfo()->numTSlices > 1)
	//	featureTable = nucSeg->table4DImage;
	//else
	//	featureTable.push_back(table);


	//bool found = false;
	//for(int col=((int)table->GetNumberOfColumns())-1; col>=0; --col)
	//{	
	//	std::string current_column = table->GetColumnName(col);
	//	if(current_column.find("prediction_active_mg") != std::string::npos )
	//	{
	//		found = true;
	//		break;			
	//	}	
	//}
	//if(!found) return false;

	//QString filename = QFileDialog::getSaveFileName(this, tr("Save Soma Image As..."),lastPath, standardImageTypes);
	//if(filename == "")
	//	return false;
	//lastPath = QFileInfo(filename).absolutePath();
	//std::string Filename = filename.toStdString();
	//string::iterator it;
	//it = Filename.end() - 4;
	//Filename.erase(it, it+4);
	//
	//typedef itk::Image<unsigned short, 3> LabelImageType;
	////typedef itk::Image<unsigned char, 3> UCharImageType;
	//typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
	//for(int t=0; t<(int)myImg->GetImageInfo()->numTSlices; ++t)
	//{
	//	std::stringstream ss;
	//	ss << t;
	//	std::string time = ss.str();
	//	std::string somasFilename = Filename + "_somas_" + time + ".tif";
	//	LabelImageType::Pointer labelImage = labImg->GetItkPtr<unsigned short>(t,0);
	//	LabelImageType::PixelType * labelArray = labelImage->GetBufferPointer();

	//	LabelImageType::Pointer somaImage = LabelImageType::New();
	//	itk::Size<3> im_size = labelImage->GetBufferedRegion().GetSize();
	//	LabelImageType::IndexType start;
	//	start[0] =   0;  // first index on X
	//	start[1] =   0;  // first index on Y    
	//	start[2] =   0;  // first index on Z  
	//	LabelImageType::PointType origin;
	//	origin[0] = 0; 
	//	origin[1] = 0;    
	//	origin[2] = 0;    
	//	somaImage->SetOrigin( origin );
	//	LabelImageType::RegionType region;
	//	region.SetSize( im_size );
	//	region.SetIndex( start );
	//	somaImage->SetRegions( region );
	//	somaImage->Allocate();
	//	somaImage->FillBuffer(0);
	//	somaImage->Update();
	//	LabelImageType::PixelType * somaArray = somaImage->GetBufferPointer();

	//	int slice_size = im_size[1] * im_size[0];
	//	int row_size = im_size[0];
	//	std::map<unsigned short, int> classMap;
	//	for(int row=0; row<(int)featureTable[t]->GetNumberOfRows(); ++row)
	//	{
	//		classMap[featureTable[t]->GetValue(row,0).ToUnsignedShort()] = featureTable[t]->GetValueByName(row, "prediction_active_mg").ToInt();
	//	}

	//	for(int i=0; i<im_size[2]; ++i)
	//	{
	//		for(int j=0; j<im_size[1]; ++j)
	//		{
	//			for(int k=0; k<im_size[0]; ++k)
	//			{
	//				unsigned long offset = (i*slice_size)+(j*row_size)+k;
	//				if(classMap[labelArray[offset]] == 1)
	//					somaArray[offset] = labelArray[offset];
	//			}
	//		}
	//	}

	//	LabelWriterType::Pointer writer = LabelWriterType::New();
	//	writer->SetFileName(somasFilename);
	//	writer->SetInput(somaImage);
	//	writer->Update();

	//	//it = Filename.end() - 12;
	//	//Filename.erase(it, it+12);
	//	std::string centroidsFilename = Filename + "_centroids_" + time + ".txt";

	//	ofstream outFile; 
	//	outFile.open(centroidsFilename.c_str(), ios::out | ios::trunc );
	//	if ( !outFile.is_open() )
	//	{
	//		std::cerr << "Failed to Load Document: " << outFile << std::endl;
	//		return false;
	//	}
	//	//Write out the features:
	//	for(int row = 0; row < (int)featureTable[t]->GetNumberOfRows(); ++row)
	//	{
	//		if(featureTable[t]->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
	//		{
	//			outFile << featureTable[t]->GetValue(row,1).ToInt() << "\t" ;
	//			outFile << featureTable[t]->GetValue(row,2).ToInt() << "\t" ;
	//			outFile << featureTable[t]->GetValue(row,3).ToInt() << "\t" ;
	//			outFile << "\n";
	//		}
	//	}
	//	outFile.close();
	//}
	
	return true;


}

bool NucleusEditor::saveNeuronImage()
{
//    if(!labImg)
//    {
//        std::cerr << "No label image detected" << std::endl;
//        return false;
//    }
    
    if (!table)
    {
        std::cerr << "No table detected" << std::endl;
        return false;
    }
    
	std::vector< vtkSmartPointer<vtkTable> > featureTableVector;
    
    //std::cerr << "numTSlice: " << myImg->GetImageInfo()->numTSlices << std::endl;

	featureTableVector.push_back(table);
    
    
	bool prediction_active_mg_column_name_found = false;
	for(vtkIdType col = table->GetNumberOfColumns() - 1; col>=0; --col)
	{
		std::string current_column_name = table->GetColumnName(col);
		if(current_column_name.find("prediction_active_mg") != std::string::npos )
		{
			prediction_active_mg_column_name_found = true;
			break;
		}
	}
	
    if(!prediction_active_mg_column_name_found)
    {
        std::cerr << "Did not find a column named \"prediction_active_mg\"" << std::endl;
        return false;
    }
	
    QString filename = QFileDialog::getSaveFileName(this, tr("Save Soma Image As..."),lastPath, standardImageTypes);
	if(filename == "")
	{
        std::cerr << "Did not get a filename from Save Soma Image As dialog" << std::endl;
        return false;
	}
    lastPath = QFileInfo(filename).absolutePath();
	std::string Filename = filename.toStdString();
    std::string::iterator it;
	it = Filename.end() - 4;
	Filename.erase(it, it+4);
	
    typedef unsigned short LabelImagePixelType;
	typedef itk::Image<LabelImagePixelType, 3> LabelImageType;
	typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
	
    for(unsigned short t = 0; t < 1; ++t)
	{
		std::stringstream ss;
		ss << t;
		std::string time = ss.str();
		
        std::string somasFilename = Filename + "_somas_" + time + ".tif";
		
        //LabelImageType::Pointer labelImage = labImg->GetItkPtr<LabelImagePixelType>(t,0);
        typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
        
        LabelImageReaderType::Pointer reader = LabelImageReaderType::New();
        reader->SetFileName("/Users/hocheung20/Desktop/nuc_test/nuc_labels.tif");
        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject &err)
        {
            std::cerr << err << std::endl;
        }
        LabelImageType::Pointer labelImage = reader->GetOutput();
        
		LabelImageType::Pointer somaImage = LabelImageType::New();
        itk::Size<3> somaImage_size = labelImage->GetLargestPossibleRegion().GetSize();
		LabelImageType::IndexType start;    start.Fill(0);
		LabelImageType::PointType origin;   origin.Fill(0);
        somaImage->SetOrigin( origin );
		LabelImageType::RegionType region(start, somaImage_size);
		somaImage->SetRegions( region );
		somaImage->Allocate();
		somaImage->FillBuffer(0);
		somaImage->Update();
        
		std::map<unsigned short, int> classMap;
		for(vtkIdType row=0; row < featureTableVector[t]->GetNumberOfRows(); ++row)
		{
            int cell_class = featureTableVector[t]->GetValueByName(row, "prediction_active_mg").ToInt();
            //std::cerr << "Cell class: " << cell_class << std::endl;
            classMap[featureTableVector[t]->GetValue(row,0).ToUnsignedShort()] = cell_class;
		}
        
        itk::ImageRegionConstIterator<LabelImageType> labelImage_iter(labelImage, labelImage->GetLargestPossibleRegion());
        itk::ImageRegionIterator<LabelImageType> soma_image_iter(somaImage, somaImage->GetLargestPossibleRegion());
        
        while(!soma_image_iter.IsAtEnd())
        {
            LabelImagePixelType label_image_pixel_value = labelImage_iter.Get();
            if (classMap[label_image_pixel_value] == 1)
                soma_image_iter.Set(label_image_pixel_value);
            
            ++soma_image_iter;
            ++labelImage_iter;
        }
        
		LabelWriterType::Pointer writer = LabelWriterType::New();
		writer->SetFileName(somasFilename);
		writer->SetInput(somaImage);
		writer->Update();
        
		//it = Filename.end() - 12;
		//Filename.erase(it, it+12);
		std::string centroidsFilename = Filename + "_centroids_" + time + ".txt";
        
		ofstream outFile;
		outFile.open(centroidsFilename.c_str(), ios::out | ios::trunc );
		if ( !outFile.is_open() )
		{
			std::cerr << "Failed to Load Document: " << outFile << std::endl;
			return false;
		}
		
        //Write out the features:
        outFile << featureTableVector[t]->GetColumnName(1) << "\t";
        outFile << featureTableVector[t]->GetColumnName(2) << "\t";
        outFile << featureTableVector[t]->GetColumnName(3) << "\t\n";
        
        
		for(vtkIdType row = 0; row < featureTableVector[t]->GetNumberOfRows(); ++row)
		{     
            if(featureTableVector[t]->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
			{
				outFile << featureTableVector[t]->GetValue(row,1).ToInt() << "\t" ;
				outFile << featureTableVector[t]->GetValue(row,2).ToInt() << "\t" ;
				outFile << featureTableVector[t]->GetValue(row,3).ToInt() << "\t" ;
				outFile << "\n";
			}
		}
		outFile.close();
	}
	
	return true;
}
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

	//if(projectFiles.input != "" && !projectFiles.inputSaved)
	//{
	//	this->saveImage();
	//}

  QString fullPathToInput =
    QString(projectFiles.path.c_str());
  fullPathToInput += QString(projectFiles.input.c_str());
	QFileInfo inputInfo(fullPathToInput);
	if(!inputInfo.exists())
	{
		projectFiles.inputSaved = false;
	}

	if(projectFiles.input != "" && !projectFiles.inputSaved && projectFiles.type != "multi" )
	{
		this->saveImage();
	}
	else
	{
		ftk::SaveImageSeries(projectFiles.input, myImg,projectFiles.path);
		this->saveImage();

	}
	if(projectFiles.output != "" && !projectFiles.outputSaved && projectFiles.type != "multi")
	{
		this->saveResult();
	}
	else
	{
		ftk::SaveLabelSeries(projectFiles.output, labImg, projectFiles.path);
		this->saveResult();
	}
	if(projectFiles.definition != "" && !projectFiles.definitionSaved && projectFiles.type != "multi")
	{
		if( projectDefinition.Write(projectFiles.GetFullDef()) )
			projectFiles.definitionSaved = true;
	}
	//if(projectFiles.table != "" && !projectFiles.tableSaved && projectFiles.type!="multi")
	if(projectFiles.table != "" && !projectFiles.tableSaved && projectFiles.type != "multi")
	{
		this->saveTable();
	}
	else
	{
		ftk::SaveTableSeries("TableSeries.xml",nucSeg->table4DImage, projectFiles.path);
		std::string tablename = projectFiles.path+"megaTable.txt";
		ftk::SaveTable(tablename,nucSeg->megaTable);
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
	QString ext = QFileInfo(fullname).suffix();

	if(ext == "lsm")
	{
		std::cerr<<"Lsm files cannot be written. Skipping save image\n";
		return false;
	}

	bool ok;
	if(myImg->GetImageInfo()->numTSlices == 1)
	{
		if(ext == "xml")
			ok = ftk::SaveXMLImage(fullname.toStdString(), myImg);
		else
		{
			QString fullExt = "." + ext;
			fullname.remove(fullExt);
			ok = myImg->SaveChannelAs(0, fullname.toStdString(), ext.toStdString());
		}
	}
	else
	{
		ftk::Image::PtrMode mode;
		mode = static_cast<ftk::Image::PtrMode>(2);
		ftk::Image::DataType dataType = myImg->GetImageInfo()->dataType;
		unsigned char databpPix = myImg->GetImageInfo()->bytesPerPix;
		unsigned short cs = myImg->GetImageInfo()->numColumns;
		unsigned short rs = myImg->GetImageInfo()->numRows;
		unsigned short zs = myImg->GetImageInfo()->numZSlices;
		std::string name;
		std::vector< std::vector <std::string> > FileNames = myImg->GetTimeChannelFilenames();
		for(int Ch = 0; Ch  <myImg->GetImageInfo()->numChannels; Ch++)
		{
			for ( int T = 0; T <myImg->GetImageInfo()->numTSlices; T++ )
			{
				ftk::Image::Pointer ImageToSave = ftk::Image::New();
				name = ftk::SetExtension(FileNames.at(T).at(Ch),"");
				ImageToSave->AppendImageFromData3D(myImg->GetItkPtr<unsigned char>(T,Ch,mode)->GetBufferPointer(), dataType, databpPix, cs, rs, zs, name, true);
				ok = ImageToSave->SaveChannelAs(0, name,"tif");
			}
		}
	}

	projectFiles.inputSaved = ok;
	return ok;
}

bool NucleusEditor::askSaveResult()
{
	if(!labImg)
		return false;

	QString filename;

	if(labImg->GetImageInfo()->numChannels == 1)
		filename = QFileDialog::getSaveFileName(this, tr("save result as..."),lastPath, tr("tiff image (*.tif)"));
	else
		filename = QFileDialog::getSaveFileName(this, tr("save result as..."),lastPath, tr("xml image definition(*.xml)"));

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
	QString ext = QFileInfo(fullname).suffix();

	bool ok = false;
	if(labImg->GetImageInfo()->numTSlices == 1)
	{
		if(ext == "xml")
			ok = ftk::SaveXMLImage(fullname.toStdString(), labImg);
		else
		{
			QString fullExt = "." + ext;
			fullname.remove(fullExt);
			ok = labImg->SaveChannelAs(0, fullname.toStdString(), ext.toStdString());
		}
	}
	else
	{
		ftk::Image::PtrMode mode;
		mode = static_cast<ftk::Image::PtrMode>(2);
		ftk::Image::DataType dataType = labImg->GetImageInfo()->dataType;
		unsigned char databpPix = labImg->GetImageInfo()->bytesPerPix;
		unsigned short cs = labImg->GetImageInfo()->numColumns;
		unsigned short rs = labImg->GetImageInfo()->numRows;
		unsigned short zs = labImg->GetImageInfo()->numZSlices;
		std::string name;
		std::vector< std::vector <std::string> > FileNames = labImg->GetTimeChannelFilenames();
		for ( int T = 0; T <labImg->GetImageInfo()->numTSlices; T++ )
		{
			ftk::Image::Pointer ImageToSave = ftk::Image::New();
		//	name = ftk::GetFilePath(FileNames.at(t).at(0))+ftk::GetFilenameFromFullPath(FileNames.at(t).at(0);
			name = ftk::SetExtension(FileNames.at(T).at(0),"");
			ImageToSave->AppendImageFromData3D(labImg->GetItkPtr<short int>(T,0,mode)->GetBufferPointer(), dataType, databpPix, cs, rs, zs, name, true);
			ok = ImageToSave->SaveChannelAs(0, name,"tif");
		}
	}
	if(ok)
		std::cout << "done saving" << std::endl;
	else
		std::cout << "saving failed" << std::endl;
  
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


//Till we fix the issue of Loading and saving time series, we will save the Active learning results
// using this slot 

void NucleusEditor::SaveActiveLearningResults(void)
{	
	ftk::SaveTableSeries(projectFiles.GetFullTable(),nucSeg->table4DImage,projectFiles.path);
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

	//clock_t startTimer = clock();
	//multi -> Multiple images 
	if(projectFiles.input != "" && projectFiles.type!="multi")
	{
		this->loadImage( QString::fromStdString(projectFiles.GetFullInput()) );
	}
	else
	{	
		myImg = ftk::LoadImageSeries( projectFiles.GetFullInput() );
		imageNames = ftk::GetSeriesPaths( projectFiles.GetFullInput() );
		segView->SetChannelImage(myImg);
	}

	if(projectFiles.output != "" )
	{
		if(projectFiles.type!="multi")
		{
			this->loadResult( QString::fromStdString(projectFiles.GetFullOutput()) );
		}
		else
		{
			labImg = ftk::LoadImageSeriesLabels(projectFiles.GetFullOutput());
			segView->SetLabelImage(labImg, selection);
			this->updateNucSeg();
			nucSeg->SetCurrentbBox(nucSeg->bBoxMap4DImage.at(segView->GetCurrentTimeVal()));
			// this next 3 lines are temporary until we figure out how to remove compute all geometries from updateNucSeg
			//table = nucSeg->table4DImage.at(segView->GetCurrentTimeVal());
			//CreateNewTableWindow();
			//CreateNewPlotWindow();
		}
	}

	if(projectFiles.table != "")
	{
		if(projectFiles.type!="multi")
		{
			this->loadTable( QString::fromStdString(projectFiles.GetFullTable()) );
		}
		else
		{
			this->loadTableSeries( QString::fromStdString(projectFiles.GetFullTable()) );		
		}
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

	if(projectFiles.adjTables != "")
	{
		this->loadAdjTables( QString::fromStdString(projectFiles.GetFullAdjTables()) );
	}

	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "project loaded" << std::endl;
	}
	#endif
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



void NucleusEditor::loadTableSeries(QString fileName)
{
	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();

	tableVector = ftk::LoadTableSeries(fileName.toStdString());
	if(tableVector.size()==0) return;

	nucSeg->updatetable4DImage(tableVector);

	projectFiles.path = lastPath.toStdString();
	projectFiles.table = QFileInfo(fileName).fileName().toStdString();
	projectFiles.tableSaved = true;

	selection->clear();

	// Get the table corresponding to the current image
	table = tableVector.at(segView->GetCurrentTimeVal());

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
	this->setEditsEnabled(true);
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

void NucleusEditor::askload5DLabelImage()
{

	if(!myImg) return;

	if(!projectFiles.inputSaved && askSaveChanges(tr("Save changes to the input image?")) )
	{
		askSaveImage();
	}
	QStringList filesTimeList = QFileDialog::getOpenFileNames(this, "Load one or more time points", lastPath, standardImageTypes);

	const ftk::Image::Info * imInfo = myImg->GetImageInfo();
	if (imInfo->numTSlices != filesTimeList.size())	return;
	if (filesTimeList.isEmpty()) return;
	load5DLabelImage(filesTimeList);
}

void NucleusEditor::askload5DImage()
{
	if(!projectFiles.inputSaved && askSaveChanges(tr("Save changes to the input image?")) )
	{
		askSaveImage();
	}

	bool ok;
	int numChann =  QInputDialog::getInt(this, tr("Number of Channels"),
		tr("Channels:"), 1, 1, 10, 1, &ok);

	// if ok button was clicked 
	if(ok)
	{
		std::vector<QStringList> filesChannTimeList;
		for (int ch = 0; ch<numChann; ++ch)
		{
			QStringList filesTimeList = QFileDialog::getOpenFileNames(this, "Load one or more time points", lastPath, standardImageTypes);
			if (filesTimeList.isEmpty()) return;
			filesChannTimeList.push_back(filesTimeList); 
		}

		for (int ch = 0; ch<numChann-1; ++ch)
		{
			if (filesChannTimeList[ch].size()!=filesChannTimeList[ch+1].size())
				return;
		}

		load5DImage(filesChannTimeList,numChann);
	}
}



void NucleusEditor::load5DImage(std::vector<QStringList> filesChannTimeList, int numChann)
{
	// Amin: I need to put more restriction on this function.
	lastPath = QFileInfo(filesChannTimeList[0][0]).absolutePath() + QDir::separator();
	// Declare necessary variables for loading:
	myImg = ftk::Image::New();
	ftk::Image::Pointer tmp4DImage = ftk::Image::New();
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(1); //RELEASE_CONTROL mode
	std::vector<std::string> channelNames;
	std::vector<unsigned char> channelColors;
	std::vector<std::string>  filesChann;
	std::vector<std::vector<std::string> >  filesChannTime;
	int numTimes = filesChannTimeList[0].size();

	// Initialize variables:
	for (int t = 0; t< numTimes; ++t)
	{
		for (int ch = 0; ch<numChann; ++ch)
		{
			filesChann.push_back(filesChannTimeList[ch][t].toStdString());
		}
		filesChannTime.push_back(filesChann);
		filesChann.clear();
	}

	getColorInfo(numChann,&channelNames,&channelColors); //Supports upto 10 channels

	// Load Images:
	if(!myImg->LoadFilesAsMultipleChannels(filesChannTime[0],channelNames,channelColors))
		myImg = NULL;

	for(int t = 1; t <numTimes; t++)
	{
		tmp4DImage->LoadFilesAsMultipleChannels(filesChannTime[t],channelNames,channelColors);
		myImg->AppendImage(tmp4DImage,mode,true);
	}

	std::vector< std::vector <std::string> > tmp_filenames;
	for (int t = 0; t< numTimes; ++t)
	{
		std::vector <std::string> tmp_file;
		for (int ch = 0; ch<numChann; ++ch)
		{			
			std::string name = filesChannTime.at(t).at(ch);
			tmp_file.push_back(name);
		}
		tmp_filenames.push_back(tmp_file);
	}
	myImg->SetTimeChannelFilenames(tmp_filenames);

	//Display Images:
	segView->SetChannelImage(myImg);
	projectFiles.path = lastPath.toStdString();
	projectFiles.inputSaved = false;
	projectFiles.input = "SeriesFileNames.xml";
	projectFiles.type = "multi";
	//ftk::SaveImageSeries(projectFiles.input, myImg,projectFiles.path);

	load5DLabelImageAction->setEnabled(true);
	DisplayChannelsMenu();
}


//Eventually needs to move to ftkImage
void NucleusEditor::getColorInfo(int numChann,std::vector<std::string> *channelNames,std::vector<unsigned char> *channelColors)
{
	std::string channame;
	for (int ch = 0; ch<numChann; ++ch)
	{
		std::stringstream out;
		out << ch;
		std::string s = out.str();
		channame="channel"+ s;
		(*channelNames).push_back(channame);
		getColor(ch,channelColors);  
	}
}

//Eventually needs to move to ftkImage
void NucleusEditor::getColor(int numChann,std::vector<unsigned char> *channelColors)
{

	switch(numChann)
	{
		//Cyan
	case 0:
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);
		(*channelColors).push_back(255);

		break;
		//Green
	case 1:
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);

		break;
		//Red
	case 2:
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);
		(*channelColors).push_back(0);

		break;
		//Blue
	case 3:
		(*channelColors).push_back(0);
		(*channelColors).push_back(0);
		(*channelColors).push_back(255);

		break;

		//Orange 	255-165-0
	case 4:
		(*channelColors).push_back(255);
		(*channelColors).push_back(165);
		(*channelColors).push_back(0);

		break;
		//Violet 	238-130-238
	case 5:
		(*channelColors).push_back(238);
		(*channelColors).push_back(130);
		(*channelColors).push_back(238);

		break;
		//Yellow 	255-255-0
	case 6:
		(*channelColors).push_back(255);
		(*channelColors).push_back(255);
		(*channelColors).push_back(0);

		break;
		//Dark Green 	0-100-0
	case 7:
		(*channelColors).push_back(0);
		(*channelColors).push_back(100);
		(*channelColors).push_back(0);

		break;
		//Royal Blue 	65-105-225
	case 8:
		(*channelColors).push_back(65);
		(*channelColors).push_back(105);
		(*channelColors).push_back(225);

		break;
		//Gray 	190-190-190
	case 9:
		(*channelColors).push_back(190);
		(*channelColors).push_back(190);
		(*channelColors).push_back(190);

		break;
	}
}


void NucleusEditor::load5DLabelImage(QStringList filesList)
{
	lastPath = QFileInfo(filesList[0]).absolutePath() + QDir::separator();

	// Declare necessary variables for loading:
	labImg = ftk::Image::New();
	ftk::Image::Pointer tmp4DImage = ftk::Image::New();
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(1); //RELEASE_CONTROL mode
	std::vector<std::string>  filesTimeList;

	// Initialize variables:
	for (int t = 0; t< filesList.size(); ++t)
	{
		filesTimeList.push_back(filesList[t].toStdString());
	}

	if(!labImg->LoadFile(filesTimeList[0]))
		labImg = NULL;
	for(int t = 1; t <filesList.size(); t++)
	{
		tmp4DImage->LoadFile(filesTimeList[t]);
		labImg->AppendImage(tmp4DImage,mode,true);
	}
	projectFiles.output = "SeriesLabelNames.xml";
	std::vector< std::vector <std::string> > tmp_filenames;
	for(int i = 0; i< filesTimeList.size(); ++i)
	{
		std::vector <std::string> tmp_file;
		tmp_file.push_back(filesTimeList.at(i));
		tmp_filenames.push_back(tmp_file);
	}
	labImg->SetTimeChannelFilenames(tmp_filenames);

	//Display Images Tables and Plots:
	segView->SetLabelImage(labImg, selection);
	this->updateNucSeg();  
	table = nucSeg->table4DImage.at(segView->GetCurrentTimeVal());
	CreateNewTableWindow();
	CreateNewPlotWindow();
	projectFiles.path = lastPath.toStdString();

	projectFiles.tableSaved = false;
	projectFiles.inputSaved = false;
	projectFiles.outputSaved = false;
	projectFiles.adjTablesSaved = true;
	projectFiles.table = "TableSeries.xml";
	projectFiles.type = "multi";


	//For storing edit information.. need to have a bBoxMap set 
	nucSeg->SetCurrentbBox(nucSeg->bBoxMap4DImage.at(segView->GetCurrentTimeVal()));

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
	this->setEditsEnabled(false);
}


//void NucleusEditor::loadImageSeries(QString fileName)
//{
//
//	lastPath = QFileInfo(fileName).absolutePath() + QDir::separator();
//	QString name = QFileInfo(fileName).fileName();
//	QString myExt = QFileInfo(fileName).suffix();
//	if(myExt == "xml")
//	{
//		myImg = ftk::LoadXMLImage(fileName.toStdString());
//	}
//	else
//	{
//		myImg = ftk::Image::New();
//		if(!myImg->LoadFile(fileName.toStdString()))
//			myImg = NULL;
//	}
//	
//	//this->updateNucSeg();
//	if(nucSeg)
//	{
//		segView->SetCenterMapPointer( 0 );
//		segView->SetBoundingBoxMapPointer( 0 );
//		delete nucSeg;
//		nucSeg = NULL;
//	}
//
//	segView->SetChannelImage(myImg);
//	projectFiles.path = lastPath.toStdString();
//	projectFiles.input = name.toStdString();
//	projectFiles.inputSaved = true;
//}

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

void NucleusEditor::getCentroids(void)
{

	ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
	iCalc->SetInputImages(myImg, labImg, 0, 0);
	iCalc->GetObjectCentroids(table);								//Create a new table
	delete iCalc;

	updateViews();

}

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
	if(table){
		if(tableVector.size()==0){
			sqlite3 *dbConn;
			dbConn = ftk::sqliteOpenConnection();
			if( dbConn ){
				std::vector<std::string> col_names;
				std::cout<<"There is one table loaded\n";
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
			else
				std::cout<<"Error opening database file. Please check if the database file is in the present working dir\n";
		}
		else{
			std::cout<<"There are "<<tableVector.size()<<" tables loaded\n";
			std::string table_name = "IMAGE_TEST";
			for (int i = 0; i< tableVector.size(); ++i){
				sqlite3 *dbConn;
				dbConn = ftk::sqliteOpenConnection();
				if( dbConn ){
					std::vector<std::string> col_names;
					for (int col = 1; col< tableVector.at(i)->GetNumberOfColumns(); ++col){
						std::string temp3=tableVector.at(i)->GetColumnName(col);
						col_names.push_back(temp3);
					}
					ftk::checkForUpdate( dbConn, col_names );
					std::string image_name = imageNames.at(i);
					char *im_nm_cstr = new char [image_name.size()+1];
					strcpy (im_nm_cstr, image_name.c_str());
					char *path_nm_cstr = new char [lastPath.toStdString().size()+1];
					strcpy (path_nm_cstr, ftk::GetFilePath(image_name).c_str());
					std::vector< double > table_array;
					for (int row = 0; row < tableVector.at(i)->GetNumberOfRows(); ++row)
						for (int col = 0; col< tableVector.at(i)->GetNumberOfColumns(); ++col)
							table_array.push_back(tableVector.at(i)->GetValue(row,col).ToDouble());
					int sql_db_img_id = ftk::GenericInsert( dbConn, im_nm_cstr, table_name.c_str(), path_nm_cstr, table_array,tableVector.at(i)->GetNumberOfColumns(), tableVector.at(i)->GetNumberOfRows(), col_names );
					std::cout << "The image ID on the database is: " << sql_db_img_id << std::endl;
					ftk::sqliteCloseConnection( dbConn );
				}
				else{
					std::cout<<"Error opening database file. Please check if the database file is in the present working dir\n";
					break;
				}
			}
		}
	}
}




void NucleusEditor::startActiveLearningwithFeat()
{	
	classify_from_model = false;
	std::vector< vtkSmartPointer<vtkTable> > VectorOfTables;
	vtkSmartPointer<vtkTable> featureTable;
	if(myImg->GetImageInfo()->numTSlices==1)
	{
		VectorOfTables.push_back(table);	
		featureTable = table;
	}
	else
	{
		VectorOfTables = nucSeg->table4DImage;	
		featureTable = nucSeg->megaTable;
	}
	AL->SetTablesToClassify(VectorOfTables);
	AL->SetTableForTraining(featureTable);
	AL->SetLabelView(segView);
	AL->RunALClassification(classify_from_model);

}	


void NucleusEditor::SaveActiveLearningModel()
{
	if(active_model.size() == 0) return;

	//Save Dialog
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Training Model As..."),lastPath, tr("TEXT(*.txt)"));
	if(filename == "")
		return;
	lastPath = QFileInfo(filename).absolutePath();	

	std::string Filename = filename.toStdString();;
	ofstream outFile; 
	outFile.open(Filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return;
	}
	//Write out the average distance:
	for(int i=0; i < (int)active_model.size(); ++i)
	{
		outFile << active_model[i].first << "\t";		
	}
	outFile << "\n";
	for(int j=0; j<(int)active_model[0].second.size(); ++j)
	{
		for(int i=0; i < (int)active_model.size(); ++i)
		{
			outFile << active_model[i].second.get(j) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
}


void NucleusEditor::classifyFromActiveLearningModel()
{
	classify_from_model = true;
	std::vector< vtkSmartPointer<vtkTable> > VectorOfTables;
	vtkSmartPointer<vtkTable> featureTable;
	if(myImg->GetImageInfo()->numTSlices==1)
	{
		VectorOfTables.push_back(table);
		featureTable = table;
	}
	else
	{
		VectorOfTables = nucSeg->table4DImage;
		featureTable = nucSeg->megaTable;
	}
	AL->SetTablesToClassify(VectorOfTables);
	AL->SetTableForTraining(featureTable);
	AL->SetLabelView(segView);
	AL->RunALClassification(classify_from_model);
}


void NucleusEditor::ExtractClassificationResult()
{
	if(myImg->GetImageInfo()->numTSlices==1)
	{
		table = AL->GetClassificationResult()[0];
	}
	else
	{		
		nucSeg->table4DImage = AL->GetClassificationResult();
	}

	if(!classify_from_model)
	{
		active_model = AL->GetALModel();
	}

	prediction_names = ftk::GetColumsWithString( "prediction_active" , table);
	selection->clear();	
	projectFiles.tableSaved = false;
	activeRun = 1;		
	this->updateViews();
}


/****************************************************************
Active Validation
*****************************************************************/
void NucleusEditor::startActiveValidation()
{
	int numbin;
	double delta;

	//Set parameters needed in validation
	SamplePercentDialog *dialog = new SamplePercentDialog(this);
	if( dialog->exec() )
	{
		delta = dialog->getDelta();
		numbin = dialog->getNumbin();		
	}
	else
		return;
	delete dialog;

	a_v = new ActiveValidation();
	a_v->Initializing(this->table, numbin, delta);
	a_v->StratifingwithoutLabel();
	a_v->CalculatenumToBeDrawn();
	queries.clear();
	queries = a_v->RandomSamplingwithoutLabel(a_v->stratas[0].numtobeDrawn, a_v->stratas[0].numleft,0);
	learned.clear();
	PreLabellearning(0, numbin);
}

void NucleusEditor::ActiveValidationQuery()
{
	a_v->CalculatehkwithoutLabel(a_v->current, learned);
	a_v->Calculatephatak(a_v->current);
	a_v->CalculateVarphatk(a_v->current);
	a_v->CalculateSigmak(a_v->current);
	a_v->UpdateBins(a_v->current);

	a_v->current += 1;
	
	if(a_v->current < a_v->numbin)
	{
		queries.clear();
		queries = a_v->RandomSamplingwithoutLabel(a_v->stratas[a_v->current].numtobeDrawn,
				  a_v->stratas[a_v->current].numleft, a_v->current);
		learned.clear();
		PreLabellearning( a_v->current,  a_v->numbin);
	}
	else
	{
		bool converge = false;
		a_v->Calculatephata();
		a_v->CalculateVarphat();
		converge = a_v->ConvergeCheck();

		if(converge)
		{
			a_v->t++;
		}
		else
		{
			a_v->t = 0;
		}

		a_v->numiteration++;
		a_v->current = 0;
		std::cout<<"iteration "<<a_v->numiteration<<std::endl;
		std::cout<<"phat = "<<a_v->phat<<std::endl;
		std::cout<<"varphat = "<<sqrt( (double)a_v->varphat)<<std::endl;

		if(a_v->t < 2)
		{
			a_v->CalculatenumToBeDrawn();
			queries.clear();
			queries = a_v->RandomSamplingwithoutLabel(a_v->stratas[a_v->current].numtobeDrawn,
				      a_v->stratas[a_v->current].numleft, a_v->current);
			learned.clear();
		    PreLabellearning(a_v->current,  a_v->numbin);
		}
		else
		{
			std::cout<<"................................................."<<std::endl;
			std::cout<<"active validation done !"<<std::endl;
			std::cout<<"number of iteration "<<a_v->numiteration<<std::endl;
			std::cout<<"number of sample "<<a_v->numsampled<<std::endl;
			std::cout<<"phat = "<<a_v->phat<<std::endl;
			std::cout<<"varphat = "<<sqrt( (double)a_v->varphat)<<std::endl;
			delete a_v;
		}
	}	
}

void NucleusEditor::Labellearning(std::vector<std::pair<int,int> > query)
{
	for(int i = 0; i < query.size(); i++)
	{
		learned.push_back(query[i].second);
	}

	if( looptime < queries.size())
	{
		std::vector<QImage> subsnapshots;
		std::vector<int > subqueries;
		int i = 0;
		for(i,looptime; (i < 5)&& (looptime < queries.size()); looptime++,i++)
		{	
			std::cout<<queries[looptime]<<std::endl;
			subsnapshots.push_back( segView->getSnapshotforID(queries[looptime] + 1) );
			subqueries.push_back(queries[looptime]);
		}	
		ActiveLearningDialog *alDialog = new ActiveLearningDialog(subsnapshots, table, a_v->current,
											subqueries, a_v->numbin);	

		connect(alDialog, SIGNAL(finishquery(std::vector<std::pair<int,int> >)), 
				this, SLOT(Labellearning (std::vector<std::pair<int,int> >)));
		connect( alDialog, SIGNAL( next(std::vector<std::pair<int,int> >) ),
				this, SLOT(Labellearning(std::vector<std::pair<int,int> >) ) );
		alDialog->show();
	}
	else
	{
		ActiveValidationQuery();
	}	
}

void NucleusEditor::PreLabellearning(int classval, int numclass)
{
	std::vector<QImage> subsnapshots;
	std::vector<int > subqueries;
	std::cout<<"queries size is "<< queries.size() <<std::endl;
	this->looptime = 0;
	int i = 0;
	for( i,looptime; (i < 5) && (looptime < queries.size()); looptime++,i++)
	{	std::cout<<queries[looptime]<<std::endl;
		subqueries.push_back(queries[looptime]);
		subsnapshots.push_back( segView->getSnapshotforID(queries[looptime] + 1) );
	}	
	ActiveLearningDialog *alDialog = new ActiveLearningDialog(subsnapshots, table, a_v->current,
									subqueries, a_v->numbin);	

	connect(alDialog, SIGNAL(finishquery(std::vector<std::pair<int,int> >)), 
			this, SLOT(Labellearning (std::vector<std::pair<int,int> >)));
	connect( alDialog, SIGNAL( next(std::vector<std::pair<int,int> >) ),
			this, SLOT( Labellearning(std::vector<std::pair<int,int> >) ) );
	alDialog->show();
}
/****************************************************************
function for heatmap of active learning result
*****************************************************************/
void NucleusEditor::HeatmapforActivelearning( vtkSmartPointer<vtkTable> table, int class_num)
{
	QMessageBox msgBox;
	//msgBox.setText("The document has been modified.");
	msgBox.setInformativeText("Do you want a heatmap for the active learning result ?");
	msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
	msgBox.setDefaultButton(QMessageBox::Cancel);
	int ret = msgBox.exec();

	switch (ret) 
	{
		case QMessageBox::Ok:
			this->ClusteringWithinLabels(table,class_num);
			break;
		case QMessageBox::Cancel:
			//nothing is done 
			break;
		default:
			// should never be reached
			break;
	}
}
void NucleusEditor::ClusteringWithinLabels( vtkSmartPointer<vtkTable> table, int num_class)
{
	std::cout<<"table row number "<<table->GetNumberOfColumns()<<std::endl;
	std::cout<<"table row number "<<table->GetNumberOfRows()<<std::endl;
	int* order = new int[table->GetNumberOfRows()];
	std::vector<std::vector<std::vector<double > > > treesdata;

	vtkSmartPointer<vtkTable> rearrangedtable = vtkSmartPointer<vtkTable>::New();
	rearrangedtable->Initialize();
	for(int col=0; col<table->GetNumberOfColumns(); ++col)
	{
		vtkSmartPointer<vtkDoubleArray> columns = vtkSmartPointer<vtkDoubleArray>::New();
		columns->SetName(table->GetColumnName(col));
		rearrangedtable->AddColumn(columns);	
	}

	std::vector< vtkSmartPointer<vtkTable> > tables;
	for(int i=0; i<=num_class; i++)
	{
		vtkSmartPointer<vtkTable> temptable = vtkSmartPointer<vtkTable>::New();
		temptable->Initialize();
		for(int col=0; col<table->GetNumberOfColumns(); ++col)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName(table->GetColumnName(col));
			temptable->AddColumn(column);	
		}
		tables.push_back(temptable);
	}
	for(int i=0; i<table->GetNumberOfRows(); i++)
	{	
		std::cout<<table->GetColumnByName("prediction_active")->GetSize()<<std::endl;
		vtkVariant cn = table->GetValueByName(i, "prediction_active" );
		std::cout<<cn<<std::endl;
		tables[cn.ToInt()]->InsertNextRow(table->GetRow(i));
	}
	for(int i=0; i<=num_class; i++)
	{
		for(int j=0; j<tables[i]->GetNumberOfRows(); j++)
		{
			rearrangedtable->InsertNextRow(tables[i]->GetRow(j));
		}
	}

	int count = 0;// for counting of rows already in the order
	int tablecount = 0;// for counting of tables already in the order
	for(int i=0; i<=num_class; i++)
	{
		std::cout<<"Num_class: " << i <<" of "<< num_class <<std::endl;
		if(tables[i]->GetNumberOfRows()==0)
		{
			std::cout<<"table is empty for this class..."<<i<<std::endl;
			continue;
		}
		pWizard = new PatternAnalysisWizard( tables[i], PatternAnalysisWizard::_CLUS, "", "", this);
		connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
		connect(pWizard, SIGNAL(enableModels()), this, SLOT(EnableModels()));
		pWizard->setWindowTitle(QString::number(i));
		pWizard->exec();
		
		if(pWizard->result())
		{
			std::cout<<"show heatmap";
			vtkSmartPointer<vtkTable> featureTable;
			featureTable = pWizard->getExtractedTable();

			//Heatmap* HeatmapWin = new Heatmap();
			//HeatmapWin->setModels(featureTable, this->selection);
			//HeatmapWin->runClus();
			//HeatmapWin->showGraph();
			////////////////////////////////////////////////////////////////////////
			double** datas;
			datas = new double*[featureTable->GetNumberOfRows()];

			for (int i = 0; i < featureTable->GetNumberOfRows(); i++)
			{
				datas[i] = new double[featureTable->GetNumberOfColumns() - 1 + 2 ];
			}

			for(int i = 0; i < featureTable->GetNumberOfRows(); i++)
			{		
				for(int j = 1; j < featureTable->GetNumberOfColumns(); j++)
				{
					vtkVariant temp = featureTable->GetValue(i, j);
					datas[i][j-1] = temp.ToDouble();
				}
			}

			clusclus* cc = new clusclus(datas, (int)featureTable->GetNumberOfRows(), (int)featureTable->GetNumberOfColumns() - 1);
			cc->RunClusClus();
			for(int i = 0; i < featureTable->GetNumberOfRows(); i++)
			{
				order[count++] = cc->optimalleaforder[i] + tablecount;
			}
			tablecount += featureTable->GetNumberOfRows();

			std::vector<std::vector<double > > treedata;
			for(int i = 0; i<featureTable->GetNumberOfRows() -1; i++)
			{	
				std::vector<double > temp;
				for(int j = 0; j<4; j++)
				{
					temp.push_back (cc->treedata[i][j]);
				}
				treedata.push_back(temp);
			}
			treesdata.push_back(treedata);
			delete datas;
			///////////////////////////////////////////////////////////////////////
		}	
	}
	Heatmap* HeatmapWin = new Heatmap();
	HeatmapWin->setModels(rearrangedtable, this->selection);
	HeatmapWin->setOrders(order);
	HeatmapWin->setMultipleTreeData(treesdata);
	HeatmapWin->creatDataForHeatmap(0.2);
	HeatmapWin->showGraph();
}






//**********************************************************************
// SLOT: start the training dialog:
//**********************************************************************
void NucleusEditor::startTraining()
{
	if(!table) return;

	TrainingDialog *d = new TrainingDialog(table, "train","", table->GetNumberOfRows(),this);
	connect(d, SIGNAL(changedTable()), this, SLOT(updateViews()));
	d->show();
}

//**********************************************************************
// SLOT: start tracking:
//**********************************************************************
#ifdef USE_TRACKING
void NucleusEditor::startTracking()
{
	if(!labImg) return;
	if(!myImg) return;
	if(myImg->GetImageInfo()->numTSlices <3) return;


		//Get Channels in current Image:
	QVector<QString> chs = getChannelStrings();
	int nucChannel = 0;
	ParamsFileDialog *dialog = new ParamsFileDialog(lastPath,getChannelStrings(),this);
	if( dialog->exec() )
	{
		nucChannel = dialog->getChannelNumber();
	}
	delete dialog;

	TrackingDialog * trackdialog= new  TrackingDialog();
	trackdialog->exec();
	if(trackdialog->result())
	{
		this->closeViews();
		mfcellTracker = new MultiFrameCellTracker();
		mfcellTracker->setTrackParameters(trackdialog->getParameters());
		mfcellTracker->settrackresultFolders(trackdialog->getFolders());
		mfcellTracker->setChannelToTrack(nucChannel);
		mfcellTracker->setTrackImages(myImg,labImg);
		labImg = mfcellTracker->getTrackImages();

		//Display Images Tables and Plots:
		segView->SetLabelImage(labImg, selection);
		this->updateNucSeg();  
		table = nucSeg->table4DImage.at(segView->GetCurrentTimeVal());
		CreateNewTableWindow();
		CreateNewPlotWindow();
		nucSeg->SetCurrentbBox(nucSeg->bBoxMap4DImage.at(segView->GetCurrentTimeVal()));
	}
	delete trackdialog;



}
//**********************************************************************
// SLOT: start kymograph view:
//**********************************************************************
void NucleusEditor::displayKymoGraph()
{
	if(!labImg) return;
	if(!myImg) return;
	kymoView = new TrackingKymoView(myImg,nucSeg->featureVector4DImage,segView,selection);
	myview = new Image3DView(myImg,labImg,segView,selection);	// Constructor

}
#endif
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
// SLOT: Clustering Heatmap
//**********************************************************************
//void NucleusEditor::runClus()
//{
//	#ifdef USE_Clusclus
//	this->HeatmapWin = new Heatmap();
//	if( table->GetNumberOfRows() <= 0)
//	{
//		QMessageBox mes;
//		mes.setText("Please compute cell features first!");
//		mes.exec();
//	}
//	else
//	{
//		
//
//		pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_CLUS, "", "", this);
//		connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
//		connect(pWizard, SIGNAL(enableModels()), this, SLOT(EnableModels()));
//		pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
//		pWizard->exec();
//		
//		if(pWizard->result())
//		{
//			vtkSmartPointer<vtkTable> featureTable;
//			featureTable = pWizard->getExtractedTable();
//
//			for(int i = 0; i < 5; i++)
//			{
//				for(int j = 0; j<featureTable->GetNumberOfColumns(); j++)
//				{
//					cout<<featureTable->GetValue(i, j)<<"\t";
//				}
//				cout<<endl;
//			}
//
//			this->HeatmapWin->setModels(featureTable, this->selection);
//			this->HeatmapWin->runClus();
//			this->HeatmapWin->showGraph();
//		}
//	}
//	#endif
//}

void NucleusEditor::runClus()
{
	//#ifdef USE_Clusclus
	this->biheatmap = new BiHeatmap();
	if( table->GetNumberOfRows() <= 0)
	{
		QMessageBox mes;
		mes.setText("Please compute cell features first!");
		mes.exec();
	}
	else
	{
		

		pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_CLUS, "", "", this);
		connect(pWizard, SIGNAL(changedTable()), this, SLOT(updateViews()));
		connect(pWizard, SIGNAL(enableModels()), this, SLOT(EnableModels()));
		pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
		pWizard->exec();
		
		if(pWizard->result())
		{
			vtkSmartPointer<vtkTable> featureTable;
			featureTable = pWizard->getExtractedTable();

			std::vector<std::vector<double > > points;
			points.resize(featureTable->GetNumberOfRows());
			for(int i = 0; i < featureTable->GetNumberOfRows(); i++)
			{
				for(int j = 1; j < featureTable->GetNumberOfColumns(); j++)
				{
					points[i].push_back(featureTable->GetValue(i,j).ToDouble());
				}
			}
			std::cout<<".............aha1"<<std::endl;
			Bicluster* bicluster = new Bicluster();
			bicluster->setDataToBicluster(points);
			bicluster->biclustering();
			bicluster->WriteFile("order1.txt", "order2.txt");
			
			std::cout<<".............aha2"<<std::endl;
			this->biheatmap->setModels(featureTable, selection);
			this->biheatmap->setDataForHeatmap(bicluster->order1, bicluster->order2);
			this->biheatmap->setDataForTree1(bicluster->levels1);
			this->biheatmap->setDataForTree2(bicluster->levels2);
			this->biheatmap->showHeatmap();
			this->biheatmap->showTree1();
			this->biheatmap->showTree2();

			delete bicluster;
		}
	}
	//#endif
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
	bool k_mutual;

	QVector<QString> classes;
	int max_class = 0;
	for(int col=0; col<(int)table->GetNumberOfColumns(); ++col)
	{	
		std::string current_column = table->GetColumnName(col);
		if(current_column.find("prediction") != std::string::npos )
		{
			for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
			{
				if(table->GetValue(row,col).ToInt() > max_class)
					max_class = table->GetValue(row,col).ToInt();
			}
			break;
		}
	}

	for(int i=0; i<max_class; ++i)
		classes.push_back(QString::number(i+1));
	

	QueryDialog *dialog = new QueryDialog(1, classes, false, this);
	if( dialog->exec() )
	{
		IDs = dialog->parseIDs();
		k = dialog->parseK();
		if(IDs.at(0) == 0)
			Class_src = dialog->getSourceClass();
		Class_dest = dialog->getDestClass();
		k_mutual = dialog->getKMutual();
	}
	delete dialog;

	//std::map<int, ftk::Object::Point> *	centerMap;
	//centerMap = segView->GetCenterMapPointer();
	//std::map<int, ftk::Object::Point>::iterator it;
	std::map< unsigned int, std::vector<double> > centroidMap;
	//for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	for ( int row=0; row<(int)table->GetNumberOfRows(); ++row )
	{
		//unsigned int id = (unsigned int)(*it).first;
		//std::vector<double> c;
		//c.push_back((*it).second.x);
		//c.push_back((*it).second.y);
		//c.push_back((*it).second.z);
		unsigned int id = table->GetValue(row,0).ToUnsignedInt();
		std::vector<double> c;
		c.push_back(table->GetValue(row,1).ToDouble());
		c.push_back(table->GetValue(row,2).ToDouble());
		c.push_back(table->GetValue(row,3).ToDouble());
		centroidMap[id] = c;
	}

	kNearestObjects<3>* KNObj = new kNearestObjects<3>(centroidMap);
	KNObj->setFeatureTable(table);
	std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	if(IDs.at(0) == 0)
		kNeighborIDs = KNObj->k_nearest_neighbors_All(k, Class_dest, Class_src);
	else
		kNeighborIDs = KNObj->k_nearest_neighbors_IDs(IDs, k, Class_dest);

	std::string full_string;
	std::stringstream ss1;
	ss1 << k;
	if(Class_dest == 0)
	{
		full_string = "D_k=" + ss1.str() + "_class=all_" ;
	}
	else
	{
		std::stringstream ss2;
		ss2 << Class_dest;
		full_string = "D_k=" + ss1.str() + "_class=" + ss2.str() + "_";
	}
	table->RemoveColumnByName(full_string.c_str());
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName(full_string.c_str());
	column->SetNumberOfValues((int)table->GetNumberOfRows());
	table->AddColumn(column);
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		table->SetValueByName(row, full_string.c_str(), 0);
	}
	for(int i=0; i < (int)kNeighborIDs.size(); ++i)
	{
		int Id = kNeighborIDs.at(i).at(0).first;
		double avg_dist = average(kNeighborIDs.at(i));
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			if(table->GetValue(row,0).ToInt() == Id)
			{
				table->SetValueByName(row, full_string.c_str(), vtkVariant(avg_dist));
				break;
			}
		}
	}
	this->updateViews();

	vtkSmartPointer<vtkTable> kNeighborTable = KNObj->vectorsToGraphTable(kNeighborIDs);
	segView->SetKNeighborTable(kNeighborTable);
	segView->SetKNeighborsVisibleOn(k_mutual);		

	#ifdef USE_QT_TESTING
	std::cout << kNeighborTable->GetNumberOfRows() << " neighbors potentially connected" << std::endl;
	#endif

	delete KNObj;
	KNObj = NULL;
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
	for(int col=0; col<(int)table->GetNumberOfColumns(); ++col)
	{	
		std::string current_column = table->GetColumnName(col);
		if(current_column.find("prediction") != std::string::npos )
		{
			for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
			{
				if(table->GetValue(row,col).ToInt() > max_class)
					max_class = table->GetValue(row,col).ToInt();
			}
			break;
		}
	}

	for(int i=0; i<max_class; ++i)
		classes.push_back(QString::number(i+1));

	QueryDialog *dialog = new QueryDialog(2, classes, false, this);
	if( dialog->exec() )
	{
		IDs = dialog->parseIDs();
		radius = dialog->parseRad();
		if(IDs.at(0) == 0)
			Class_src = dialog->getSourceClass();
		Class_dest = dialog->getDestClass();
	}
	delete dialog;

	//std::map<int, ftk::Object::Point> *	centerMap;
	//centerMap = segView->GetCenterMapPointer();
	//std::map<int, ftk::Object::Point>::iterator it;
	std::map< unsigned int, std::vector<double> > centroidMap;
	//for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	for ( int row=0; row<(int)table->GetNumberOfRows(); ++row )
	{
		//unsigned int id = (unsigned int)(*it).first;
		//std::vector<double> c;
		//c.push_back((*it).second.x);
		//c.push_back((*it).second.y);
		//c.push_back((*it).second.z);
		unsigned int id = table->GetValue(row,0).ToUnsignedInt();
		std::vector<double> c;
		c.push_back(table->GetValue(row,1).ToDouble());
		c.push_back(table->GetValue(row,2).ToDouble());
		c.push_back(table->GetValue(row,3).ToDouble());
		centroidMap[id] = c;
	}

	kNearestObjects<3>* KNObj = new kNearestObjects<3>(centroidMap);
	KNObj->setFeatureTable(table);
	std::vector<std::vector< std::pair<unsigned int, double> > > radNeighborIDs;
	if(IDs.at(0) == 0)
		radNeighborIDs = KNObj->neighborsWithinRadius_All(radius, Class_dest, Class_src);
	else
		radNeighborIDs = KNObj->neighborsWithinRadius_IDs(IDs, radius, Class_dest);

	std::cout << radNeighborIDs.size() << "\n";

	std::string full_string;
	std::stringstream ss1;
	ss1 << radius;
	if(Class_dest == 0)
	{
		full_string = "D_rad=" + ss1.str() + "_class=all_";
	}
	else
	{
		std::stringstream ss2;
		ss2 << Class_dest;
		full_string = "D_rad=" + ss1.str() + "_class=" + ss2.str() + "_";
	}
	table->RemoveColumnByName(full_string.c_str());
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName(full_string.c_str());
	column->SetNumberOfValues((int)table->GetNumberOfRows());
	table->AddColumn(column);
	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		table->SetValueByName(row, full_string.c_str(), 0);
	}
	for(int i=0; i < (int)radNeighborIDs.size(); ++i)
	{
		int Id = radNeighborIDs.at(i).at(0).first;
		//double avg_dist = average(radNeighborIDs.at(i));
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			if(table->GetValue(row,0).ToInt() == Id)
			{
				table->SetValueByName(row, full_string.c_str(), vtkVariant(radNeighborIDs[i].size()-1));
				break;
			}
		}
	}
	this->updateViews();

	vtkSmartPointer<vtkTable> radNeighborTable = KNObj->vectorsToGraphTable(radNeighborIDs);	
	segView->SetRadNeighborTable(radNeighborTable);
	segView->SetRadNeighborsVisibleOn();

	#ifdef USE_QT_TESTING
	std::cout << radNeighborTable->GetNumberOfRows() << " neighbors in radius" << std::endl;
	#endif
	delete KNObj;

}


//**********************************************************************
// SLOT: Compute Diffusion Map
//**********************************************************************
void NucleusEditor::computeDistanceMap()
{
	if(!table) return;

	pWizard = new PatternAnalysisWizard( table, PatternAnalysisWizard::_ACTIVE,"","", this);
	pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
	pWizard->exec();

	if(pWizard->result())
	{	
		// Extracted Table containing features of trained model 
		pawTable = pWizard->getExtractedTable();
		//if(diffusion_map)
		//	delete diffusion_map;
		
		diffusion_map = new DiffusionMap();
		diffusion_map->Initialize(table, true);
		diffusion_map->ComputeDiffusionMap();		
	}
}

//**********************************************************************
// SLOT: Compute Diffusion Map
//**********************************************************************
void NucleusEditor::queryKDiffusion()
{
	if(!table) return;
	if(!diffusion_map)
	{
		std::cout<< "Diffusion Map not computed !!\n";
		return;
	}

	std::vector<unsigned int> IDs;
	unsigned int k;
	bool k_mutual;

	QVector<QString> classes;
	classes.clear();
	QueryDialog *dialog = new QueryDialog(1, classes, true, this);
	if( dialog->exec() )
	{
		IDs = dialog->parseIDs();
		k = dialog->parseK();
		k_mutual = dialog->getKMutual();
	}
	delete dialog;
	std::cout << k;
	vtkSmartPointer<vtkTable> kNeighborTable = diffusion_map->GetKDiffusionNeighbors(IDs, k);
	segView->SetKNeighborTable(kNeighborTable);
	segView->SetKNeighborsVisibleOn(k_mutual);
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
			delete *table_it;
			tblWin.erase(table_it);
			return;
		}
	}

	std::vector<PlotWindow *>::iterator plot_it;
	for ( plot_it = pltWin.begin(); plot_it < pltWin.end(); plot_it++ )
	{
		if( *plot_it == view )
		{
			delete *plot_it;
			pltWin.erase(plot_it);
			return;
		}
	}

	std::vector<HistoWindow *>::iterator hist_it;
	for ( hist_it = hisWin.begin(); hist_it < hisWin.end(); hist_it++ )
	{
		if( *hist_it == view )
		{
			delete *hist_it;
			hisWin.erase(hist_it);
			return;
		}
	}

	std::vector<FTKRenderWindow *>::iterator rend_it;
	for ( rend_it = renWin.begin(); rend_it < renWin.end(); rend_it++ )
	{
		if( *rend_it == view )
		{
			delete *rend_it;
			renWin.erase(rend_it);
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

	for(int p=0; p<(int)hisWin.size(); ++p)
		hisWin.at(p)->close();
}

//Call this slot when the table has been modified (new rows or columns) to update the views:
void NucleusEditor::updateViews()
{
	//Show colored seeds after kPLS has run
	if( kplsRun )
	{
		segView->SetClassMap(table, prediction_names);
		showCentroidsAction->setChecked(true);
		segView->SetCentroidsVisible(true);
		kplsRun = 0;
		#ifdef USE_QT_TESTING
		if(this->TestInputFile != "")
		{
			std::cout << "KPLS classification complete" << std::endl;
		}
		#endif
	}

	//Show colored seeds after Active Learning Classification has finished
	//if( activeRun && prediction_names.size()>0 )
	if( prediction_names.size()>0 )
	{
		segView->SetClassMap(table, prediction_names);
		showCentroidsAction->setChecked(true);
		segView->SetCentroidsVisible(true);
	}

	segView->update();

	for(int p=0; p<(int)tblWin.size(); ++p)
		tblWin.at(p)->update();

	for(int p=0; p<(int)pltWin.size(); ++p)
		pltWin.at(p)->update();

	for(int p=0; p<(int)hisWin.size(); ++p)
		hisWin.at(p)->update();

	for(int p=0; p<(int)renWin.size(); ++p)
		renWin.at(p)->update();

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


//*******************************************************************************
// Create new Render Window
//*******************************************************************************
void NucleusEditor::CreateNewRenderWindow(void)
{
	if(!table) return;

	renWin.push_back(new FTKRenderWindow());
	connect(renWin.back(), SIGNAL(closing(QWidget *)), this, SLOT(viewClosing(QWidget *)));
	renWin.back()->show();
	renWin.back()->setModels(table,selection);
	//renWin.back()->show();
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
	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "toggleIDs" << std::endl;
	}
	#endif
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
		action->setObjectName(tr(channel_names.at(i).c_str()));
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
		if( !projectDefinition.inputs.size() )
			projectDefinition.MakeDefaultNucleusSegmentation(nucChannel);
	}

	nucSeg = new ftk::NuclearSegmentation();
	nucSeg->SetInput(myImg,"nuc_img", nucChannel);
	nucSeg->SetLabelImage(labImg,"lab_img");
	const ftk::Image::Info * imInfo = labImg->GetImageInfo();

#ifdef USE_TRACKING
	if ((imInfo->numTSlices > 1)&& mfcellTracker)
		nucSeg->SetTrackFeatures(mfcellTracker->getTrackFeatures());// needs to be called before ComputeAllGeometries.
#endif

	if (imInfo->numTSlices > 1)
	{
		nucSeg->ComputeAllGeometries(imInfo->numTSlices);
		nucSeg->createMegaTable();
		nucSeg->AddTimeToMegaTable();
		segView->SetCenterMapVectorPointer(nucSeg->centerMap4DImage);
		segView->SetBoundingBoxMapVectorPointer(nucSeg->bBoxMap4DImage);
	}
	else
	{
		nucSeg->ComputeAllGeometries();
	}
	segView->SetCenterMapPointer( nucSeg->GetCenterMapPointer() );
	segView->SetBoundingBoxMapPointer( nucSeg->GetBoundingBoxMapPointer() );
}

void NucleusEditor::startEditing(void)
{
	if(!myImg || !labImg)
		return;

	std::string log_entry = "NUCLEAR_SEGMENTATION , ";
	log_entry += ftk::NumToString(nucSeg->GetNumberOfObjects()) + " , ";
	log_entry += ftk::TimeStamp();
	ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);

	projectFiles.nucSegValidated = false;
	//if(projectFiles.type == "multi")
	//	setEditsEnabled(false);
	//else
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

	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	std::vector<int> allIds;

//	Get the new class number:
	prediction_names = ftk::GetColumsWithString( "prediction" , table);
	QVector<QString> classifiers;
	for(int i=0; i<(int)prediction_names.size(); ++i)
	{
		classifiers.push_back(QString::fromStdString(prediction_names[i]));
	}

	std::string classifier;
	unsigned short new_class;
	ChangeClassDialog *dialog = new ChangeClassDialog (classifiers,this);
	if( dialog->exec() )
	{
		classifier = dialog->getClassColumn();
		new_class = dialog->getClass();		
		std::cout << classifier << "_" << new_class << "\n";
	}
	else
		return;
	delete dialog;

	vtkAbstractArray * columnData = table->GetColumn(0);
	for (vtkIdType row = 0; row != table->GetNumberOfRows(); row++)
	{
		allIds.push_back(columnData->GetVariantValue(row).ToTypeInt64());
	}

	
	for(int i = 0 ; i<ids.size() ;++i)
	{
		std::vector<int>::iterator posn1 = std::find(allIds.begin(), allIds.end(), ids.at(i));
		table->SetValueByName(posn1-allIds.begin(),classifier.c_str(),new_class);
	}

	this->updateViews();
    
	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "class changed for " << ids.size() << " objects." << std::endl;
	}
	#endif

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
		if(NucAdjTable)
		{
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
		}
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		selection->clear();
		this->updateViews();

		segView->SetNucAdjTable(NucAdjTable);

		std::string log_entry = "DELETE , ";
		for(int i=0; i<(int)ids.size(); ++i)
			log_entry += " " + ftk::NumToString(ids.at(i));
		log_entry += " , ";
		log_entry += ftk::TimeStamp();
		ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
	}
}

void NucleusEditor::updateMultiLabels(void)
{
	if(!nucSeg) return;
	if(!selection) return;
	std::vector<ObjectSelection::Point> * selected = selection->GetSelectedPoints();
	std::vector<int> times;
	std::vector<int> ids;
	std::vector<int> new_ids;

	for(int i=0; i< selected->size(); ++i)
	{
		times.push_back(selected->at(i).time);
		ids.push_back(selected->at(i).id);
		new_ids.push_back(selected->at(i).new_id);
	}
	nucSeg->ReassignLabels(times,ids,new_ids);
	segView->SetLabelImage(nucSeg->GetLabelImage(),selection);
	segView->SetCenterMapVectorPointer(nucSeg->centerMap4DImage);
	segView->SetBoundingBoxMapVectorPointer(nucSeg->bBoxMap4DImage);
	this->update5DTable();

}


void NucleusEditor::update5DTable(void)
{
	//if(!table) return;

	int currentTime = segView->GetCurrentTimeVal();
	table = nucSeg->table4DImage.at(currentTime);
	nucSeg->SetCurrentbBox(nucSeg->bBoxMap4DImage.at(segView->GetCurrentTimeVal()));

	//Update the table and scatter plot views
	if (tblWin.size()!=0)
		tblWin.back()->setModels(table,selection);
	if (pltWin.size()!=0)
		pltWin.back()->setModels(table,selection);

	this->clearSelections();
	this->updateViews();

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
		if(NucAdjTable)
			segView->SetNucAdjTable(NucAdjTable);

		for(unsigned i=0; i<new_grps.size(); ++i)
		{
			std::string log_entry = "MERGE , ";
			for(unsigned j=0; j<new_grps.at(i).size()-1; ++j)
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

void NucleusEditor::batchSplitCells(void)
{
	if(!nucSeg) return;

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());
	
	Split_Params_Dialog *dialog = new Split_Params_Dialog(this);
	int numSplitObjects;
	if( dialog->exec() )
	{
		numSplitObjects = dialog->getSplitNumber();		
	}
	else
	{
		delete dialog;
		return;
	}	

	std::vector< std::vector<int> > new_grps = nucSeg->BatchSplit(ids, numSplitObjects, table, NucAdjTable);
	if(new_grps.size() != 0)
	{
		projectFiles.outputSaved = false;
		projectFiles.tableSaved = false;
		projectFiles.adjTablesSaved = false;
		selection->clear();
		this->updateViews();

		for(int i=0; i<(int)new_grps.size(); ++i)
		{
			std::string log_entry = "SPLIT , ";
			log_entry += ftk::NumToString(new_grps[i][0]) + " , ";
			for(int j=1; j<(int)new_grps[i].size(); ++j)
			{
				log_entry += ftk::NumToString(new_grps[i][j]) + " ";
			}
			log_entry += " , ";
			log_entry += ftk::TimeStamp();
			ftk::AppendTextFile(projectFiles.GetFullLog(), log_entry);
		}
		#ifdef USE_QT_TESTING
		std::cout << "batch split" << std::endl;
		#endif
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
    #ifdef USE_QT_TESTING
    if(this->TestInputFile != "")
    {
      std::cout << "cell split" << std::endl;
    }
    #endif
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
	for ( std::set<long int>::iterator it=sels.begin(); it != sels.end(); it++ )
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
			#ifdef USE_QT_TESTING
			std::cout << "split along Z" << std::endl;
			#endif
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
	#ifdef USE_QT_TESTING
	if(this->TestInputFile != "")
	{
		std::cout << "cells filled" << std::endl;
	}
	#endif
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
		#ifdef USE_QT_TESTING
		if(this->TestInputFile != "")
		{
			std::cout << "margins excluded" << std::endl;
		}
	#endif
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
	editMenu->setEnabled(true);
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

void NucleusEditor::segmentByActiveContour()
{
	if(!myImg) return;
	editMenu->setEnabled(true);

	SomaExtractor *Somas = new SomaExtractor();	

	QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Param Files"),
                                                    QDir::currentPath(), tr("Option files (*.opt)"));
    if (!fileName.isEmpty())
    {
		Somas->LoadOptions( fileName.toStdString().c_str()); 
    }
	else
	{
		delete Somas;
		return;
	}
	
	std::vector< itk::SizeValueType > size = myImg->Size();
	unsigned char* ucharImagePtr = static_cast<unsigned char*> (myImg->GetDataPtr(0,0));
	
	std::cout<< "Generating seeds and get binary image.  "<<std::endl;
	std::cout<<(int)size[2]<<"\t"<<(int)size[3]<<"\t"<<(int)size[1]<<std::endl;

	std::vector< itk::Index<3> > seedVector;
	SomaExtractor::ProbImageType::Pointer binImagePtr = Somas->GenerateSeedPoints(ucharImagePtr, (int)size[2], (int)size[3], (int)size[1], seedVector);
	std::cout<< "Segmenting..."<<std::endl;
	/// SegmentSoma: Active Contour without GVF, eliminate small objects
	SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(seedVector, binImagePtr);
	if(segImage)
	{
		std::cout<< "Writing SomaLabeledImage.nrrd"<<std::endl;
		Somas->writeImage("_somaLabeledImage.nrrd", segImage);
		Somas->writeCentroids( "_somaCentroids.txt" ,seedVector);
	}

	labImg = ftk::Image::New();
	if(!labImg->LoadFile(std::string("_somaLabeledImage.nrrd")))
	{
		labImg = NULL;
	}
	else
	{
		std::cout<< "Visualize labeled image"<<std::endl;
		selection->clear();
		segView->SetLabelImage(labImg, selection);
		this->updateNucSeg();
	}
	delete Somas;
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
	//if(nucSeg->GetCenterMapPointer())
	//	pProc->SetCenterMapPointer(nucSeg->GetCenterMapPointer());
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
		//binImg = pProc->GetBinaryImage();
		segView->SetLabelImage(labImg,selection);
		this->updateNucSeg();
		
		if(labImg->GetImageInfo()->numTSlices>1)	// if multi time points get the table from nucseg after computing all geometries in updateNucSeg			
			table = nucSeg->table4DImage.at(segView->GetCurrentTimeVal());
		else
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
// A dialog to get the parameters for batch split:
//***********************************************************************************
Split_Params_Dialog::Split_Params_Dialog(QWidget *parent)
: QDialog(parent)
{
	//splitNumberLabel = new QLabel("Split Into ? ");
	//splitNumber = new QLineEdit();
	//splitNumber->setMinimumWidth(20);
	//splitNumber->setFocusPolicy(Qt::StrongFocus);
	//splitNumberLayout = new QHBoxLayout;
	//splitNumberLayout->addWidget(splitNumberLabel);
	//splitNumberLayout->addWidget(splitNumber);

	splitNumberLabel = new QLabel("Split Into ? ");
	splitNumber = new QComboBox();
	for(int v=2; v<=5; ++v)
	{
		splitNumber->addItem(QString::number(v));
	}
	splitNumberLayout = new QHBoxLayout;
	splitNumberLayout->addWidget(splitNumberLabel);
	splitNumberLayout->addWidget(splitNumber);
	
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	layout = new QVBoxLayout;
	layout->addLayout(splitNumberLayout);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Number Of Objects After Split"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

int Split_Params_Dialog::getSplitNumber()
{
	int num;
	QString input = splitNumber->currentText();
	num = input.toInt();
	return num;
}



//***************************************************************************
//***********************************************************************************
//***********************************************************************************
// A dialog to get the number of samples the user is willing to validate.
// Reusing the variables of the Confidence Dialog.
//***********************************************************************************
SamplePercentDialog::SamplePercentDialog(int no_of_samples,int numberOfClasses,QWidget *parent)
: QDialog(parent)
{
	sampleNumberLabel = new QLabel("The total number of cells present : "); 
	numberLabel = new QLabel(QString::number(no_of_samples));

	sampleLabel = new QLabel(" % of cells you are willing to validate : "); 
	samples = no_of_samples;
	class_number = numberOfClasses;

	QGridLayout * layout = new QGridLayout;

	QHBoxLayout *sampleLayout1 = new QHBoxLayout;
	sampleLayout1->addWidget(sampleNumberLabel);
	sampleLayout1->addWidget(numberLabel);


	QLabel *sampleLabel2 = new QLabel(" Number of cells you are willing to validate : "); 
	vSpinNumber = new QSpinBox(this);
	vSpinNumber->resize(vSpinNumber->minimumSizeHint());
	vSpinNumber->setRange(1,no_of_samples);
	vSpinNumber->setEnabled(true);


	QHBoxLayout *sampleLayout2 = new QHBoxLayout;
	sampleLayout2->addWidget(sampleLabel2);
	sampleLayout2->addWidget(vSpinNumber);

	vSpinPercent = new QDoubleSpinBox(this);
	vSpinPercent->resize( vSpinPercent->minimumSizeHint() );
	vSpinPercent->setRange(0,100);
	vSpinPercent->setEnabled(true);

	QHBoxLayout *sampleLayout = new QHBoxLayout;
	sampleLayout->addWidget(sampleLabel);
	sampleLayout->addWidget(vSpinPercent);

	layout->addLayout(sampleLayout1,1,0,0);
	layout->addLayout(sampleLayout,2,0,0);
	layout->addLayout(sampleLayout2,3,0,0);


	connect(vSpinPercent, SIGNAL(valueChanged(double)), this, SLOT(setNumber(double)));
	connect(vSpinNumber, SIGNAL(valueChanged(int)),this, SLOT(setPercent(int)));

	//Default values
	vSpinNumber->setValue(30);	


	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	layout->addLayout(bLayout,4,0,0);
	this->setLayout(layout);
	this->setWindowTitle(tr("Validation Sample %"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

SamplePercentDialog::SamplePercentDialog(QWidget *parent) 
: QDialog(parent)
{
	deltaLabel = new QLabel("Enter delta : ");
	delta = new QLineEdit();
	delta->setMinimumWidth(30);
	delta->setFocusPolicy(Qt::StrongFocus);
	deltaLayout = new QHBoxLayout;
	deltaLayout->addWidget(deltaLabel);
	deltaLayout->addWidget(delta);
	
	numbinLabel = new QLabel("Specify numbin: ");
	numbin = new QLineEdit();
	numbin->setMinimumWidth(30);
	numbin->setFocusPolicy(Qt::StrongFocus);
	numbinLayout = new QHBoxLayout;
	numbinLayout->addWidget(numbinLabel);
	numbinLayout->addWidget(numbin);
	
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	layout = new QVBoxLayout;
	layout->addLayout(deltaLayout);
	layout->addLayout(numbinLayout);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Parameter Setting"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

double SamplePercentDialog::getDelta()
{
	QString input = delta->displayText();
	double delta;
	delta = atof((input.toStdString()).c_str());
	return delta;
}

int SamplePercentDialog::getNumbin()
{
	QString input = numbin->displayText();
	int numbin;
	numbin = atoi((input.toStdString()).c_str());
	return numbin;
}
void SamplePercentDialog::setNumber(double x)
{
	disconnect(vSpinNumber, SIGNAL(valueChanged(int)),this, SLOT(setPercent(int)));
	//vSpinNumber->setValue((int)((x*samples/100)/class_number)*class_number);// Divide and multiply by class_number to ensure consistency
	vSpinNumber->setValue((int)(x*samples/100));// Divide and multiply by class_number to ensure consistency
	connect(vSpinNumber, SIGNAL(valueChanged(int)),this, SLOT(setPercent(int)));
}

void SamplePercentDialog::setPercent(int x)
{
	disconnect(vSpinPercent, SIGNAL(valueChanged(double)), this, SLOT(setNumber(double)));
	vSpinPercent->setValue(((double)x*100/samples));
	connect(vSpinPercent, SIGNAL(valueChanged(double)), this, SLOT(setNumber(double)));
}


//***************************************************************************
//***********************************************************************************
//***********************************************************************************
// A dialog to get the paramaters file to use and specify the channel if image has
// more than one:
//***********************************************************************************
QueryDialog::QueryDialog(int QueryType, QVector<QString> classes, bool diffusion_graph, QWidget *parent)
: QDialog(parent)
{
	layout = new QVBoxLayout;

	idLabel = new QLabel("Choose Object IDs ('0' for all IDs) :");
	IDs = new QLineEdit();
	IDs->setMinimumWidth(100);
	IDs->setFocusPolicy(Qt::StrongFocus);
	idLayout = new QHBoxLayout;
	idLayout->addWidget(idLabel);
	idLayout->addWidget(IDs);
	layout->addLayout(idLayout);

	if(QueryType == 1)
	{
		kLabel = new QLabel("Choose Number Of Nearest Neighbors K: ");
		K = new QLineEdit();
		K->setMinimumWidth(50);
		K->setFocusPolicy(Qt::StrongFocus);
		kLayout = new QHBoxLayout;
		kLayout->addWidget(kLabel);
		kLayout->addWidget(K);
		layout->addLayout(kLayout);
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
		layout->addLayout(radLayout);
	}	

	if(!diffusion_graph)
	{
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
		layout->addLayout(classLayout1);
		
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
		layout->addLayout(classLayout2);
	}

	//autoButton = new QRadioButton(tr("Automatic Parameter Selection"),this);
	//autoButton->setChecked(true);

	//fileButton = new QRadioButton(tr("Use Parameter File..."),this);

	//fileCombo = new QComboBox();
	//fileCombo->addItem(tr(""));
	//fileCombo->addItem(tr("Browse..."));
	//connect(fileCombo, SIGNAL(currentIndexChanged(QString)),this,SLOT(ParamBrowse(QString)));

	if(QueryType == 1)
	{
		check = new QCheckBox("Show k mutual graph");
		check->setChecked(false);
		layout->addWidget(check);
	}

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);	
	layout->addLayout(bLayout);

	//layout = new QVBoxLayout;
	//layout->addLayout(idLayout);
	//if(QueryType == 1)
	//	layout->addLayout(kLayout);
	//if(QueryType == 2)
	//	layout->addLayout(radLayout);
	//if(!diffusion_graph)
	//{
	//	layout->addLayout(classLayout1);
	//	layout->addLayout(classLayout2);
	//}
	//if(QueryType == 1)
	//	layout->addWidget(check);
	////layout->addWidget(autoButton);
	////layout->addWidget(quitButton);
	////layout->addWidget(okButton);
	//layout->addLayout(bLayout);
	this->setLayout(layout);
	if(QueryType == 1)
	{
		if(!diffusion_graph)
			this->setWindowTitle(tr("K Nearest Neighbors"));
		else
			this->setWindowTitle(tr("K Nearest Diffusion Neighbors"));
	}
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

bool QueryDialog::getKMutual()
{
	return check->isChecked();
}


//***************************************************************************
//***********************************************************************************
//***********************************************************************************
// A dialog to get the paramaters file to use and specify the channel if image has
// more than one:
//***********************************************************************************
ChangeClassDialog::ChangeClassDialog(QVector<QString> cls, QWidget *parent)
: QDialog(parent)
{

	classifiers = cls;
	classifierLabel = new QLabel("Choose Classifier: ");
	classifierCombo = new QComboBox();
	for(int v = 0; v<classifiers.size(); ++v)
	{
		classifierCombo->addItem(classifiers.at(v));
	}
	QHBoxLayout *classifierLayout = new QHBoxLayout;
	classifierLayout->addWidget(classifierLabel);
	classifierLayout->addWidget(classifierCombo);

	classLabel = new QLabel("Choose New Class Value: ");
	classCombo = new QComboBox();
	for(int v = 0; v<=10; ++v)
	{
		classCombo->addItem(QString::number(v));
	}
	QHBoxLayout *classLayout = new QHBoxLayout;
	classLayout->addWidget(classLabel);
	classLayout->addWidget(classCombo);

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
	layout->addLayout(classifierLayout);
	layout->addLayout(classLayout);
	//layout->addWidget(autoButton);
	//layout->addWidget(quitButton);
	//layout->addWidget(okButton);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Change Class"));
	
	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}



std::string ChangeClassDialog::getClassColumn()
{
	return classifiers[classifierCombo->currentIndex()].toStdString();
}

unsigned short ChangeClassDialog::getClass()
{
	return classCombo->currentIndex();
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
//*******************************************************************************************************************************
//Perform Spectral Unmixing of Channels
//*******************************************************************************************************************************
void NucleusEditor::unmixChannels(void)
{
	if(!myImg)	return;
	int nchannels = (int)myImg->GetImageInfo()->numChannels;
	if(nchannels<2) return;
	// define the unmixing mode
	ftk::SpectralUnmixing::UnmixMode umode = ftk::SpectralUnmixing::LINEAR_UNMIX;
	SpecUnmix = new ftk::SpectralUnmixing();
	SpecUnmix->SetInputImage(myImg);
	SpecUnmix->SetNumberOfChannels(nchannels);
	SpecUnmix->SetUnmixMode(umode);
	SpecUnmix->Update();
	myImg = SpecUnmix->GetOutput();
	segView->SetChannelImage(myImg);
	printf("done with unmixing\n");

}
	
//*******************************************************************************************************************************
//Load files from command line arguments
//*******************************************************************************************************************************
void NucleusEditor::parseArguments(QStringList args)
{
	for(int i = 1; i < args.size(); ++i)
	{
		QString arg = args[i];
		QString fileName;

		//worry about relative paths...
		if(QString::compare(arg, "--image", Qt::CaseInsensitive) == 0)
		{
			++i;
			fileName = args[i];
			this->loadImage(fileName);
		}
		if(QString::compare(arg, "--result", Qt::CaseInsensitive) == 0)
		{
			++i;
			fileName = args[i];
			this->loadResult(fileName);
		}
		if(QString::compare(arg, "--table", Qt::CaseInsensitive) == 0)
		{
			++i;
			fileName = args[i];
			this->loadTable(fileName);
		}
		if(QString::compare(arg, "--adjtables", Qt::CaseInsensitive) == 0)
		{
			++i;
			fileName = args[i];
			this->loadAdjTables(fileName);
		}
#ifdef USE_QT_TESTING
		if(QString::compare(arg, "--test", Qt::CaseInsensitive) == 0)
		{
			++i;
			this->TestInputFile = args[i];
		}
		if(QString::compare(arg, "--baseline", Qt::CaseInsensitive) == 0)
		{
			++i;
			this->TestBaselineImageFileName = args[i];
		}
#endif
    /*
    else if(QString::compare(arg, "--5dimage", Qt::CaseInsensitive) == 0)
    {
    }
    else if(QString::compare(arg, "--labelimage", Qt::CaseInsensitive) == 0)
    {
    }
    if(QString::compare(arg, "--project", Qt::CaseInsensitive) == 0)
    {
    }
    */
  }
}

//******************************************************************************
int NucleusEditor::runTest()
{
  #ifdef USE_QT_TESTING
  if( this->TestInputFile == "" )
  {
    return -1;
  }

  //setup test utility
  this->Tester->SetBaselineImage(
    this->TestBaselineImageFileName.toStdString().c_str() );

  //playback the test recording
  this->Tester->playTestFile( this->TestInputFile );
  
  //if this is an image comparison test
  if( this->TestBaselineImageFileName != "" )
  {
    //write our QImage to disk so the testing framework can read it back in
    QString tmpFile = TEST_OUTPUT_DIR;
    tmpFile += "/";
    tmpFile += "nucleus_test_output.png";
    segView->SaveDisplayImageToFile(tmpFile);

    //compare displayed image to baseline, then print & return the results of
    //this test
    if(this->Tester->compareResults(tmpFile) == false)
    {
      std::cout << "ERROR: test failed" << std::endl;
      return 1;
    }
    else
    {
      std::cout << "test passed" << std::endl;
    }
  }
  return 0;

  #endif
  return -1;
}

