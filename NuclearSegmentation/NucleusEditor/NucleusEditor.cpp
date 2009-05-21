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

	segResult = NULL;
	currentModel = NULL;

	lastPath = ".";

	tblWin.clear();
	pltWin.clear();

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

	xmlAction = new QAction(tr("Load Result..."), this);
	xmlAction->setStatusTip(tr("Open an XML result file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu->addAction(xmlAction);

	segmentAction = new QAction(tr("Segment Image..."), this);
	segmentAction->setStatusTip(tr("Starts the Nuclear Segmenation Wizard"));
	connect(segmentAction,SIGNAL(triggered()),this,SLOT(segmentImage()));
	fileMenu->addAction(segmentAction);

	saveAction = new QAction(tr("Save Result"), this);
	saveAction->setStatusTip(tr("Save Changes (Edits, etc)"));
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(saveResult()));
	fileMenu->addAction(saveAction);

    exitAction = new QAction(tr("Exit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
	fileMenu->addSeparator();
    fileMenu->addAction(exitAction);

	//VIEW MENU
	viewMenu = menuBar()->addMenu(tr("&View"));
	newScatterAction = new QAction(tr("New Scatter"), this);
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

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
				this->saveResult();
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

void NucleusEditor::saveResult()
{
	if(segResult)
	{
		if(segResult->editsNotSaved)
		{
			segResult->WriteToXML( QFileDialog::getSaveFileName(this, tr("Save File As"), lastPath, tr("XML (*.xml)") ).toStdString() );
			segResult->SaveLabel();
		}
	}
}

//***************************************************************************
//  THIS FUNCTION CLEARS THE MODEL/SEGMENTATION BY CLOSING OPEN WINDOW 
//  GROUPS ASSOCIATED WITH THEM, THEN DELETING THE MODEL AND SEGMENTATION
//***************************************************************************
void NucleusEditor::clearModel(void)
{
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

	if(segResult)
		currentModel = new SegmentationModel(segResult);
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

	//Load up a generic segmentationResult and display it!
	//segResult = new ftk::SegmentationResult( path.toStdString(), name.toStdString() );
	segResult = new ftk::NuclearSegmentation();
	
	if ( !segResult->RestoreFromXML(filename.toStdString()) )
	{
		std::cerr << segResult->GetErrorMessage() << std::endl;
		return;
	}

	//Now I have objects stored in memory - put the features into the model
	newModel();
	CreateNewTableWindow();
	CreateNewPlotWindow();
	segWin = CreateNewSegmentationWindow();
	segWin->SetChannelImage(segResult->getDataImage());
	segWin->SetLabelImage(segResult->getLabelImage());
	this->setCentralWidget(segWin);
	this->show();
	//segwin->show();
}

void NucleusEditor::segmentImage()
{
	NuclearSegmentationWizard *wizard = new NuclearSegmentationWizard(this);
	wizard->show();
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
// Create a new segmentation window.  Shows Images and Segmentation results.
//******************************************************************************
SegmentationWindow* NucleusEditor::CreateNewSegmentationWindow(void)
{
	SegmentationWindow *sWin = new SegmentationWindow();
	sWin->setWindowTitle(tr("Segmentation Results Viewer"));
	
	if(currentModel)
		sWin->SetModels(currentModel);

	return sWin;
}