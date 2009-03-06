#include "ControlBar.h"

using namespace std;
//*******************************************************************************
// ControlBar
//
// This widget is the main window.  It basically consists of a menu bar and a 
// status bar.  From these menus, all of the modules can be called.
// 
// The control bar also keeps track of some things like open images and the 
// data model that is being used.  Only one model is allowed at one time.  
// This is to help conserve memory, screen space, and allow for easier programming.
//********************************************************************************

ControlBar::ControlBar(const char *c)
{
	this->argv0 = new char[strlen(c)+1];
	strcpy(this->argv0, c);
	createMenus();
	createStatusBar();
	this->settings = new QSettings("RPI", "Farsight");
	this->pythonProcess = new QProcess(); 
	this->InitializePreferencesDialog();

	setWindowTitle(tr("Farsight Explorer"));

	Qt::WindowFlags flags = windowFlags();
	flags &= ~Qt::WindowMaximizeButtonHint;
	setWindowFlags(flags);

	setFixedSize (350, 50); 
	move(0,0);

	segResult = NULL;
	currentModel = NULL;

	module = NULL;

	loadedImages.clear();
	lastPath = ".";

	tblWin.clear();
	pltWin.clear();

	//Crashes when this is enabled!
	//setAttribute ( Qt::WA_DeleteOnClose );
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
void ControlBar::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));

	openAction = new QAction(tr("&Open Image..."), this);
    openAction->setShortcut(tr("Ctrl+O"));
    openAction->setStatusTip(tr("Open an existing image file"));
    connect(openAction, SIGNAL(triggered()), this, SLOT(loadImage()));
	fileMenu->addAction(openAction);

	openSeriesAction = new QAction(tr("Open Time Series..."), this);
	openSeriesAction->setShortcut(tr("Ctrl+T"));
	openSeriesAction->setStatusTip(tr("Open multi-dimensional images as a time series"));
	connect(openSeriesAction, SIGNAL(triggered()), this, SLOT(loadImageSeries()));
	fileMenu->addAction(openSeriesAction);

	xmlAction = new QAction(tr("Load Result..."), this);
	xmlAction->setStatusTip(tr("Open an XML result file"));
	connect(xmlAction,SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu->addAction(xmlAction);

	loadOutliersAction = new QAction(tr("Load Outliers"), this);
	loadOutliersAction->setStatusTip(tr("Load Outliers in outliers.txt"));
	connect(loadOutliersAction, SIGNAL(triggered()), this, SLOT(loadOutliers()));
	fileMenu->addAction(loadOutliersAction);

	saveAction = new QAction(tr("Save Result"), this);
	saveAction->setStatusTip(tr("Save Changes (Edits, etc)"));
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(saveResult()));
	fileMenu->addAction(saveAction);

	testAction = new QAction(tr("Test Function"),this);
	connect(testAction,SIGNAL(triggered()), this, SLOT(test()));
	//fileMenu->addAction(testAction);

    exitAction = new QAction(tr("Exit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
	fileMenu->addSeparator();
    fileMenu->addAction(exitAction);

  //EDIT MENU
	editMenu = menuBar()->addMenu(tr("&Edit"));
	preferencesAction = new QAction(tr("Preferences"), this);
	preferencesAction->setStatusTip(tr("Edit Farsight Preferences"));
	connect(preferencesAction,SIGNAL(triggered()),this,SLOT(EditPreferences()));
	editMenu->addAction(preferencesAction);

	//VIEW MENU
	viewMenu = menuBar()->addMenu(tr("&View"));
	newScatterAction = new QAction(tr("New Scatter"), this);
	newScatterAction->setStatusTip(tr("Open a new Scatterplot Window"));
	connect(newScatterAction,SIGNAL(triggered()),this,SLOT(CreateNewPlotWindow()));
	viewMenu->addAction(newScatterAction);

	viewMenu->addSeparator();

	outlierAction = new QAction(tr("Show Outliers"), this);
	outlierAction->setStatusTip(tr("Show outliers in views"));
	outlierAction->setCheckable(true);
	connect(outlierAction, SIGNAL(triggered()), this, SLOT(showOutliers()));
	viewMenu->addAction(outlierAction);

	viewMenu->addSeparator();

	refreshAction = new QAction(tr("refresh/reload"), this);
	connect(refreshAction, SIGNAL(triggered()), this, SLOT(refreshViews()));
	viewMenu->addAction(refreshAction);

	viewMenu->addSeparator();

	pythonAction = new QAction(tr("Open Python Window"), this);
	pythonAction->setStatusTip(tr("Start your favorite python interpreter"));
	pythonAction->setShortcut(tr("Ctrl+P"));
	connect(pythonAction, SIGNAL(triggered()), this, SLOT(OpenPythonWindow()));
	viewMenu->addAction(pythonAction);

	//SEGMENT MENU
	segmentMenu = menuBar()->addMenu(tr("Segment"));

//*******************************************
// ADD NEW MODULES HERE
// FOLLOW BUILD_YOUSEF EXAMPLE
//*******************************************
#ifdef BUILD_NUCLEI
	yousefAction = new QAction(tr("Nuclei"),this);
	yousefAction->setStatusTip(tr("Click to start a segmenting nuclei"));
	connect(yousefAction,SIGNAL(triggered()), this, SLOT(newYousefSeg()));
	segmentMenu->addAction(yousefAction);
#endif

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
void ControlBar::about()
{
	    QMessageBox::about(this, tr("About FARSIGHT"),
            tr("<h2>FARSIGHT</h2>"
               "Renssalear Polytechnic Institute"
               ));
}

//******************************************************************************
// Because ControlBar inherits from QMainWindow, a status bar is provided
// Here we just show a message in the status bar when loading
//******************************************************************************
void ControlBar::createStatusBar()
{
    statusLabel = new QLabel(tr(" Ready"));
    statusBar()->addWidget(statusLabel, 1);
}

//******************************************************************************
// SLOT: loadImage()
//  This function uses QFileDialog to get the filename of the image to open
//  If a filename is selected, a new imagewindow is created and the image is 
//  loaded and displayed
//******************************************************************************
void ControlBar::loadImage()
{
	QString filename  = QFileDialog::getOpenFileName(this,"Choose an image",lastPath, 
			tr("TIF Files (*.tif *.tiff)\n" 
			   "PIC Files (*.pic)\n" 
			   "LSM Files (*.lsm)\n"
			   "All Files (*.*)"));

    if(filename != "")
    {  
		//The segmentationWindow requires a model of some type to be created.
		//If we already have a model then we can use it.  Otherwise we need a new one.

		lastPath = QFileInfo(filename).absolutePath();

		//vector<string> oneImage(0);
		//oneImage.push_back(filename.toStdString());
		ftk::Image *newImg = NewFTKImage(filename.toStdString());
		if (newImg)
		{
			SegmentationWindow *segWin = new SegmentationWindow();
			connect(segWin, SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
			segWin->SetChannelImage(newImg);
			segWin->show();
		}
		else
		{
			std::cerr << "Couldn't load Image" << std::endl;
		}
    }
}

//******************************************************************************
// SLOT: loadImageSeries()
// This function creates a new dialog to get the image series information,
// verifies that the info is true, and then attempts to load all of the images
// into a single ftkImage
//******************************************************************************
void ControlBar::loadImageSeries(void)
{
	QStringList fileNames = QFileDialog::getOpenFileNames(
                             this, "Select files to open as time series", lastPath,
                             tr("Images (*.tif *.tiff *.pic)"));
	
	if (fileNames.size() > 0 )
	{
		lastPath = QFileInfo(fileNames[0]).absolutePath();

		vector<string> images(0);

		QStringList list = fileNames;
		QStringList::Iterator it = list.begin();
		while(it != list.end()) 
		{
			images.push_back( (*it).toStdString() );
			++it;
		}

		ftk::Image *newImg = NewFTKImage(images);
		if (newImg)
		{
			SegmentationWindow *segWin = new SegmentationWindow();
			connect(segWin, SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
			segWin->SetChannelImage(newImg);
			segWin->show();
		}
		else
		{
			std::cerr << "Couldn't load Image" << std::endl;
		}
	}
}
//******************************************************************************
// Creates a new FTKImage and returns a pointer to it
//******************************************************************************
ftk::Image * ControlBar::NewFTKImage(std::vector<std::string> filenames)
{
	ftk::Image *img = new ftk::Image();

	if ( /*img->load(filenames)*/0 )
	{
		loadedImages.push_back(img);
	}
	else
	{
		delete img;
		img = NULL;
	}
	return img;
}

ftk::Image * ControlBar::NewFTKImage(std::string filename)
{
	ftk::Image *img = new ftk::Image();
	bool forDisplay = true;
	if( img->LoadFile(filename, forDisplay) )
		return img;
	else
	{
		delete img;
		return NULL;
	}
}

//******************************************************************************
// Checks to see if the image with filename "filename" is already loaded
// return the index if it is
//******************************************************************************
int ControlBar::isLoaded(std::string filename)
{
	for (int i = 0; i<loadedImages.size(); ++i)
	{
		vector<string> filenames;// = loadedImages.at(i)->GetFilenames();
		for (int f = 0; f < filenames.size(); f++)
		{
			if(filenames[f] == filename)
			{
				return i;
			}
		}
	}
	return int(-1);
}

void ControlBar::displayHelp()
{/*
	QTextBrowser text;
	text.setSource(QUrl("file:///C:/FARSIGHT/html/aboutfarsight.html"));
	text.show();*/	
}

//******************************************************************************
//A function (slot) for implementing test utilities
//******************************************************************************
void ControlBar::test(void)
{
	/*
	QHelpEngineCore helpEngine("C:\FARSIGHT\html\ftkhelp.qch");

	// get all file references for the identifier
	QMap<QString, QUrl> links = helpEngine.linksForIdentifier(QLatin1String("FTK::about"));
	//QMap<QString, QUrl> links = 

	// If help is available for this keyword, get the help data
	// of the first file reference.
	if (links.count()) 
	{
		QByteArray helpData = helpEngine.fileData(links.constBegin().value());
		// show the documentation to the user
		if (!helpData.isEmpty());
			displayHelp(helpData);
	}
	*/
	displayHelp();
}

//******************************************************************************
//Reimplement closeEvent to also close all other windows in the application
//******************************************************************************
void ControlBar::closeEvent(QCloseEvent *event)
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
				segResult->SaveAll();
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

void ControlBar::saveResult()
{
	if(segResult)
	{
		if(segResult->editsNotSaved)
		{
			segResult->SaveAll();
		}
	}
}

//*************************************************************************
// THIS SLOT IS CALLED WHEN THE MODULE OR RESULT VIEWER IS CLOSED
//*************************************************************************
void ControlBar::closeWidget(QWidget* widget)
{
	if(widget == module)
	{
		//delete module;
		module = NULL;
	}
}

//***************************************************************************
//  THIS FUNCTION CLEARS THE MODEL/SEGMENTATION BY CLOSING OPEN WINDOW 
//  GROUPS ASSOCIATED WITH THEM, THEN DELETING THE MODEL AND SEGMENTATION
//***************************************************************************
void ControlBar::clearModel(void)
{
	if(module)
	{	//This should close all open windows for existing module
		module->close();
		delete module;
		module = NULL;
		closeWidget(module);
	}

	if(currentModel)
	{
		delete currentModel;
		currentModel = NULL;
	}
}

//*********************************************************************************
// This function initializes the model and selection model
//*********************************************************************************
void ControlBar::newModel(void)
{
	if(currentModel)
		clearModel();

	if(segResult)
		currentModel = new SegmentationModel(segResult);

}

//*****************************************************************
// THIS FUNCTION STARTS A MODULE. 
// IT SHOULD BE CALLED AFTER INITIALIZING THE MODULE
//*****************************************************************
void ControlBar::startModule(void)
{
	newModel();

	//Create a list of loaded images(by filename)
	vector<string> filenames(0);
	for (int i = 0; i<loadedImages.size(); ++i)
	{
		vector<string> myNames;// = loadedImages.at(i)->GetFilenames();
		for (int f = 0; f < myNames.size(); f++)
		{
			filenames.push_back(myNames[f]);
		}
		
	}

	if(module)
		delete module;

	module = new ModuleWidget(filenames);
	
	connect(module, SIGNAL(closing(QWidget*)),this, SLOT(closeWidget(QWidget*)));
	module->show();
}

//******************************************************************************
// This function loads a segmentation result from XML
// The XML file should tell where to find the original image/data
//   It also contains all of the feature and segmentation information
//   As well as what type of segmentation has been performed
//******************************************************************************
void ControlBar::loadResult(void)
{
	QString filename  = QFileDialog::getOpenFileName(this,"Choose a Result",lastPath, 
			tr("XML Files (*.xml)\n"));

    if(filename == "")
		return;

	if(module)
	{
		delete module;
		module = NULL;
	}
	QString path = QFileInfo(filename).absolutePath();
	QString name = QFileInfo(filename).baseName();

	lastPath = path;

	//Load up a generic segmentationResult and display it!
	//segResult = new ftk::SegmentationResult( path.toStdString(), name.toStdString() );
	segResult = new ftk::NuclearSegmentation( path.toStdString(), name.toStdString() );
	
	if ( !segResult->RestoreFromXML() )
		std::cerr << segResult->GetErrorMessage() << std::endl;

	//Now I have objects stored in memory - put the features into the model
	newModel();
	CreateNewTableWindow();
	CreateNewPlotWindow();
	SegmentationWindow *segwin = CreateNewSegmentationWindow();

	if (segResult->LoadData())
	{
		segwin->SetChannelImage(segResult->getDataImage());
	}
	if(segResult->LoadLabel())
	{
		segwin->SetLabelImage(segResult->getLabelImage());
	}
	segwin->show();
}

void ControlBar::refreshViews()
{
	if( !currentModel )
		return;

	if( outlierAction->isChecked() && currentModel->HasOutliers() )
	{
		loadOutliers();
	}
}

void ControlBar::loadOutliers()
{
	if( !currentModel )
		return;

	string outlierfilename = currentModel->SegResult()->PrependProjectPath("outliers.txt");

	//Check to see if I have outliers and if new outliers are available
	if( currentModel->HasOutliers() )
	{
		QDateTime lastModTime = QFileInfo(QString::fromStdString(outlierfilename)).lastModified();
		QMessageBox::information(this, tr("Times"), lastModTime.toString() + "\n" + outlierModTime.toString() );
		if(lastModTime <= outlierModTime)	//Load is not required
		{
			outlierAction->setChecked(true);
			return;
		}
	}

	//Load is required:
	outlierModTime = QFileInfo(QString::fromStdString(outlierfilename)).lastModified();

	ifstream inFile; 
	inFile.open(outlierfilename.c_str() );
	if ( !inFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outlierfilename << std::endl;
		return;
	}

	const int MAXLINESIZE = 512;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];

	vector< int > outliers(0); 
	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() ) //Get all outliers
	{
		char * pch = strtok (line," ");
		int num = (int)atof(pch);
		pch = strtok (NULL, " ");
		int id = (int)atof(pch);
		pch = strtok (NULL, " ");
		float val = atof(pch);

		outliers.push_back( id );
	
		inFile.getline(line, MAXLINESIZE);
	}
	inFile.close();

	//Now keep 1/10 of them
	int n = (int)outliers.size();
	outliers.resize(n/10);

	currentModel->SetOutliers( outliers );
	outlierAction->setChecked(true);
	currentModel->ShowOutliers(true);

}

void ControlBar::showOutliers()
{
	if( !currentModel )
	{
		outlierAction->setChecked(false);
		return;
	}

	if( outlierAction->isChecked() )
	{
		if( currentModel->HasOutliers() )
		{
			 currentModel->ShowOutliers(true);
		}
		else
		{
			QMessageBox::information(this, tr("MESSAGE"), tr("Outlier Have not been Loaded"));
			outlierAction->setChecked(false);
		}
	}
	else
	{
		//std::cerr << "Turning off outliers" << std::endl;
		currentModel->ShowOutliers(false);
	}
}

//******************************************************************************
//******************************************************************************
// ADD SLOTS TO START MODULES/PACKAGES HERE
//******************************************************************************
//******************************************************************************
#ifdef BUILD_NUCLEI
void ControlBar::newYousefSeg(void)
{
	clearModel();
	startModule();
}
#endif

//******************************************************************************
// Set up the edit preferences dialog
//******************************************************************************
void ControlBar::InitializePreferencesDialog()
  {
  this->preferencesDialog = new QDialog(this);
  this->preferencesDialog->setModal(true);
	this->preferencesDialog->setWindowTitle(tr("Edit preferences"));

	this->preferencesLayout = new QGridLayout();

	this->pythonLabel = new QLabel(tr("Python interpreter"));
	this->preferencesLayout->addWidget(pythonLabel,0,0);

  this->currentPythonLabel =
    new QLabel(settings->value("python/window").toString());
  this->preferencesLayout->addWidget(currentPythonLabel,0,1);

	this->browseForPythonButton = new QPushButton(tr("Browse"));
	connect(browseForPythonButton,SIGNAL(clicked()),this,
          SLOT(BrowseForPythonExecutable()));
	this->preferencesLayout->addWidget(browseForPythonButton,0,2);

  this->submitPreferencesButton = new QPushButton(tr("OK"));
	connect(submitPreferencesButton,SIGNAL(clicked()),
          preferencesDialog, SLOT(accept()));
  this->preferencesLayout->addWidget(submitPreferencesButton,1,0);

  this->cancelPreferencesButton = new QPushButton(tr("Cancel"));
	connect(cancelPreferencesButton,SIGNAL(clicked()),
          preferencesDialog, SLOT(reject()));
  this->preferencesLayout->addWidget(cancelPreferencesButton,1,1);

  this->preferencesDialog->setLayout(preferencesLayout);
  }

//******************************************************************************
// Open a dialog that allows the user to change application preferences
//******************************************************************************
void ControlBar::EditPreferences()
{
  this->preferencesDialog->exec();
}

//******************************************************************************
// Ask the user to select a python interpreter to use with Farsight
//******************************************************************************
bool ControlBar::BrowseForPythonExecutable()
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
bool ControlBar::ConfirmClosePython()
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

//******************************************************************************
// Open a python interpreter.  Ask the user where one is located if this
// information hasn't been specified yet.
//******************************************************************************
void ControlBar::OpenPythonWindow()
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

	QString path1Cmd = QString("import sys;sys.path.append('") + pythonFiles + QString("');");
	QString path2Cmd = QString("import os;os.environ['PATH'] = '") + exeFiles + QString("';");
	//QString path2Cmd = QString("sys.path.append('") + exeFiles + QString("');");
	QString importCmds = QString("from farsightutils import *;");
	QString printCmd = QString("print 'FARSIGHT ENVIRONMENT';");
	QString arg = QString(" -i -c ") + quote + path1Cmd + path2Cmd + importCmds + printCmd + quote;
	cmd.append(arg);
  
	this->pythonProcess->startDetached(cmd);
}

/*
//-----------------------------------------------------------------------------
void ControlBar::initPythonInterpretor()
{
  QString initStr = QString(
    "import sys\n"
    "sys.path.append('%1')\n"
    "import os\n"
    "os.environ['PATH'] = '%2'\n"
    "sys.path.append('%2')\n"
    "from farsightutils import *\n")
    .arg(this->pythonFiles)
    .arg(this->exeFiles);
  this->PythonDialog->print("FARSIGHT ENVIRONMENT");
  this->PythonDialog->runString(initStr);
  this->PythonDialog->setAttribute(Qt::WA_QuitOnClose, false);
}
*/

//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void ControlBar::CreateNewPlotWindow(void)
{
	if(!currentModel)
		return;

	pltWin.push_back(new PlotWindow(currentModel));
	connect(pltWin.back(), SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	pltWin.back()->show();
}

//******************************************************************************
// Create a new table window
//******************************************************************************
void ControlBar::CreateNewTableWindow(void)
{
	if(!currentModel)
		return;

	tblWin.push_back(new TableWindow(currentModel));
	connect(tblWin.back(), SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	tblWin.back()->ResizeToOptimalSize();
	tblWin.back()->show();
}

//******************************************************************************
// Create a new segmentation window.  Shows Images and Segmentation results.
//******************************************************************************
SegmentationWindow* ControlBar::CreateNewSegmentationWindow(void)
{
	segWin.push_back( new SegmentationWindow() );
	connect(segWin.back(), SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	segWin.back()->setWindowTitle(tr("Segmentation Results Viewer"));
	
	if(currentModel)
		segWin.back()->SetModels(currentModel);

	return segWin.back();
}

/*


void ControlBar::startPython(void)
{
	MyThread *t = new MyThread;
	t->start();
	//int i = system("C:/FARSIGHT/bin64/tinyxml/debug/XMLtoMETA C:/FARSIGHT/documents/SupplementE1.xml");
	//system("gui");
	//std::cerr << "return system value is: " << i << std::endl;
}

class MyThread : public QThread
{
    public:
      void run();
};

void MyThread::run()
{
	//system("gui");
	system("python");
    exec();
}
*/
