#include "SampleEditor.h"

//*******************************************************************************
// SampleEditor
//********************************************************************************

SampleEditor::SampleEditor(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{
	model = NULL;
	selModel = NULL;
	plot = NULL;
	histo = NULL;

	createMenus();
	createStatusBar();

	table = new QTableView();
	//table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	//table->setSelectionBehavior( QAbstractItemView::SelectRows );
	setCentralWidget(table);

	setWindowTitle(tr("Sample Editor"));

	this->resize(500,500);
}

//******************************************************************************
// Here we just show a message in the status bar when loading
//******************************************************************************
void SampleEditor::createStatusBar()
{
    QLabel *statusLabel = new QLabel(tr(" Ready"));
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
void SampleEditor::createMenus()
{
	//FIRST HANDLE FILE MENU
	fileMenu = menuBar()->addMenu(tr("&File"));

	loadAction = new QAction(tr("Load File..."), this);
	loadAction->setStatusTip(tr("Load to table from text file"));
	connect(loadAction, SIGNAL(triggered()), this, SLOT(loadFile()));
	fileMenu->addAction(loadAction);
}

void SampleEditor::loadFile()
{
	//QString filename  = QFileDialog::getOpenFileName(this,"Choose a File", ".", 
	//		tr("TXT Files (*.txt)\n"
	//		   "SWC Files (*.swc)\n"
	//		   "All Files (*.*)"));

	//if(filename == "")
	//	return;

	//READ FILE HERE!!!!!
	std::vector<QString> headers;
	for(int i=0; i<10; ++i)
	{
		headers.push_back("Header " + QString::number(i));
	}
	
	int numRows = 100;
	std::vector< std::vector< double > > data;
	for(int i=0; i<numRows; ++i)
	{
		std::vector<double> row;
		row.assign(10,i);
		data.push_back(row);
	}

	//Now put data into Model:
	if(selModel)
		delete selModel;

	if(model)
		delete model;

	model = new QStandardItemModel;
	selModel = new QItemSelectionModel(model);

	model->clear();
	model->setColumnCount(headers.size());

	for(int i=0; i<(int)headers.size(); ++i)
	{
		model->setHeaderData(i, Qt::Horizontal, headers.at(i));
	}

	for (int row=0; row<(int)numRows; ++row)
	{
		model->insertRow(row);
		for(int col=0; col<(int)headers.size(); ++col)
		{
			model->setData(model->index(row, col), data.at(row).at(col));
		}
	}

	table->setModel(model);
	table->setSelectionModel(selModel); 
	table->update();

	if(plot)
	{
		plot->close();
		delete plot;
	}
	plot = new PlotWindow(selModel);
	plot->show();

	if(histo)
	{
		histo->close();
		delete histo;
	}
	histo = new HistoWindow(selModel);
	histo->show();

}
