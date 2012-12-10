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

//////////////////////////////////////////////////////////////////////////////////////////
/// Plot window uses scatterview to display a scatterplot of the data in a model.
/// In order for this view to work the model must be set.
/// The window allows the user to change the columns used in the scatterplot.
//////////////////////////////////////////////////////////////////////////////////////////
#include "PlotWindow.h"

//for DBL_MAX...
#include <float.h>

//!Constructor
PlotWindow::PlotWindow(QWidget *parent)
: QMainWindow(parent)
{
	this->setupUI();
}

void PlotWindow::setupUI(void)
{
	/**
	* Set up the GUI for the Scatterplot
	* This provides signals and slots for changing what is displayed
	*/
	resize(500, 500);

	QWidget *centralWidget = new QWidget();
	QVBoxLayout *vlayout = new QVBoxLayout();
	scatter = new ScatterView(this);
	vlayout->addWidget(scatter);
	centralWidget->setLayout(vlayout);
	this->setCentralWidget(centralWidget);

	//!Setup menu:
	optionsMenu = menuBar()->addMenu(tr("&Options"));

	normalizeAction = new QAction(tr("&Normalize"),this);
	normalizeAction->setCheckable(true);
	normalizeAction->setChecked(false);
	normalizeAction->setShortcut(tr("Ctrl+N"));
	connect(normalizeAction,SIGNAL(toggled(bool)), scatter, SLOT(SetNormalize(bool)));
	optionsMenu->addAction(normalizeAction);
//!action to change x
	optionsMenu->addSeparator();
	xMenu = new QMenu(tr("Set X Axis"));
	connect(xMenu, SIGNAL(triggered(QAction *)), this, SLOT(xChange(QAction *)));
	menuBar()->addMenu(xMenu);
//! action to change y
	yMenu = new QMenu(tr("Set Y Axis"));
	connect(yMenu, SIGNAL(triggered(QAction *)), this, SLOT(yChange(QAction *)));
	menuBar()->addMenu(yMenu);
//! action to change color
	colorMenu = new QMenu(tr("Set Color Column"));
	connect(colorMenu, SIGNAL(triggered(QAction *)), this, SLOT(colorChange(QAction *)));
	optionsMenu->addMenu(colorMenu);

	optionsMenu->addSeparator();

	clearAction = new QAction(tr("&Clear Selections"), this);
	clearAction->setShortcut(tr("Ctrl+C"));
	connect(clearAction, SIGNAL(triggered()), scatter, SLOT(clearSelections()));
	optionsMenu->addAction(clearAction);

	windowAction = new QAction(tr("Window Operation"), this);
	windowAction->setShortcut(tr("Ctrl+W"));
	connect(windowAction, SIGNAL(triggered()), scatter, SLOT(windowChange()));
	optionsMenu->addAction(windowAction);

	setWindowTitle(tr("Scatter Plot"));
	// If we do the following, the program crashes when we first close the scatterplot
	// and then close the  image viewer. 
	//setAttribute ( Qt::WA_DeleteOnClose );
}

void PlotWindow::setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels)
{
	scatter->setModels(tbl,sels);
	updateOptionMenus(tbl);
}

//!Redefine update:
void PlotWindow::update()
{
	//!detects changes in the linked space 
	vtkSmartPointer<vtkTable> tbl = scatter->GetTable();
	if(tbl)
		this->updateOptionMenus( tbl );
	if(scatter)
		scatter->update();

	QWidget::update();
}

void PlotWindow::closeEvent(QCloseEvent *event)
{
	//! Detects window closing and emits signal 
	emit closing(this);
	event->accept();
}

void PlotWindow::xChange(QAction *action)
{
	//! Sets the data column for the X axis
	scatter->SetColForX( action->toolTip().toInt() );
	action->setChecked(true);
}

void PlotWindow::yChange(QAction *action)
{
	//! Sets the data column for the Y axis
	scatter->SetColForY( action->toolTip().toInt() );
	action->setChecked(true);
}

void PlotWindow::colorChange(QAction *action)
{
	//! sets the color coding of the scatter plot
	scatter->SetColForColor( action->toolTip().toInt() );
	action->setChecked(true);
}

void PlotWindow::updateOptionMenus(vtkSmartPointer<vtkTable> tbl)
{
	//! populates the menues  
	if(!scatter) return;

	int xc = scatter->ColForX();
	int yc = scatter->ColForY();
	int cc = scatter->ColForColor();

	//!Add a new Action for each column for each menu item:
	xMenu->clear();
	yMenu->clear();
	colorMenu->clear();

	QActionGroup *xGroup = new QActionGroup(this);
	QActionGroup *yGroup = new QActionGroup(this);
	QActionGroup *cGroup = new QActionGroup(this);

	for (int c=1; c<tbl->GetNumberOfColumns(); ++c)
	{
		QString name = QString(tbl->GetColumnName(c));
		QAction *xAct = new QAction( name, this );
		xAct->setToolTip( QString::number(c) );
		xAct->setCheckable(true);
		QAction *yAct = new QAction( name, this );
		yAct->setToolTip( QString::number(c) );
		yAct->setCheckable(true);
		QAction *cAct = new QAction( name, this );
		cAct->setToolTip( QString::number(c) );
		cAct->setCheckable(true);

		xMenu->addAction(xAct);
		xGroup->addAction(xAct);
		yMenu->addAction(yAct);
		yGroup->addAction(yAct);
		colorMenu->addAction(cAct);
		cGroup->addAction(cAct);

		if(c==xc)
			xAct->setChecked(true);
		if(c==yc)
			yAct->setChecked(true);
		if(c==cc)
			cAct->setChecked(true);
	}
}

/*
void PlotWindow::modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
	updateOptionMenus();
}
*/
