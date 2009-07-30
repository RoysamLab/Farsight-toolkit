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

//****************************************************************************************
// Plot window uses scatterview to display a scatterplot of the data in a model.
// In order for this view to work the model must be set.
// The window allows the user to change the columns used in the scatterplot.
//****************************************************************************************
#include "PlotWindow.h"

//for DBL_MAX...
#include <float.h>
//Constructor
PlotWindow::PlotWindow(QItemSelectionModel *mod, QWidget *parent)
: QMainWindow(parent)
{
	this->setupUI();

	this->scatter->setModel( (QAbstractItemModel*)mod->model() );
	this->scatter->setSelectionModel(mod);

	this->updateOptionMenus(true);

	connect(mod->model(), SIGNAL(dataChanged(const QModelIndex &, const QModelIndex &)), this, SLOT(modelChange(const QModelIndex &, const QModelIndex &)));
	svmWidget = NULL;
}

void PlotWindow::setupUI(void)
{
	resize(500, 500);

	//Setup menu:
	optionsMenu = menuBar()->addMenu(tr("&Options"));
	xMenu = new QMenu(tr("Set X Axis"));
	connect(xMenu, SIGNAL(triggered(QAction *)), this, SLOT(xChange(QAction *)));
	optionsMenu->addMenu(xMenu);
	yMenu = new QMenu(tr("Set Y Axis"));
	connect(yMenu, SIGNAL(triggered(QAction *)), this, SLOT(yChange(QAction *)));
	optionsMenu->addMenu(yMenu);
	colorMenu = new QMenu(tr("Set Color Column"));
	connect(colorMenu, SIGNAL(triggered(QAction *)), this, SLOT(colorChange(QAction *)));
	optionsMenu->addMenu(colorMenu);

	toolsMenu = menuBar()->addMenu(tr("&Tools"));
	svmAction = new QAction(tr("SVM"), this);
	connect(svmAction, SIGNAL(triggered()), this, SLOT(startSVM()));
	toolsMenu->addAction(svmAction);

	QWidget *centralWidget = new QWidget();
	QVBoxLayout *vlayout = new QVBoxLayout();
	scatter = new ScatterView();
	vlayout->addWidget(scatter);
	centralWidget->setLayout(vlayout);
	this->setCentralWidget(centralWidget);

	setWindowTitle(tr("Scatter Plot"));
	// If we do the following, the program crashes when we first close the scatterplot
	// and then close the  image viewer. 
	//setAttribute ( Qt::WA_DeleteOnClose );
}

void PlotWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
}

void PlotWindow::xChange(QAction *action)
{
	scatter->SetColForX(action->toolTip().toInt(), action->text().toStdString());
	action->setChecked(true);
}

void PlotWindow::yChange(QAction *action)
{
	scatter->SetColForY(action->toolTip().toInt(), action->text().toStdString());
	action->setChecked(true);
}

void PlotWindow::colorChange(QAction *action)
{
	scatter->SetColForColor( action->toolTip().toInt() );
	action->setChecked(true);
}

void PlotWindow::updateOptionMenus(bool first = false)
{
	QAbstractItemModel *model = scatter->model();

	int xc = 0;
	int yc = 1;
	int cc = 0;
	if(!first)
	{
		xc = scatter->ColForX();
		yc = scatter->ColForY();
		cc = scatter->ColForColor();
	}

	//Add a new Action for each column for each menu item:
	xMenu->clear();
	yMenu->clear();
	colorMenu->clear();

	QActionGroup *xGroup = new QActionGroup(this);
	QActionGroup *yGroup = new QActionGroup(this);
	QActionGroup *cGroup = new QActionGroup(this);

	for (int c=0; c<model->columnCount(); ++c)
	{
		QString name = model->headerData(c,Qt::Horizontal).toString();
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

		if(first)
		{
			if(name == "class")
				cc = c;
		}
		if(c==xc)
			xChange(xAct);
		if(c==yc)
			yChange(yAct);
		if(c==cc)
			colorChange(cAct);
	}
}

void PlotWindow::modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
	updateOptionMenus();
}

void PlotWindow::startSVM()
{
	if(!svmWidget)
		svmWidget = new LibSVMWidget( scatter->model() );

	svmWidget->show();
}

/*
void PlotWindow::keyPressEvent(QKeyEvent *event)
 {	
     //switch (event->key()) {
	 //case Qt::Key_D:	//For delete
	//	 resultModel->deleteTrigger();
	//	 break;
     //default:
     //    QWidget::keyPressEvent(event);
     //}
	 
	 QWidget::keyPressEvent(event);
 }
*/
