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

//***********************************************************************************
// TableWindow provides a classic table view into a model
//***********************************************************************************
#include "TableWindow.h"

//Constructor
TableWindow::TableWindow(QItemSelectionModel *selectionModel, QWidget *parent)
: QMainWindow(parent)
{
	this->table = new QTableView();
	this->table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	this->setCentralWidget(this->table);
	this->setWindowTitle(tr("Table"));
	// The following causes the program to get crashed if we close
	// the table first and the main window afterwards
	//this->setAttribute ( Qt::WA_DeleteOnClose );

	this->table->setModel( (QAbstractItemModel*)selectionModel->model() );
	this->table->setSelectionModel(selectionModel);
	this->table->setSelectionBehavior( QAbstractItemView::SelectRows );
	connect(selectionModel->model(), SIGNAL(dataChanged(const QModelIndex &, const QModelIndex &)), this, SLOT(modelChange(const QModelIndex &, const QModelIndex &)));

	this->createMenus();
}

void TableWindow::createMenus()
{
	viewMenu = menuBar()->addMenu(tr("View"));
	sortByAction = new QAction(tr("Sort by..."),this);
	connect(sortByAction, SIGNAL(triggered()), this, SLOT(sortBy()));
	viewMenu->addAction(sortByAction);

	visibleColumnsAction = new QAction(tr("Visible Columns..."), this);
	connect(visibleColumnsAction, SIGNAL(triggered()), this, SLOT(changeColumns()));
	viewMenu->addAction(visibleColumnsAction);
}

void TableWindow::sortBy()
{
	//Get the Currently Selected Features
	QStringList features;
	for( int i=0; i < this->table->model()->columnCount(); ++i)
	{	
		if( !this->table->isColumnHidden(i) )
		{
			features << this->table->model()->headerData(i,Qt::Horizontal).toString();
		}
	}
	
	//Let user choose one using popup:
	ChooseItemDialog *dialog = new ChooseItemDialog(features, this);
	if(!dialog->exec())
	{
		delete dialog;
		return;
	}
	QString feat = dialog->getSelectedItem();
	delete dialog;
	//Which column is this feature?
	for( int i=0; i < this->table->model()->columnCount(); ++i)
	{	
		if( this->table->model()->headerData(i,Qt::Horizontal).toString() == feat )
		{
			//this->table->sortByColumn(i,Qt::AscendingOrder);

			emit sorted();
			break;
		}
	}
}

//Pop-up a window that allows the used to chang the columns that are visible in the table
void TableWindow::changeColumns()
{
	//Get the Currently Features in model:
	QStringList features;
	QList<bool> visible;
	for( int i=0; i < this->table->model()->columnCount(); ++i)
	{	
		features << this->table->model()->headerData(i,Qt::Horizontal).toString();
		visible << !(this->table->isColumnHidden(i));
	}
	//Let user choose one's to display using popup:
	ChooseItemsDialog *dialog = new ChooseItemsDialog(features, &visible, this);
	if(!dialog->exec())
	{
		delete dialog;
		return;
	}
	delete dialog;

	for( int i=0; i < this->table->model()->columnCount(); ++i)
	{	
		this->table->setColumnHidden( i, !visible.at(i) );
	}
	this->ResizeToOptimalSize();
}

void TableWindow::modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
	table->update();
	int l_column = topLeft.column();
	int t_row = topLeft.row();
	int r_column = bottomRight.column();
	int b_row = bottomRight.row();

	for(int r=t_row; r<=b_row; r++)
	{
		table->resizeRowToContents(r);
		table->verticalHeader()->resizeSection(r,18);
	}
	for(int c=l_column; c<=r_column; c++)
	{
		table->resizeColumnToContents(c);
	}
	QWidget::update();
}

void TableWindow::update()
{
	table->update();
	table->resizeRowsToContents();
	table->resizeColumnsToContents();
	//Resize Rows to be as small as possible
	for (int i=0; i<table->model()->rowCount(); i++)
	{
		table->verticalHeader()->resizeSection(i,18);
	}
	QWidget::update();
}

void TableWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
} 

/*
void TableWindow::SetModels(QItemSelectionModel *selectionModel)
{
	table->setModel( (QAbstractItemModel*)selectionModel->model() );
	table->setSelectionModel(selectionModel);
}
*/

//***********************************************************************************
//This function calculates the optimal size of the window so that all columns can be
//displayed.  It also figures the a good height depending on the number of rows in the
//model.
//***********************************************************************************
void TableWindow::ResizeToOptimalSize(void)
{
	int screenWidth = qApp->desktop()->width();
	//Resize rows to minimum height
	table->resizeRowsToContents();
	table->resizeColumnsToContents();

	//Resize Rows to be as small as possible
	for (int i=0; i<table->model()->rowCount(); i++)
	{
		table->verticalHeader()->resizeSection(i,18);
	}

	//assumes all rows have the same height
	int rowHeight = table->rowHeight(0);
	int numRows = table->model()->rowCount();
	if (numRows > 15) numRows = 15;
	int bestHeight =( numRows + 1 ) * rowHeight;

	int bestWidth = 0;
	for (int i=0; i<table->model()->columnCount(); i++)
	{
		bestWidth = bestWidth + table->columnWidth(i);
		if(bestWidth > 600)
			break;
	}
	bestWidth = bestWidth + 50;

	resize(bestWidth,bestHeight+5);

	if (this->frameGeometry().width() > screenWidth)
		resize(screenWidth-10,bestHeight+5);
}

ChooseItemDialog::ChooseItemDialog(QStringList items, QWidget *parent)
: QDialog(parent)
{
	itemCombo = new QComboBox();
	for (int i = 0; i < items.size(); ++i)
	{
		itemCombo->addItem( items.at(i) );
	}

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	QHBoxLayout *bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addWidget(itemCombo);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Choose Item"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

ChooseItemsDialog::ChooseItemsDialog(QStringList items, QList<bool> * _selected, QWidget *parent)
: QDialog(parent)
{
	this->selected = _selected;

	//First set up the items choices:
	QGroupBox *groupBox = new QGroupBox();
	QVBoxLayout *itemLayout = new QVBoxLayout;

	QHBoxLayout *selLayout = new QHBoxLayout;
	QLabel *selLabel = new QLabel(tr("Select:"));
	QPushButton *allButton = new QPushButton(tr("All"));
	QPushButton *nonButton = new QPushButton(tr("None"));
	connect(allButton, SIGNAL(clicked()), this, SLOT(selectAll()));
	connect(nonButton, SIGNAL(clicked()), this, SLOT(selectNone()));
	selLayout->addWidget(selLabel);
	selLayout->addWidget(allButton);
	selLayout->addWidget(nonButton);
	selLayout->addStretch(10);
	itemLayout->addLayout(selLayout);

	QWidget *groupWidget = new QWidget;
	QVBoxLayout *vLayout = new QVBoxLayout;
	itemGroup = new QButtonGroup;
	itemGroup->setExclusive(false);

	for (int i = 0; i < items.size(); ++i)
	{
		QString name = items.at(i);
		QCheckBox *check = new QCheckBox(name);
		if(selected->at(i))
			check->setChecked(true);
		else
			check->setChecked(false);
		vLayout->addWidget(check);
		itemGroup->addButton(check, i);
	}
	connect(itemGroup, SIGNAL(buttonClicked(int)), this, SLOT(selectionChanged(int)));
	groupWidget->setLayout(vLayout);

	QScrollArea *scrollArea = new QScrollArea;
	scrollArea->setWidget(groupWidget);
	scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	itemLayout->addWidget(scrollArea);
	groupBox->setLayout(itemLayout);
	
	//The ok button:
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	QHBoxLayout *bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	//The final layout:
	QVBoxLayout *layout = new QVBoxLayout;
	layout->addWidget(groupBox);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Choose Items"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

void ChooseItemsDialog::selectNone()
{
	QList<QAbstractButton *> buttons = itemGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			buttons.at(b)->setChecked(false);
			selected->replace(b,false);
		}
	}
}

void ChooseItemsDialog::selectAll()
{
	QList<QAbstractButton *> buttons = itemGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( !buttons.at(b)->isChecked() )
		{
			buttons.at(b)->setChecked(true);
			selected->replace(b,true);
		}
	}
}

void ChooseItemsDialog::selectionChanged(int id)
{
	QList<QAbstractButton *> buttons = itemGroup->buttons();
	selected->replace( id, buttons.at(id)->isChecked() );
}