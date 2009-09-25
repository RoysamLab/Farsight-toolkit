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
	this->filters = new FilterRowsDialog( this->table, this);
	this->filters->hide();
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

	filterRowsAction = new QAction(tr("Filters..."), this);
	filterRowsAction->setStatusTip(tr("Hide/Show rows by applying filters"));
	connect(filterRowsAction, SIGNAL(triggered()), this, SLOT(showFilters()));
	viewMenu->addAction(filterRowsAction);
}

void TableWindow::showFilters()
{
	this->filters->show();
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

//Pop-up a window that allows the used to change the columns that are visible in the table
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


FilterRowsDialog::FilterRowsDialog(QTableView *table, QWidget *parent)
: QDialog(parent)
{
	this->setVisible(false);

	mTable = table;

	smaller = tr("<=");
	bigger = tr(">");

	//Add Filter button on left always:
	addButton = new QPushButton(tr("Add Filter..."));
	connect(addButton,SIGNAL(clicked()), this, SLOT(AddEquation()));
	delButton = new QPushButton(tr("Remove"));
	connect(delButton,SIGNAL(clicked()), this, SLOT(RemoveEquation()));

	minVal1 = new QDoubleSpinBox;
	minVal2 = new QDoubleSpinBox;
	minVal3 = new QDoubleSpinBox;
	maxVal1 = new QDoubleSpinBox;
	maxVal2 = new QDoubleSpinBox;
	maxVal3 = new QDoubleSpinBox;
	this->InitRanges();
	minComp1 = this->NewCompButton(1);
	minComp2 = this->NewCompButton(2);
	minComp3 = this->NewCompButton(3);
	maxComp1 = this->NewCompButton(1);
	maxComp2 = this->NewCompButton(2);
	maxComp3 = this->NewCompButton(3);
	feature1 = this->NewFeatureCombo();
	feature2 = this->NewFeatureCombo();
	feature3 = this->NewFeatureCombo();
	bool1 = this->NewBoolCombo();
	bool2 = this->NewBoolCombo();

	fLayout = new QGridLayout;
	this->AddWidget(minVal1,0,0);
	this->AddWidget(minComp1,0,1);
	this->AddWidget(feature1,0,2);
	this->AddWidget(maxComp1,0,3);
	this->AddWidget(maxVal1,0,4);
	this->AddWidget(addButton,1,0);

	numEquations = 1;

	groupBox = new QGroupBox(tr("Create Filters"));
	groupBox->setLayout(fLayout);

	//Hide Button always on the right:
	updateButton = new QPushButton(tr("UPDATE"));
	connect(updateButton, SIGNAL(clicked()), this, SLOT(DoFilter()));
	QHBoxLayout *updateLayout = new QHBoxLayout;
	updateLayout->addStretch(50);
	updateLayout->addWidget(updateButton);

	//topLayout manages the filters and the hide button:
	QVBoxLayout *topLayout = new QVBoxLayout();
	topLayout->addWidget(groupBox);
	topLayout->addStretch(50);
	topLayout->addLayout(updateLayout);
	this->setLayout(topLayout);
	this->setWindowTitle(tr("Filters"));
}

//Set initial ranges of the spin boxes from 0 to the highest ID!
void FilterRowsDialog::InitRanges()
{
	double max = 0;
	for( int row=0; row < this->mTable->model()->rowCount(); ++row)
	{
		QModelIndex index = this->mTable->model()->index(row, 0);
		double val = this->mTable->model()->data(index).toDouble();
		if( val > max ) max = val;
	}
	maxVal1->setMaximum(max);
	maxVal2->setMaximum(max);
	maxVal3->setMaximum(max);
}

QPushButton * FilterRowsDialog::NewCompButton(int n)
{
	QPushButton *button = new QPushButton(smaller);
	switch(n)
	{
	case 1:
		connect(button, SIGNAL(clicked()), this, SLOT(ToggleComp1()));
		break;
	case 2:
		connect(button, SIGNAL(clicked()), this, SLOT(ToggleComp2()));
		break;
	case 3:
		connect(button, SIGNAL(clicked()), this, SLOT(ToggleComp3()));
		break;
	}
	return button;
}

void FilterRowsDialog::ToggleComp(int n)
{
	switch(n)
	{
	case 1:
		if( minComp1->text() == smaller )
		{
			minComp1->setText(bigger);
			maxComp1->setText(bigger);
		}
		else
		{
			minComp1->setText(smaller);
			maxComp1->setText(smaller);
		}
		break;
	case 2:
		if( minComp2->text() == smaller )
		{
			minComp2->setText(bigger);
			maxComp2->setText(bigger);
		}
		else
		{
			minComp2->setText(smaller);
			maxComp2->setText(smaller);
		}
		break;
	case 3:
		if( minComp3->text() == smaller )
		{
			minComp3->setText(bigger);
			maxComp3->setText(bigger);
		}
		else
		{
			minComp3->setText(smaller);
			maxComp3->setText(smaller);
		}
		break;
	}
}

QComboBox * FilterRowsDialog::NewFeatureCombo()
{
	QStringList features = GetVisibleFeatures();
	QComboBox *combo = new QComboBox;
	combo->addItems(features);
	return combo;
}

QStringList FilterRowsDialog::GetVisibleFeatures()
{
	QStringList features;
	for( int i=0; i < this->mTable->model()->columnCount(); ++i)
	{
		if ( !(this->mTable->isColumnHidden(i)) )
		{
			features << this->mTable->model()->headerData(i,Qt::Horizontal).toString();
		}
	}
	return features;
}

QComboBox * FilterRowsDialog::NewBoolCombo()
{
	QComboBox *combo = new QComboBox;
	combo->addItem(tr("AND"));
	combo->addItem(tr("OR"));
	return combo;
}

void FilterRowsDialog::AddEquation()
{
	switch(numEquations)
	{
	case 1:
		numEquations = 2;
		this->RemoveWidget(addButton);
		this->AddWidget(bool1,0,5);
		this->AddWidget(minVal2,1,0);
		this->AddWidget(minComp2,1,1);
		this->AddWidget(feature2,1,2);
		this->AddWidget(maxComp2,1,3);
		this->AddWidget(maxVal2,1,4);
		this->AddWidget(delButton,1,5);
		this->AddWidget(addButton,2,0);	
		break;
	case 2:
		numEquations = 3;
		this->RemoveWidget(delButton);
		this->RemoveWidget(addButton);
		this->AddWidget(bool2,1,5);
		this->AddWidget(minVal3,2,0);
		this->AddWidget(minComp3,2,1);
		this->AddWidget(feature3,2,2);
		this->AddWidget(maxComp3,2,3);
		this->AddWidget(maxVal3,2,4);
		this->AddWidget(delButton,2,5);
		break;
	default:
		break;
	}
}

void FilterRowsDialog::RemoveEquation()
{
	switch(numEquations)
	{
	case 2:
		numEquations = 1;
		this->RemoveWidget(delButton);
		this->RemoveWidget(addButton);
		this->RemoveWidget(bool1);
		this->RemoveWidget(minVal2);
		this->RemoveWidget(minComp2);
		this->RemoveWidget(feature2);
		this->RemoveWidget(maxComp2);
		this->RemoveWidget(maxVal2);
		this->AddWidget(addButton,2,0);	
		break;
	case 3:
		numEquations = 2;
		this->RemoveWidget(delButton);
		this->RemoveWidget(bool2);
		this->RemoveWidget(minVal3);
		this->RemoveWidget(minComp3);
		this->RemoveWidget(feature3);
		this->RemoveWidget(maxComp3);
		this->RemoveWidget(maxVal3);
		this->AddWidget(delButton,1,5);
		this->AddWidget(addButton,2,0);
		break;
	default:
		break;
	}
}

void FilterRowsDialog::RemoveWidget(QWidget *widget)
{
	widget->setVisible(false);
	fLayout->removeWidget(widget);
}

void FilterRowsDialog::AddWidget(QWidget *widget, int r, int c)
{
	fLayout->addWidget(widget, r, c);
	widget->setVisible(true);
}

void FilterRowsDialog::DoFilter(void)
{
	//First find out which column numbers I care about!
	int featureColumns[3];
	for( int i=0; i < this->mTable->model()->columnCount(); ++i)
	{	
		QString f = this->mTable->model()->headerData(i,Qt::Horizontal).toString();
		if( f == feature1->currentText() )
		{
			featureColumns[0] = i;
		}
		else if ( f == feature2->currentText() )
		{
			featureColumns[1] = i;
		}
		else if ( f == feature2->currentText() )
		{
			featureColumns[2] = i;
		}
	}

	//************************************************************************************
	// Do each equation
	// NOTE: smaller means in between two values (including ends)
	//       bigger means outside of two values
	//*************************************************************************************
	for( int row=0; row < this->mTable->model()->rowCount(); ++row)
	{
		bool ok[3];
		for(int c=0; c<numEquations; c++)
		{
			QModelIndex index = this->mTable->model()->index(row, featureColumns[c]);
			double val = this->mTable->model()->data(index).toDouble();
			if(c==0)
			{
				if( minComp1->text() == smaller )	//I want to be inside the range
				{
					if( val >= minVal1->value() && val <= maxVal1->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
				else		// I want to be outside the range
				{
					if( val < minVal1->value() || val > maxVal1->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
			}
			else if(c==1)
			{
				if( minComp2->text() == smaller )	//I want to be inside the range
				{
					if( val >= minVal2->value() && val <= maxVal2->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
				else		// I want to be outside the range
				{
					if( val < minVal2->value() || val > maxVal2->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
			}
			else if(c==2)
			{
				if( minComp3->text() == smaller )	//I want to be inside the range
				{
					if( val >= minVal3->value() && val <= maxVal3->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
				else		// I want to be outside the range
				{
					if( val < minVal3->value() || val > maxVal3->value() )
						ok[c] = true;
					else
						ok[c] = false;
				}
			}
		}

		//Now check the equations:
		if(numEquations == 1)
		{
			this->mTable->setRowHidden(row, !ok[0]);
		}
		else if( numEquations == 2)
		{
			if(bool1->currentText() == tr("AND"))
				this->mTable->setRowHidden(row, !( ok[0] && ok[1] ) );
			else if( bool1->currentText() == tr("OR"))
				this->mTable->setRowHidden(row, !( ok[0] || ok[1] ) );

		}
		else if( numEquations == 3)
		{
			if(bool1->currentText() == tr("AND") && bool2->currentText() == tr("AND"))
				this->mTable->setRowHidden(row, !( ok[0] && ok[1] && ok[2] ) );
			else if(bool1->currentText() == tr("AND") && bool2->currentText() == tr("OR"))
				this->mTable->setRowHidden(row, !( (ok[0] && ok[1]) || ok[2] ) );
			else if(bool1->currentText() == tr("OR") && bool2->currentText() == tr("AND"))
				this->mTable->setRowHidden(row, !( (ok[0] || ok[1]) && ok[2] ) );
			else if(bool1->currentText() == tr("OR") && bool2->currentText() == tr("OR"))
				this->mTable->setRowHidden(row, !( (ok[0] || ok[1]) || ok[2] ) );
		}
	}
}