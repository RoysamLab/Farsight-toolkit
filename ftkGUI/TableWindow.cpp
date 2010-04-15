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

//Constructors:
TableWindow::TableWindow(QWidget * parent)
: QMainWindow(parent)
{
	this->setup();
}

void TableWindow::setQtModels(QItemSelectionModel * mod)
{
	this->tableView->setModel( (QAbstractItemModel*)mod->model() );
	this->tableView->setSelectionModel(mod);
	this->tableView->setSelectionBehavior( QAbstractItemView::SelectRows );
	this->update();
}
	
void TableWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels)
{
	if(this->modAdapter)
		delete this->modAdapter;	
	this->modAdapter = new vtkQtTableModelAdapter();
	this->modAdapter->setTable(  table );
	this->tableView->setModel( modAdapter );

	QItemSelectionModel *mod = new QItemSelectionModel(modAdapter);
	this->tableView->setSelectionModel(mod);
	this->tableView->setSelectionBehavior( QAbstractItemView::SelectRows );

	if(sels)
	{
		this->selection = sels;
		if(selAdapter) delete selAdapter;
		selAdapter = new SelectionAdapter( this->tableView );
		selAdapter->SetPair(selection,mod);
	}
		
	//Resize Rows to be as small as possible
	for (int i=0; i<tableView->model()->rowCount(); i++)
	{
		tableView->verticalHeader()->resizeSection(i,rowHeight);
	}
}

void TableWindow::setup()
{
	this->tableView = new QTableView();
	this->tableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
	this->setCentralWidget(this->tableView);
	this->setWindowTitle(tr("Table"));
	this->createMenus();
	this->resize(800,300);

	selAdapter = NULL;
	modAdapter = NULL;
	// The following causes the program to get crashed if we close
	// the table first and the main window afterwards
	//this->setAttribute ( Qt::WA_DeleteOnClose );
}

void TableWindow::createMenus()
{
	viewMenu = menuBar()->addMenu(tr("View"));

	//sortByAction = new QAction(tr("Sort by..."),this);
	//connect(sortByAction, SIGNAL(triggered()), this, SLOT(sortBy()));
	//viewMenu->addAction(sortByAction);

	visibleColumnsAction = new QAction(tr("Visible Columns..."), this);
	connect(visibleColumnsAction, SIGNAL(triggered()), this, SLOT(changeColumns()));
	viewMenu->addAction(visibleColumnsAction);

	filterRowsAction = new QAction(tr("Filters..."), this);
	filterRowsAction->setStatusTip(tr("Hide/Show rows by applying filters"));
	connect(filterRowsAction, SIGNAL(triggered()), this, SLOT(showFilters()));
	viewMenu->addAction(filterRowsAction);

	testAction = new QAction(tr("Test"), this);
	testAction->setStatusTip(tr("Run Test Function"));
	connect(testAction, SIGNAL(triggered()), this, SLOT(test()));
	//viewMenu->addAction(testAction);

}

void TableWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
} 

void TableWindow::update()
{
	if(modAdapter)
		modAdapter->setTable( modAdapter->table() );

	//I took these out because sometimes they were really, really slow:
	//tableView->resizeRowsToContents();
	//tableView->resizeColumnsToContents();

	//Resize Rows to be as small as possible
	for (int i=0; i<tableView->model()->rowCount(); i++)
	{
		tableView->verticalHeader()->resizeSection(i,rowHeight);
	}
	QWidget::update();
}

void TableWindow::showFilters()
{
	if(!this->tableView->model())
		return;

	FilterRowsDialog * filters = new FilterRowsDialog( this->tableView, this->selection, this);
	filters->exec();
	delete filters;
}

void TableWindow::sortBy()
{
	if(!this->tableView->model())
		return;

	//Get the Currently Selected Features
	QStringList features;
	for( int i=0; i < this->tableView->model()->columnCount(); ++i)
	{	
		if( !this->tableView->isColumnHidden(i) )
		{
			features << this->tableView->model()->headerData(i,Qt::Horizontal).toString();
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
	for( int i=0; i < this->tableView->model()->columnCount(); ++i)
	{	
		if( this->tableView->model()->headerData(i,Qt::Horizontal).toString() == feat )
		{
			this->tableView->sortByColumn(i,Qt::AscendingOrder);
			emit sorted();
			break;
		}
	}
}

//Pop-up a window that allows the used to change the columns that are visible in the table
void TableWindow::changeColumns()
{
	if( !this->tableView->model() )
		return;

	//Get the Currently Visible Features in model:
	QStringList features;
	QList<bool> visible;
	for( int i=0; i < this->tableView->model()->columnCount(); ++i)
	{	
		features << this->tableView->model()->headerData(i,Qt::Horizontal).toString();
		visible << !(this->tableView->isColumnHidden(i));
	}
	//Let user choose one's to display using popup:
	ChooseItemsDialog *dialog = new ChooseItemsDialog(features, &visible, this);
	if(!dialog->exec())
	{
		delete dialog;
		return;
	}
	delete dialog;

	for( int i=0; i < this->tableView->model()->columnCount(); ++i)
	{	
		this->tableView->setColumnHidden( i, !visible.at(i) );
	}
	//this->update();
}


void TableWindow::test()
{
	//Are rows hidden?
	/*
	for (int i=0; i<tableView->model()->rowCount(); i++)
	{
		bool hidden = tableView->isRowHidden(i);
		std::cerr << i << " " << hidden << std::endl;
	}
	*/

	int y = tableView->rowViewportPosition(200);
	std::cerr << "Row 200 Viewport position is: " << y << std::endl;

	QRegion visRegion = tableView->visibleRegion();
	std::cerr << "top = " << visRegion.boundingRect().top() << " bottom = " << visRegion.boundingRect().bottom() << std::endl;

	//tableView->scroll(10,10); //This does not have desired effect.
	tableView->scrollTo( tableView->model()->index(200, 0) ); //This works great!!


}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// OTHER USEFUL CLASSES:
//***************************************************************************************************************
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


FilterRowsDialog::FilterRowsDialog(QTableView *tableView, ObjectSelection * sel, QWidget *parent)
: QDialog(parent)
{
	mTable = tableView;
	selection = sel;

	smaller = tr("<=");
	bigger = tr(">");

	//Add Filter button on left always:
	addButton = new QPushButton(tr("Add Filter..."));
	addButton->setDefault(false);
	addButton->setAutoDefault(false);
	connect(addButton,SIGNAL(clicked()), this, SLOT(AddEquation()));
	delButton = new QPushButton(tr("Remove"));
	delButton->setDefault(false);
	delButton->setAutoDefault(false);
	connect(delButton,SIGNAL(clicked()), this, SLOT(RemoveEquation()));

	minVal1 = new QDoubleSpinBox;
	minVal2 = new QDoubleSpinBox;
	minVal3 = new QDoubleSpinBox;
	maxVal1 = new QDoubleSpinBox;
	maxVal2 = new QDoubleSpinBox;
	maxVal3 = new QDoubleSpinBox;
	minComp1 = this->NewCompButton(1);
	minComp2 = this->NewCompButton(2);
	minComp3 = this->NewCompButton(3);
	maxComp1 = this->NewCompButton(1);
	maxComp2 = this->NewCompButton(2);
	maxComp3 = this->NewCompButton(3);
	feature1 = this->NewFeatureCombo();
	connect(feature1, SIGNAL(currentIndexChanged(QString)), this, SLOT(SetF1Ranges(QString)));
	feature2 = this->NewFeatureCombo();
	connect(feature2, SIGNAL(currentIndexChanged(QString)), this, SLOT(SetF2Ranges(QString)));
	feature3 = this->NewFeatureCombo();
	connect(feature3, SIGNAL(currentIndexChanged(QString)), this, SLOT(SetF3Ranges(QString)));
	this->InitRanges();
	bool1 = this->NewBoolCombo();
	bool2 = this->NewBoolCombo();

	fLayout = new QGridLayout;
	fLayout->addWidget(minVal1,0,0);
	fLayout->addWidget(minComp1,0,1);
	fLayout->addWidget(feature1,0,2);
	fLayout->addWidget(maxComp1,0,3);
	fLayout->addWidget(maxVal1,0,4);
	fLayout->addWidget(addButton,1,0);

	numEquations = 1;

	groupBox = new QGroupBox(tr("Create Filters"));
	groupBox->setLayout(fLayout);

	//Hide Button always on the right:
	updateButton = new QPushButton(tr("UPDATE"));
	connect(updateButton, SIGNAL(clicked()), this, SLOT(DoFilter()));
	QHBoxLayout *updateLayout = new QHBoxLayout;
	updateLayout->addStretch(50);
	updateLayout->addWidget(updateButton);
	addButton->setDefault(true);

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
	SetF1Ranges(feature1->currentText());
	minVal1->setSingleStep(.01);
	maxVal1->setSingleStep(.01);
	SetF2Ranges(feature2->currentText());
	minVal2->setSingleStep(.01);
	maxVal2->setSingleStep(.01);
	SetF3Ranges(feature3->currentText());
	minVal3->setSingleStep(.01);
	maxVal3->setSingleStep(.01);
}

void FilterRowsDialog::SetF1Ranges(QString text)
{
	int c = this->GetColumnFor(text);
	double max, min;
	this->GetMinMaxFor(c, &min, &max);
	minVal1->setRange(min,max);
	minVal1->setValue(min);
	maxVal1->setRange(min,max);
	maxVal1->setValue(max);
}

void FilterRowsDialog::SetF2Ranges(QString text)
{
	int c = this->GetColumnFor(text);
	double max, min;
	this->GetMinMaxFor(c, &min, &max);
	minVal2->setRange(min,max);
	minVal2->setValue(min);
	maxVal2->setRange(min,max);
	maxVal2->setValue(max);
}

void FilterRowsDialog::SetF3Ranges(QString text)
{
	int c = this->GetColumnFor(text);
	double max, min;
	this->GetMinMaxFor(c, &min, &max);
	minVal3->setRange(min,max);
	minVal3->setValue(min);
	maxVal3->setRange(min,max);
	maxVal3->setValue(max);
}

QPushButton * FilterRowsDialog::NewCompButton(int n)
{
	QPushButton *button = new QPushButton(smaller);
	button->setDefault(false);
	button->setAutoDefault(false);
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

int FilterRowsDialog::GetColumnFor(QString headerText)
{
	for( int i=0; i < this->mTable->model()->columnCount(); ++i)
	{	
		QString f = this->mTable->model()->headerData(i,Qt::Horizontal).toString();
		if( f == headerText )
		{
			return i;
		}
	}
	return -1;
}

void FilterRowsDialog::GetMinMaxFor(int c, double *min, double *max)
{
	//Initialize
	QModelIndex index = this->mTable->model()->index(0, c);
	double val = this->mTable->model()->data(index).toDouble();
	(*min) = val;
	(*max) = val;

	for( int row=1; row < this->mTable->model()->rowCount(); ++row)
	{
		QModelIndex index = this->mTable->model()->index(row, c);
		double val = this->mTable->model()->data(index).toDouble();
		if( val > (*max) ) (*max) = val;
		else if( val < (*min) ) (*min) = val;
	}
}

void FilterRowsDialog::DoFilter(void)
{
	//First find out which column numbers I care about!
	int featureColumns[3];
	switch(numEquations)
	{
	case 3:
		featureColumns[2] = GetColumnFor( feature3->currentText() );
	case 2:
		featureColumns[1] = GetColumnFor( feature2->currentText() );
	case 1:
		featureColumns[0] = GetColumnFor( feature1->currentText() );
		break;
	}

	std::set<long int> matchIds;	//Ids that match the filter

	//************************************************************************************
	// Do each equation
	// NOTE: smaller means in between two values (including ends)
	//       bigger means outside of two values
	//*************************************************************************************
	for( int row=0; row < this->mTable->model()->rowCount(); ++row)
	{
		bool ok[3] = {false, false, false};
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

		bool mOK = false;				//Master OK

		//Now check the equations:
		if(numEquations == 1)
		{
			mOK = ok[0];
		}
		else if( numEquations == 2)
		{
			if(bool1->currentText() == tr("AND"))
				mOK = ok[0] && ok[1];
			else if( bool1->currentText() == tr("OR"))
				mOK = ok[0] || ok[1];
		}
		else if( numEquations == 3)
		{
			if(bool1->currentText() == tr("AND") && bool2->currentText() == tr("AND"))
				mOK = ok[0] && ok[1] && ok[2];
			else if(bool1->currentText() == tr("AND") && bool2->currentText() == tr("OR"))
				mOK = (ok[0] && ok[1]) || ok[2];
			else if(bool1->currentText() == tr("OR") && bool2->currentText() == tr("AND"))
				mOK = (ok[0] || ok[1]) && ok[2];
			else if(bool1->currentText() == tr("OR") && bool2->currentText() == tr("OR"))
				mOK = ok[0] || ok[1] || ok[2];
		}

		//this->mTable->setRowHidden(row, !mOK);
		if(mOK)
		{
			QModelIndex index = this->mTable->model()->index(row, 0);
			int val = this->mTable->model()->data(index).toInt();
			matchIds.insert(val);
		}
	}

	selection->select(matchIds);
	accept();
}

//******************************************************************************************************
//******************************************************************************************************
//******************************************************************************************************
// This class will convert between the two selection models (make sure they are syncronized)
//******************************************************************************************************
SelectionAdapter::SelectionAdapter()
: QObject()
{
	m_obj = NULL;
	m_qmod = NULL;
	okToChange = TRUE;
	m_tableView = NULL;
}

SelectionAdapter::SelectionAdapter( QTableView * tableView)
{
	m_obj = NULL;
	m_qmod = NULL;
	okToChange = TRUE;
	m_tableView = tableView;
}

void SelectionAdapter::SetPair(ObjectSelection * obj, QItemSelectionModel * qmod)
{
	m_obj = obj;
	m_qmod = qmod;

	connect(m_obj, SIGNAL(changed()), this, SLOT(updateQMOD()));
	connect(m_qmod, SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)), \
			this,SLOT(updateOBJ(const QItemSelection &, const QItemSelection &)));
}

//This class is when an object gets selected from outside the table - so table needs to reflect the selection
void SelectionAdapter::updateQMOD(void)
{
	if(!okToChange) return;
	okToChange = FALSE;
	const QAbstractItemModel * model = m_qmod->model();
	int rows = model->rowCount();
	QItemSelection selection;
	selection.clear();

	for (int row = 0; row < rows; ++row) 
	{
		QModelIndex index1 = model->index(row, 0);
		QModelIndex index2 = model->index(row, model->columnCount()-1);
		int id = model->data(index1).toInt();
		if (m_obj->isSelected(id))
		{
			selection.merge(QItemSelection(index1,index2),QItemSelectionModel::Select);
		}
    }
	m_qmod->select(selection, QItemSelectionModel::ClearAndSelect);

	//Make sure first selection is visible
	if( selection.size() > 0 && m_tableView )
	{
		int x = m_tableView->horizontalScrollBar()->value();
		m_tableView->scrollTo( selection.indexes().at(0) ); //This works great!!
		m_tableView->horizontalScrollBar()->setValue(x);
	}

	okToChange = TRUE;
}

//The class gets called when a row in the table is clicked/selected
void SelectionAdapter::updateOBJ(const QItemSelection & selected, const QItemSelection & deselected)
{
	if(!okToChange) return;
	okToChange = FALSE;
	std::set<long int> ids;
	const QAbstractItemModel * model = m_qmod->model();
	//int rows = model->rowCount();
	//int columns = model->columnCount();
	
	QModelIndexList sels = m_qmod->selectedRows();
	for(int i=0; i<sels.size(); ++i)
	{
		QModelIndex index = sels.at(i);
		long int id = model->data(index).toInt();
		ids.insert(id);
	}
	//for (int row=0; row<rows; ++row)
	//{
	//	QModelIndex index = model->index(row, 0);
	//	long int id = model->data(index).toInt();
	//	if (m_qmod->isSelected(index))
	//	{
	//		ids.insert(id);
	//	}
	//}
	m_obj->select(ids);
	okToChange = TRUE;
}
