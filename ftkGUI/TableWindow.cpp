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
	this->modAdapter = NULL;
	this->setup();
	this->rowsSelected = false;
}

//Destructor
TableWindow::~TableWindow()
{
	if(this->modAdapter)
    {
		delete this->modAdapter;
		this->modAdapter = NULL;
    }
}

void TableWindow::setQtModels(QItemSelectionModel * mod)
{
	this->tableView->setModel( (QAbstractItemModel*)mod->model() );
	this->tableView->setSelectionModel(mod);
	this->tableView->setSelectionBehavior( QAbstractItemView::SelectRows);
	this->update();
}
	
void TableWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)///////////////////////////////////////////
{
	if(this->modAdapter)
	{
		delete this->modAdapter;	
		this->modAdapter = NULL;
	}
	this->modAdapter = new vtkQtTableModelAdapter();
	this->modAdapter->setTable(  table );
	//this->modAdapter->SetKeyColumn(0);	//Key column is used as the row headers
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

	if(sels2)///////////////////////////////////////////////////////////
	{
		this->selection2 = sels2;
		connect(selection2, SIGNAL(changed()), this, SLOT(selectColumns()));
	}
		
	//Resize Rows to be as small as possible
	for (int i=0; i<tableView->model()->rowCount(); i++)
	{
		tableView->verticalHeader()->resizeSection(i,rowHeight);
	}

	Pointer2Table = table;
}

void TableWindow::setup()
{
	this->tableView = new QTableView();
	this->tableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
	this->setCentralWidget(this->tableView);
	this->setWindowTitle(tr("Features Table"));
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
	exportMenu = menuBar()->addMenu(tr("Export"));

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

	exportAction = new QAction(tr("Save As..."), this);
	connect(exportAction, SIGNAL(triggered()), this, SLOT(exportTable()));
	exportMenu->addAction(exportAction);

	exportIDAction = new QAction(tr("Export selected IDs"), this);
	connect(exportIDAction, SIGNAL(triggered()), this, SLOT(exportSelectedIDs()));
	exportMenu->addAction(exportIDAction);
}

// Exports the features data table to file
void TableWindow::exportTable()
{
	if(!this->tableView->model())
		return;

	QString Filename = QFileDialog::getSaveFileName( this, tr("Save Data Table"), QDir::currentPath(), tr("Microsoft Excel Document (*.xls);; Plain Text Document (*.txt)") );
    
	ofstream TableOutputStream;
	QFile DumpedTable(Filename);
	TableOutputStream.open(DumpedTable.fileName().toAscii());

	for( int i=0; i < this->tableView->model()->columnCount(); ++i)
	{	
		QString QHeaders = this->tableView->model()->headerData(i,Qt::Horizontal).toString();
		QHeaders.replace(" ", "_");
		std::string SHeaders = QHeaders.toStdString();
		TableOutputStream << SHeaders.c_str() << "\t";
	}

	TableOutputStream << "\n" ;
	
	if(rowsSelected == false)
	{
		for(vtkIdType row = 0; row < Pointer2Table->GetNumberOfRows(); row++ )
		{
			for(vtkIdType col = 0; col < Pointer2Table->GetNumberOfColumns(); col++ )
			{	
				TableOutputStream << Pointer2Table->GetValue(row,col).ToString() << "\t";
			}
			TableOutputStream << "\n";
		}
	}

	DumpedTable.close();
}

void TableWindow::exportSelectedIDs()
{
	if(!this->tableView->model())
		return;

	QString Filename = QFileDialog::getSaveFileName( this, tr("Save Selected IDs"), QDir::currentPath(), tr("Plain Text Document (*.txt)") );
    
	ofstream TableOutputStream(Filename.toStdString().c_str());
	std::set< long int> selectedIDs = this->selection->getSelections();
	std::set< long int>::iterator iter;
	for( iter = selectedIDs.begin(); iter != selectedIDs.end(); iter++)
	{
		TableOutputStream<< *iter <<"\t";
	}
	TableOutputStream.close();
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

//Sort does not work so well with the change to vtkTable. Removed for now!
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
	for( int i=1; i < this->tableView->model()->columnCount(); ++i) //skip root trace
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

	QString qstr1("Trace File"); //heatmap cannot take chars so skip
	QString qstr2("Distance to Device");
	int num_visible = 0;
	for (int i = 0; i < visible.size(); i++)
	{
		int equal1 = QString::compare(features[i],qstr1);
		int equal2 = QString::compare(features[i],qstr2);
		if ( equal1 != 0 && equal2 != 0)
			num_visible += visible.at(i);
	}

	if(num_visible == 0)
	{
		std::cout<< "No features were selected" << std::endl;
		return;
	}

	for( int i=1; i < this->tableView->model()->columnCount(); ++i)
	{	
		this->tableView->setColumnHidden( i, !visible.at(i-1) );
	}

	if(this->selection2)
	{
		std::set<long int> selectedIDs;
		
		selectedIDs.insert(0); //always show root trace

		for( int i=1; i < this->tableView->model()->columnCount(); ++i)
		{
			int equal1 = QString::compare(features[i-1],qstr1);
			int equal2 = QString::compare(features[i-1],qstr2);
			if ( equal1 == 0 || equal2 == 0)
				continue;
			if(visible.at(i-1))
				selectedIDs.insert(i);
		}

		this->selection2->select(selectedIDs);
	}
}

void TableWindow::selectColumns()//////////////////////////////////////////////////////
{
	std::set<long int> selectedIDs = this->selection2->getSelections();

	for(int i = 0; i < this->tableView->model()->columnCount(); i++)   // always showing the index and distance 
	{
		if (selectedIDs.find(i) != selectedIDs.end())
		{
			this->tableView->setColumnHidden( i+1, false);
		}
		else
		{
			this->tableView->setColumnHidden( i+1, true);
		}
	}
	this->tableView->setColumnHidden( 0, false);
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

	/*
	int y = tableView->rowViewportPosition(200);
	std::cerr << "Row 200 Viewport position is: " << y << std::endl;

	QRegion visRegion = tableView->visibleRegion();
	std::cerr << "top = " << visRegion.boundingRect().top() << " bottom = " << visRegion.boundingRect().bottom() << std::endl;

	//tableView->scroll(10,10); //This does not have desired effect.
	tableView->scrollTo( tableView->model()->index(200, 0) ); //This works great!!
	*/

	//tableView->verticalHeader()->hide();
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

	buttonSignalMapper = new QSignalMapper(this);
	featurSignalMapper = new QSignalMapper(this);
	for(int i=0; i<tests; ++i)
	{
		minVal[i] = new QDoubleSpinBox;
		maxVal[i] = new QDoubleSpinBox;
		minComp[i] = this->NewCompButton(i);
		maxComp[i] = this->NewCompButton(i);
		feature[i] = this->NewFeatureCombo(i);
	}
	connect(buttonSignalMapper, SIGNAL(mapped(int)), this, SLOT(ToggleComp(int)));
	connect(featurSignalMapper, SIGNAL(mapped(int)), this, SLOT(SetRanges(int)));

	this->InitRanges();

	for(int i=0; i<tests-1; ++i)
	{
		bools[i] = this->NewBoolCombo();
	}

	fLayout = new QGridLayout;

	numEquations = 0;
	AddEquation();

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
	for(int i=0; i<tests; ++i)
	{
		SetRanges(i);
	}
}

void FilterRowsDialog::SetRanges(int i)
{
	int c = this->GetColumnFor( feature[i]->currentText() );
	double max, min;
	this->GetMinMaxFor(c, &min, &max);
	minVal[i]->setRange(min,max);
	minVal[i]->setValue(min);
	minVal[i]->setSingleStep(.01);
	maxVal[i]->setRange(min,max);
	maxVal[i]->setValue(max);
	maxVal[i]->setSingleStep(.01);
}

QPushButton * FilterRowsDialog::NewCompButton(int n)
{
	QPushButton *button = new QPushButton(smaller);
	button->setDefault(false);
	button->setAutoDefault(false);
	connect(button, SIGNAL(clicked()), buttonSignalMapper, SLOT(map()));
	buttonSignalMapper->setMapping( button, n );
	return button;
}

void FilterRowsDialog::ToggleComp(int i)
{
	if( minComp[i]->text() == smaller )
	{
		minComp[i]->setText(bigger);
		maxComp[i]->setText(bigger);
	}
	else
	{
		minComp[i]->setText(smaller);
		maxComp[i]->setText(smaller);
	}
}

QComboBox * FilterRowsDialog::NewFeatureCombo(int n)
{
	QStringList features = GetVisibleFeatures();

	QComboBox *combo = new QComboBox;
	combo->addItems(features);
	connect(combo, SIGNAL(currentIndexChanged(QString)), featurSignalMapper, SLOT(map()));
	featurSignalMapper->setMapping( combo, n );
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
	if(numEquations+1 > tests)
		return;

	numEquations++;

	if(numEquations > 1)
		this->RemoveWidget(addButton);
	if(numEquations > 2)
		this->RemoveWidget(delButton);
	if(numEquations > 1)
		this->AddWidget(bools[numEquations-2], numEquations-2, 5);

	this->AddWidget(minVal[numEquations-1],numEquations-1,0);
	this->AddWidget(minComp[numEquations-1],numEquations-1,1);
	this->AddWidget(feature[numEquations-1],numEquations-1,2);
	this->AddWidget(maxComp[numEquations-1],numEquations-1,3);
	this->AddWidget(maxVal[numEquations-1],numEquations-1,4);

	if(numEquations > 1)
		this->AddWidget(delButton,numEquations-1,5);
	if(numEquations < tests)
		this->AddWidget(addButton,numEquations,0);
}

void FilterRowsDialog::RemoveEquation()
{
	if(numEquations-1 == 0)
		return;

	numEquations--;

	if(numEquations < tests-1)
		this->RemoveWidget(addButton);

	this->RemoveWidget(delButton);

	this->RemoveWidget(minVal[numEquations]);
	this->RemoveWidget(minComp[numEquations]);
	this->RemoveWidget(feature[numEquations]);
	this->RemoveWidget(maxComp[numEquations]);
	this->RemoveWidget(maxVal[numEquations]);

	this->RemoveWidget(bools[numEquations-1]);

	if(numEquations > 0)
		this->AddWidget(delButton,numEquations-1, 5);
	
	this->AddWidget(addButton, numEquations, 0);
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
	for(int i=0; i<numEquations; ++i)
	{
		featureColumns[i] = GetColumnFor( feature[i]->currentText() );
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

			if( minComp[c]->text() == smaller )	//I want to be inside the range
			{
				if( val >= minVal[c]->value() && val <= maxVal[c]->value() )
					ok[c] = true;
				else
					ok[c] = false;
			}
			else		// I want to be outside the range
			{
				if( val < minVal[c]->value() || val > maxVal[c]->value() )
					ok[c] = true;
				else
					ok[c] = false;
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
			if(bools[0]->currentText() == tr("AND"))
				mOK = ok[0] && ok[1];
			else if( bools[0]->currentText() == tr("OR"))
				mOK = ok[0] || ok[1];
		}
		else if( numEquations == 3)
		{
			if(bools[0]->currentText() == tr("AND") && bools[1]->currentText() == tr("AND"))
				mOK = ok[0] && ok[1] && ok[2];
			else if(bools[0]->currentText() == tr("AND") && bools[1]->currentText() == tr("OR"))
				mOK = (ok[0] && ok[1]) || ok[2];
			else if(bools[0]->currentText() == tr("OR") && bools[1]->currentText() == tr("AND"))
				mOK = (ok[0] || ok[1]) && ok[2];
			else if(bools[0]->currentText() == tr("OR") && bools[1]->currentText() == tr("OR"))
				mOK = ok[0] || ok[1] || ok[2];
		}

		//this->mTable->setRowHidden(row, !mOK);	//Removed because selection will select hidden rows!
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
	emit newRowsSelected();
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
	
	QModelIndexList sels = m_qmod->selectedRows();	//Gets Indexes for column=0
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
	emit newRowsSelected();
}
