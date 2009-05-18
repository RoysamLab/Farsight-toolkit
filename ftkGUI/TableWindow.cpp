//***********************************************************************************
// TableWindow provides a classic table view into a model
//***********************************************************************************
#include "TableWindow.h"

//Constructor
TableWindow::TableWindow(QItemSelectionModel *selectionModel, QWidget *parent)
: QWidget(parent)
{
	this->table = new QTableView();
	this->table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	this->layout = new QVBoxLayout();
	this->layout->addWidget(table);
	this->layout->setContentsMargins(2,2,2,2); 
	this->setLayout(layout);
	this->setWindowTitle(tr("Table"));
	this->setAttribute ( Qt::WA_DeleteOnClose );

	this->table->setModel( (QAbstractItemModel*)selectionModel->model() );
	this->table->setSelectionModel(selectionModel);
	
	this->table->setSelectionBehavior( QAbstractItemView::SelectRows );

	//QWidget::connect(mod, SIGNAL(modelChanged()), this, SLOT(update()));
}

/*
void TableWindow::setup()
{
	table = new QTableView();
	table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	layout = new QVBoxLayout();
	layout->addWidget(table);
	layout->setContentsMargins(2,2,2,2); 
	setLayout(layout);
	setWindowTitle(tr("Table"));
	setAttribute ( Qt::WA_DeleteOnClose );
}
*/
/*
TableWindow::TableWindow(QWidget *parent) 
 : QWidget(parent) 
{
	this->setup();
	//QFont f = font();
	//f.setPointSize(8);
	//table->setFont(f);
}
*/

/*
TableWindow::TableWindow(SegmentationModel *mod, QWidget *parent) 
 : QWidget(parent)
{
	this->setup();

	table->setModel( mod->GetModel() );
	table->setSelectionModel( mod->GetSelectionModel() );
	connect(mod, SIGNAL(modelChanged()), this, SLOT(update()));
	visibleRows = mod->NumFeatures()+2;
	//for(int i = visibleRows; i < table->model()->rowCount(); ++i)
	//	table->setColumnHidden(i,true);
}
*/

void TableWindow::update()
{
	QWidget::update();
	table->resizeRowsToContents();
	table->resizeColumnsToContents();
	//Resize Rows to be as small as possible
	for (int i=0; i<table->model()->rowCount(); i++)
	{
		table->verticalHeader()->resizeSection(i,18);
	}
	//for(int i = visibleRows; i < table->model()->rowCount(); ++i)
	//{
	//	table->setColumnHidden(i,true);
	//}
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
	if (numRows > 5) numRows = 5;
	int bestHeight =( numRows + 1 ) * rowHeight;

	int bestWidth = 0;
	for (int i=0; i<table->model()->columnCount(); i++)
	{
		bestWidth = bestWidth + table->columnWidth(i);
	}
	bestWidth = bestWidth + 100;

	resize(bestWidth,bestHeight+5);

	if (this->frameGeometry().width() > screenWidth)
		resize(screenWidth-10,bestHeight+5);
}
