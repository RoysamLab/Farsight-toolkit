//****************************************************************************************
// Plot window uses scatterview to display a scatterplot of the data in a model.
// In order for this view to work the model must be set.
// The window allows the user to change the columns used in the scatterplot.
//****************************************************************************************
#include "PlotWindow.h"

//Constructor
PlotWindow::PlotWindow(QWidget *parent)
  : QWidget(parent)
{
	this->setupUI();
	resultModel = NULL;
}

PlotWindow::PlotWindow(SegmentationModel *rModel,QWidget *parent)
  : QWidget(parent)
{
	this->setupUI();

	resultModel = rModel;
	connect(resultModel, SIGNAL(modelChanged()), scatter, SLOT(dataChanged()));
	connect(resultModel, SIGNAL(modelChanged()), this, SLOT(updateColumnForColor()));

	scatter->setModel( resultModel->GetModel() );
	scatter->setSelectionModel( resultModel->GetSelectionModel() );
	scatter->SetColForColor(resultModel->ColumnForColor(), resultModel->ColorMap());
	updateCombos();
}

void PlotWindow::setupUI(void)
{
	resize(500, 500);

	comboX = new QComboBox();
	comboY = new QComboBox();
	comboSelMode = new QComboBox();
    scatter = new ScatterView();
	vlayout = new QVBoxLayout();
	hlayout = new QHBoxLayout();
	hlayoutT = new QHBoxLayout();

	selectButton = new QPushButton("Select", this);
	clearButton = new QPushButton("Clear",this);

	ylabel = new QLabel("y: ");
	xlabel = new QLabel("x: ");
	colorlabel = new QLabel("Color by: ");

	hlayout->addStretch(10);
	hlayout->addWidget(ylabel);
	hlayout->addWidget(comboY);
	hlayout->addStretch(50);
	hlayout->addWidget(xlabel);
	hlayout->addWidget(comboX);
	hlayout->addStretch(50);

	hlayoutT->addWidget(comboSelMode);
	hlayoutT->addStretch(50);
	hlayoutT->addWidget(selectButton);
	hlayoutT->addWidget(clearButton);

	vlayout->addLayout(hlayoutT);
	vlayout->addWidget(scatter);
	vlayout->addLayout(hlayout);
	setLayout(vlayout);

	comboX->setCurrentIndex(1);
	comboY->setCurrentIndex(2);
	connect(comboX, SIGNAL(currentIndexChanged(int)),this,SLOT(comboXChange(int)));
	connect(comboY, SIGNAL(currentIndexChanged(int)),this,SLOT(comboYChange(int)));
	comboX->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	comboY->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	
	setupSelectionModes();
	connect(comboSelMode, SIGNAL(currentIndexChanged(int)), scatter, SLOT(selModeChanged(int)));
	connect(selectButton, SIGNAL(clicked()), scatter, SLOT(selectClicked()));
	connect(clearButton, SIGNAL(clicked()), scatter, SLOT(clearClicked()));

	setWindowTitle(tr("Scatter Plot"));
	setAttribute ( Qt::WA_DeleteOnClose );
}

void PlotWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
}

void PlotWindow::keyPressEvent(QKeyEvent *event)
 {
     switch (event->key()) {
	 case Qt::Key_D:	//For delete
		 resultModel->deleteTrigger();
		 break;
     default:
         QWidget::keyPressEvent(event);
     }
 }

void PlotWindow::updateColumnForColor()
{
	scatter->SetColForColor(resultModel->ColumnForColor(), resultModel->ColorMap());
}

void PlotWindow::comboXChange(int c)
{
	scatter->SetColForX(c+1);
}

void PlotWindow::comboYChange(int c)
{
	scatter->SetColForY(c+1);
}

void PlotWindow::updateCombos()
{
	QStandardItemModel *model = resultModel->GetModel();

	//Setup x/y combos
	comboX->clear();
	comboY->clear();
	for (int c = 1; c <= resultModel->NumFeatures(); ++c )
	{
		comboX->addItem(model->headerData(c,Qt::Horizontal).toString());
		comboY->addItem(model->headerData(c,Qt::Horizontal).toString());
	}
	scatter->SetColForX(1);
	scatter->SetColForY(2);
	comboX->setCurrentIndex(0);
	comboY->setCurrentIndex(1);
}

void PlotWindow::setupSelectionModes(void)
{
	comboSelMode->clear();
	comboSelMode->addItem("SingleSelectionMode");
	comboSelMode->addItem("RegionSelectionMode");
	comboSelMode->setCurrentIndex(0);
}
