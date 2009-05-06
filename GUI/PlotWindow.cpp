//****************************************************************************************
// Plot window uses scatterview to display a scatterplot of the data in a model.
// In order for this view to work the model must be set.
// The window allows the user to change the columns used in the scatterplot.
//****************************************************************************************
#include "PlotWindow.h"

#include "libsvm/svm.h"

//for DBL_MAX...
#include <float.h>
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
	outlierButton = new QPushButton("Find Outliers",this);

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

	hlayoutT->addWidget(outlierButton);
	hlayoutT->addStretch(80);
	hlayoutT->addWidget(comboSelMode);
	hlayoutT->addStretch(20);
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
	connect(outlierButton, SIGNAL(clicked()), this, SLOT(findOutliers()));

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

void PlotWindow::findOutliers()
{
	//First create two vectors with the features:
	int c1 = comboX->currentIndex()+1;
	int c2 = comboY->currentIndex()+1;
	int cID = resultModel->ColumnForID();	//and one with ids

	std::vector<int> ids;
	std::vector<double> f1;
	std::vector<double> f2;
	
	//Init these to the extreme opposite ends of double
	double f1_max = -DBL_MAX;
	double f1_min = DBL_MAX;
	double f2_max = -DBL_MAX;
	double f2_min = DBL_MAX;

	QModelIndex index;
	for(int i=0; i<(int)resultModel->NumObjects(); ++i)
	{
		index = resultModel->GetModel()->index(i, cID);
		double id = resultModel->GetModel()->data(index).toDouble();
		index = resultModel->GetModel()->index(i, c1);
		double f_1 = resultModel->GetModel()->data(index).toDouble();
		index = resultModel->GetModel()->index(i, c2);
		double f_2 = resultModel->GetModel()->data(index).toDouble();

		ids.push_back( id );
		f1.push_back( f_1 );
		f2.push_back( f_2 );

		//Also gather min/max
		f1_max = (f_1) > (f1_max) ? f_1 : f1_max;
		f2_max = (f_2) > (f2_max) ? f_2 : f2_max;
		f1_min = (f_1) < (f1_min) ? f_1 : f1_min;
		f2_min = (f_2) < (f2_min) ? f_2 : f2_min;
	}

	//Normalize/Scale the Features Data:
	double upper = 1;
	double lower = -1;
	double desired_range = upper-lower;
	double f1_act_range = f1_max - f1_min;
	double f2_act_range = f2_max - f2_min;

	for(int i=0; i<(int)ids.size(); ++i)
	{
		f1.at(i) = lower + desired_range * (f1.at(i) - f1_min) / f1_act_range;
		f2.at(i) = lower + desired_range * (f2.at(i) - f2_min) / f2_act_range;
	}


	//Set the Parameters
	struct svm_parameter param;
	param.svm_type = ONE_CLASS;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 1.0/double(ids.size());	// 1/k
	//param.gamma = 1;
	param.coef0 = 0;
	param.nu = 0.1;
	param.cache_size = 100;
	param.C = 1;
	param.eps = .001;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

	//Create the problem:
	struct svm_problem prob;

	prob.l = resultModel->NumObjects();				//Number of objects
	prob.y = Malloc(double,prob.l);					//Array Containing target values (unknowns)
	prob.x = Malloc(struct svm_node *,prob.l);		//Array of Pointers to Nodes

	for(int i=0; i<prob.l; i++)
	{
		prob.y[i] = ids.at(i);				//This is the label (target) and it is unknown

		struct svm_node *x_space = Malloc(struct svm_node,3);			//Individual node
		x_space[0].index = 1;
		x_space[0].value = f1.at(i);
		x_space[1].index = 2;
		x_space[1].value = f2.at(i);
		x_space[2].index = -1;

		prob.x[i] = &x_space[0];	//Point to this new set of nodes.
	}

	//Now train
	struct svm_model *model;
	model = svm_train(&prob,&param);

	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);

	//Predict:
	vector<int> outliers;

	for(int i=0; i<prob.l; i++)
	{
		struct svm_node *x = Malloc(struct svm_node,3);			//Individual node
		x[0].index = 1;
		x[0].value = f1.at(i);
		x[1].index = 2;
		x[1].value = f2.at(i);
		x[2].index = -1;

		double v = svm_predict(model,x);
		free(x);

		if( v == -1 )
		{
			outliers.push_back( ids.at(i) );
		}
	}
	svm_destroy_model(model);

	//Now set the outliers:
	resultModel->SetOutliers(outliers);

}
