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

	//this->scatter->SetColForColor(resultModel->ColumnForColor(), resultModel->ColorMap());
	this->updateOptionMenus();
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

	QWidget *centralWidget = new QWidget();
	QVBoxLayout *vlayout = new QVBoxLayout();
	scatter = new ScatterView();
	vlayout->addWidget(scatter);
	centralWidget->setLayout(vlayout);
	this->setCentralWidget(centralWidget);

	setWindowTitle(tr("Scatter Plot"));
	setAttribute ( Qt::WA_DeleteOnClose );
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

void PlotWindow::updateOptionMenus()
{
	//QStandardItemModel *model = resultModel->GetModel();
	QAbstractItemModel *model = scatter->model();

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

		if(c==0)
			xChange(xAct);
		if(c==1)
			yChange(yAct);
		if(name == "class")
			colorChange(cAct);
	}
}


/*
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
*/

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