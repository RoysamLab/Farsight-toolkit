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
#include "LibSVMWidget.h"

//for DBL_MAX...
#include <float.h>

//Constructor
LibSVMWidget::LibSVMWidget(QAbstractItemModel *mod, QWidget *parent)
: QWidget(parent)
{
	model = mod;
	QHBoxLayout *allLayout = new QHBoxLayout;

	allLayout->addWidget( initFeatureBox() );


	QVBoxLayout *optionLayout = new QVBoxLayout;
	QGroupBox *optionBox = initOptionBox();
	optionLayout->addWidget(optionBox);

	optionLayout->addStretch(50);

	goButton = new QPushButton(tr("GO"));
	connect(goButton, SIGNAL(clicked()), this, SLOT(go()));
	optionLayout->addWidget(goButton);

	allLayout->addLayout(optionLayout);

	this->setLayout(allLayout);
	this->setWindowTitle(tr("SVM"));

	columnForSVM = model->columnCount();
}

QGroupBox * LibSVMWidget::initFeatureBox()
{
	QGroupBox *groupBox = new QGroupBox(tr("Features"));
	QVBoxLayout *featureLayout = new QVBoxLayout;

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
	featureLayout->addLayout(selLayout);

	QWidget *groupWidget = new QWidget;
	QVBoxLayout *vLayout = new QVBoxLayout;
	featureGroup = new QButtonGroup;
	featureGroup->setExclusive(false);

	for (int c=1; c<model->columnCount()-1; ++c)
	{
		QString name = model->headerData(c,Qt::Horizontal).toString();
		QCheckBox *check = new QCheckBox(name);
		check->setChecked(false);
		vLayout->addWidget(check);
		featureGroup->addButton(check, c);
	}

	groupWidget->setLayout(vLayout);

	QScrollArea *scrollArea = new QScrollArea;
	scrollArea->setWidget(groupWidget);
	scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	featureLayout->addWidget(scrollArea);

	groupBox->setLayout(featureLayout);
	return groupBox;
}

QGroupBox * LibSVMWidget::initOptionBox()
{
	QGroupBox *groupBox = new QGroupBox(tr("Options"));
	QVBoxLayout *vLayout = new QVBoxLayout;

	normalizeBox = new QCheckBox(tr("Normalize"));
	normalizeBox->setToolTip(tr("Normalize features to values between -1 and 1 before classification"));
	normalizeBox->setChecked(true);
	vLayout->addWidget(normalizeBox);

	QHBoxLayout *nuLayout = new QHBoxLayout;

	QLabel *nuLabel = new QLabel(tr("Nu:"));
	nuLabel->setToolTip(tr("Modify radius of basis function. Increasing Nu will decrease the radius"));
	nuLayout->addWidget(nuLabel);

	nuSpin = new QDoubleSpinBox;
	nuSpin->setDecimals(2);
	nuSpin->setSingleStep(0.01);
	nuSpin->setRange(0,1);
	nuSpin->setValue(0.10);
	nuLayout->addWidget(nuSpin);

	vLayout->addLayout(nuLayout);

	groupBox->setLayout(vLayout);
	return groupBox;
}

void LibSVMWidget::selectNone()
{
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			buttons.at(b)->setChecked(false);
		}
	}
}

void LibSVMWidget::selectAll()
{
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( !buttons.at(b)->isChecked() )
		{
			buttons.at(b)->setChecked(true);
		}
	}
}

void LibSVMWidget::go()
#ifdef USE_KPLS
{
	//Find out which features are checked (which columns to use).
	std::vector<int> columnsToUse;
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			columnsToUse.push_back( featureGroup->id(buttons.at(b)) );
		}
	}

	//Find out which column is for class

	//Setup up the kpls:
	KPLS *kpls = new KPLS();
	kpls->SetLatentVars(5);
	kpls->SetSigma(20);

	int num_rows = (int)model->rowCount();
	int num_cols = (int)columnsToUse.size();

	MATRIX data = kpls->GetDataPtr(num_rows, num_cols);
	VECTOR ids = kpls->GetIDPtr();
	VECTOR training = kpls->GetTrainingPtr();

	//extract data from the model
	QModelIndex index;
	for(int r=0; r<num_rows; ++r)
	{
		index = model->index(r,0);
		ids[r] = model->data(index).toDouble();
		for(int c=0; c<num_cols; ++c)
		{
			index = model->index(r, columnsToUse.at(c));
			double val = model->data(index).toDouble();
			data[r][c] = val;
		}
		index = model->index(r,model->columnCount()-1);
		training[r] = model->data(index).toDouble();
	}

	kpls->InitVariables();
	kpls->ScaleData();
	kpls->Train();
	kpls->Classify();

	VECTOR predictions = kpls->GetPredictions();

	//Add the outliers to the model!!
	if(columnForSVM >= model->columnCount())
	{
		model->insertColumn(columnForSVM);			//Add Column for svm result
		model->setHeaderData(columnForSVM, Qt::Horizontal, tr("predict"));
	}

	//stop signalling:
	model->blockSignals(true);
	for(int row = 1; (int)row < model->rowCount(); ++row)  //Set all values to 0
	{
		model->setData(model->index(row, columnForSVM), (int)predictions[row]);
	}
	//turn signals back on & change one more piece of data to force dataChanged signal
	model->blockSignals(false);
	model->setData(model->index(0, columnForSVM), (int)predictions[0]);

	delete kpls;
}
#else
{
	//Find out which featueres are checked (which columns to use).
	std::vector<int> columnsToUse;

	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			columnsToUse.push_back( featureGroup->id(buttons.at(b)) );
		}
	}

	//Setup the scaling values
	std::vector<double> f_max;						//The maximum value of each feature
	f_max.assign(columnsToUse.size(), -DBL_MAX);
	std::vector<double> f_min;						//The minimum value of each feature
	f_min.assign(columnsToUse.size(), DBL_MAX);

	std::vector< std::vector< double > > features;
	features.resize( model->rowCount() );

	//exract data from the model and get min/max values:
	QModelIndex index;
	for(int r=0; r<(int)model->rowCount(); ++r)
	{
		for(int c=0; c<(int)columnsToUse.size(); ++c)
		{
			index = model->index(r, columnsToUse.at(c));
			double val = model->data(index).toDouble();
			features.at(r).push_back( val );

			if( normalizeBox->isChecked() )
			{
				//Also gather min/max
				f_max.at(c) = val > f_max.at(c) ? val : f_max.at(c);
				f_min.at(c) = val < f_min.at(c) ? val : f_min.at(c);
			}
		}
	}

	if( normalizeBox->isChecked() )
	{
		//Normalize/Scale the Features Data:
		double upper = 1;
		double lower = -1;
		double desired_range = upper-lower;

		for(int r=0; r<(int)features.size(); ++r)
		{
			for(int c=0; c<(int)columnsToUse.size(); ++c)
			{
				double oldval = features.at(r).at(c);
				features.at(r).at(c) = lower + desired_range * (oldval - f_min.at(c)) / (f_max.at(c) - f_min.at(c));
			}
		}
	}

	//Create the libSVM problem:
	#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
	struct svm_problem prob;

	prob.l = model->rowCount();							//Number of objects
	prob.y = Malloc(double,prob.l);					//Array Containing target values (unknowns)
	prob.x = Malloc(struct svm_node *,prob.l);		//Array of Pointers to Nodes

	for(int r=0; r<prob.l; ++r)
	{
		prob.y[r] = 1;								//This is the label (target) and it is unknown

		struct svm_node *x_space = Malloc(struct svm_node, columnsToUse.size()+1);			//Individual node

		for(int c=0; c<(int)columnsToUse.size(); ++c)
		{
			x_space[c].index = c+1;
			x_space[c].value = features.at(r).at(c);
		}
		x_space[ columnsToUse.size() ].index = -1;
		prob.x[r] = &x_space[0];	//Point to this new set of nodes.
	}

	//Set the Parameters
	struct svm_parameter param;
	param.svm_type = ONE_CLASS;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 1.0/double(prob.l);	// 1/k
	//param.gamma = 1;
	param.coef0 = 0;
	//param.nu = 0.1;
	param.nu = nuSpin->value();
	param.cache_size = 100;
	param.C = 1;
	param.eps = .001;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	//Now train
	struct svm_model *m_svm_model;
	m_svm_model = svm_train(&prob,&param);

	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);

	//Predict:
	std::vector<int> outliers;

	for(int r=0; r<prob.l; r++)
	{
		struct svm_node *x = Malloc(struct svm_node,columnsToUse.size()+1);			//Individual node
		for(int c=0; c<(int)columnsToUse.size(); ++c)
		{
			x[c].index = c+1;
			x[c].value = features.at(r).at(c);
		}
		x[ columnsToUse.size() ].index = -1;

		double v = svm_predict(m_svm_model,x);
		free(x);

		if( v == -1 )
		{
			outliers.push_back( r );
		}
	}
	svm_destroy_model(m_svm_model);

	//Add the outliers to the model!!
	if(columnForSVM >= model->columnCount())
	{
		model->insertColumn(columnForSVM);			//Add Column for svm result
		model->setHeaderData( columnForSVM, Qt::Horizontal, tr("outlier?") );
	}

	//stop signalling:
	model->blockSignals(true);

	int z = 0;
	int o = 1;
	for(int row = 0; (int)row < model->rowCount(); ++row)  //Set all values to 0
	{
		model->setData(model->index(row, columnForSVM), z);
	}

	for(int i = 0; i < (int)outliers.size()-1; ++i)							//Set outliers to 1
	{
		model->setData(model->index(outliers.at(i), columnForSVM), o);
	}

	//turn signals back on & change one more piece of data to force dataChanged signal
	model->blockSignals(false);
	model->setData(model->index(outliers.back(), columnForSVM), o);
}
#endif
