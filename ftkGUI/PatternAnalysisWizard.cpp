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

#include "PatternAnalysisWizard.h"

PatternAnalysisWizard::PatternAnalysisWizard(QAbstractItemModel *mod, int outputColumn, QWidget *parent)
	: QWizard(parent)
{
	model = mod;
	this->columnForPrediction = outputColumn;
	if(columnForPrediction >= model->columnCount())
	{
		model->insertColumn(columnForPrediction);			//Add Column for svm result
		model->setHeaderData(columnForPrediction, Qt::Horizontal, tr("pattern") );
	}
	
	optionGroup = new QButtonGroup;
	initOptionGroup();

	featureGroup = new QButtonGroup;
	initFeatureGroup();

	this->setPage(Page_Start, new StartPage(optionGroup));
	this->setPage(Page_Features, new FeaturesPage(featureGroup));
	//this->setPage(Page_Execute, new ExecutePage();
	//this->setPage(Page_Binarize, new BinarizePage);
	//this->setPage(Page_Seeds, new SeedsPage);
	//this->setPage(Page_Cluster, new ClusterPage);
	//this->setPage(Page_Finalize, new FinalizePage);

	this->setStartId(Page_Start);
	this->setModal(true);

	//this->setOption(QWizard::HaveCustomButton1);
	//this->setButtonText(QWizard::CustomButton1,"Execute");
	//connect(this, SIGNAL(customButtonClicked(int)), this, SLOT(executeNextStep(int)));

	this->setOption(QWizard::NoBackButtonOnStartPage,true);
	//setOption(HaveHelpButton, true);
	//setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	//connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));

	this->setWindowTitle(tr("Pattern Analysis Wizard"));
 }

void PatternAnalysisWizard::initFeatureGroup(void)
{
	if(!model) return;

	featureGroup->setExclusive(false);
	for (int c=1; c<model->columnCount()-1; ++c)
	{
		QString name = model->headerData(c,Qt::Horizontal).toString();
		QCheckBox *check = new QCheckBox(name);
		check->setChecked(false);
		featureGroup->addButton(check, c);
	}
}

void PatternAnalysisWizard::initOptionGroup(void)
{
	optionGroup->setExclusive(true);
	QRadioButton *outlierButton = new QRadioButton(tr("Detect Outliers (using SVM)"));
	QRadioButton *classifyButton = new QRadioButton(tr("Classify (using KPLS)"));
	optionGroup->addButton(outlierButton, 0);
	optionGroup->addButton(classifyButton, 1);
	outlierButton->setChecked(true);
}

/*
//is called by QWizard to prepare page id just before it is shown as a result of the user clicking Next
void PatternAnalysisWizard::initializePage(int id)
{
	switch(id)
	{
		case Page_Start:
			break;
	}
}

//is called by QWizard to clean up page id just before the user leaves it by clicking Back
void PatternAnalysisWizard::cleanupPage(int id)
{
	switch(id)
	{
		case Page_Start:
			break;
	}
}
*/

int PatternAnalysisWizard::nextId() const
{
	switch( this->currentId() )
	{
		case Page_Start:
			return Page_Features;
			break;
		//case Page_Features:
			//return Page_Execute;
			//break;
	}
	return -1;
}

bool PatternAnalysisWizard::validateCurrentPage()
{
	switch( this->currentId() )
	{
		//********************************************************************
		// Run the desired module upon exiting the wizard:
		case Page_Features:
			int id = optionGroup->checkedId();
			if(id == 0)
				runSVM();
			else if(id == 1)
				runKPLS();
			return true;
		break;
		//********************************************************************
	}
	return true;
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
// PAGES:
//****************************************************************************

//****************************************************************************
// StartPage
//****************************************************************************
StartPage::StartPage(QButtonGroup *oGroup, QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Select Module"));
	//setPixmap(QWizard::WatermarkPixmap, QPixmap(":/images/watermark.png"));
	
	QVBoxLayout *layout = new QVBoxLayout;

	QList<QAbstractButton *> buttons = oGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		layout->addWidget( buttons.at(b) );
	}
	setLayout(layout);
}
//****************************************************************************

//****************************************************************************
// FeaturesPage
//****************************************************************************
FeaturesPage::FeaturesPage(QButtonGroup *fGroup, QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Choose Features"));
	
	QVBoxLayout *layout = new QVBoxLayout;
	layout->addWidget( initFeatureBox(fGroup) );
	setLayout(layout);
}

QGroupBox * FeaturesPage::initFeatureBox(QButtonGroup *fGroup)
{
	featureGroup = fGroup;

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

	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		vLayout->addWidget( buttons.at(b) );
	}
	groupWidget->setLayout(vLayout);

	QScrollArea *scrollArea = new QScrollArea;
	scrollArea->setWidget(groupWidget);
	scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	featureLayout->addWidget(scrollArea);

	groupBox->setLayout(featureLayout);
	return groupBox;
}

void FeaturesPage::selectNone()
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

void FeaturesPage::selectAll()
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

bool FeaturesPage::isComplete() const
{
	//if(imageFileCombo->currentText() != "")
		return true;
	//else
	//	return false;
}
//****************************************************************************

/*
//****************************************************************************
// ExecutePage
//****************************************************************************
ExecutePage::ExecutePage(QWidget *parent)
	: QWizardPage(parent)
{
	setTitle(tr("Execute"));
	
	QVBoxLayout *layout = new QVBoxLayout;

	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			buttons.at(b)->setChecked(false);
		}
	}

	QLabel * text = new QLabel(tr("You have chosen to
	layout->addWidget(  );
	setLayout(layout);
}
*/

//****************************************************************************

//****************************************************************************
// .....Page
//****************************************************************************

//****************************************************************************

//****************************************************************************
// .....Page
//****************************************************************************

//****************************************************************************
//****************************************************************************
//****************************************************************************
// RunSVM - run support vector machine in one class mode to find outliers.
//****************************************************************************
void PatternAnalysisWizard::runSVM()
{
	bool normalize = true; //ALWAYS TRUE!!
	double nu = 0.10;

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

			if( normalize )
			{
				//Also gather min/max
				f_max.at(c) = val > f_max.at(c) ? val : f_max.at(c);
				f_min.at(c) = val < f_min.at(c) ? val : f_min.at(c);
			}
		}
	}

	if( normalize )
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
	param.nu = nu;
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

	model->setHeaderData( columnForPrediction, Qt::Horizontal, tr("outlier?") );

	//stop signalling:
	model->blockSignals(true);

	int z = 0;
	int o = 1;
	for(int row = 0; (int)row < model->rowCount(); ++row)  //Set all values to 0
	{
		model->setData(model->index(row, columnForPrediction), z);         
	}

	for(int i = 0; i < (int)outliers.size()-1; ++i)							//Set outliers to 1
	{
		model->setData(model->index(outliers.at(i), columnForPrediction), o);
	}

	//turn signals back on & change one more piece of data to force dataChanged signal
	model->blockSignals(false);
	model->setData(model->index(outliers.back(), columnForPrediction), o);
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
// RunKPLS - run kpls to classify objects
//****************************************************************************
void PatternAnalysisWizard::runKPLS()
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
		index = model->index(r,model->columnCount()-2);
		training[r] = model->data(index).toDouble();
	}

	kpls->InitVariables();
	kpls->ScaleData();
	kpls->Train();
	kpls->Classify();

	VECTOR predictions = kpls->GetPredictions();

	model->setHeaderData(columnForPrediction, Qt::Horizontal, tr("predict"));

	//stop signalling:
	model->blockSignals(true);
	for(int row = 1; (int)row < model->rowCount(); ++row)  //Set all values to 0
	{
		model->setData(model->index(row, columnForPrediction), (int)predictions[row]);         
	}
	//turn signals back on & change one more piece of data to force dataChanged signal
	model->blockSignals(false);
	model->setData(model->index(0, columnForPrediction), (int)predictions[0]);
	
	delete kpls;
}