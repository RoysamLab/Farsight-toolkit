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

PatternAnalysisWizard::PatternAnalysisWizard(
  vtkSmartPointer<vtkTable> table, Module mod, const char * trainColumn,
  const char * resultColumn, QWidget *parent)
	: QWizard(parent)
{
	this->m_table = table;

	this->columnForTraining = trainColumn;
	this->columnForPrediction = resultColumn;
	this->m_module = mod;

	optionGroup = new QButtonGroup;
	optionGroup->setObjectName("optionGroup");
	initOptionGroup();

	featureGroup = new QButtonGroup;
	featureGroup->setObjectName("featureGroup");
	initFeatureGroup();
	StartPage *startPage = new StartPage(optionGroup);
	startPage->setObjectName("startPage");
	this->setPage(Page_Start, startPage);
	FeaturesPage *featuresPage = new FeaturesPage(featureGroup);
	featuresPage->setObjectName("featuresPage");
	this->setPage(Page_Features, featuresPage);
	//this->setPage(Page_Execute, new ExecutePage();

	this->button(QWizard::FinishButton)->setObjectName("doneButton");
	this->button(QWizard::FinishButton)->parent()->setObjectName("wizardWidget");

	this->setStartId(Page_Features);
	this->setModal(false);

	//this->setOption(QWizard::HaveCustomButton1);
	//this->setButtonText(QWizard::CustomButton1,"Execute");
	//connect(this, SIGNAL(customButtonClicked(int)), this, SLOT(executeNextStep(int)));

	this->setOption(QWizard::NoBackButtonOnStartPage,true);
	//setOption(HaveHelpButton, true);
	//setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	//connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));
	extractedTable = false;
	
 }
//******************************************************************************************
PatternAnalysisWizard::PatternAnalysisWizard(
  vtkSmartPointer<vtkTable> table, vtkSmartPointer<vtkTable> model_table , QString fileName, Module mod, const char * trainColumn,
  const char * resultColumn, QWidget *parent)
	: QWizard(parent)
{
	this->m_table = table;
	this->mod_table = model_table;
	this->filename = fileName;

	this->columnForTraining = trainColumn;
	this->columnForPrediction = resultColumn;
	this->m_module = mod;

	optionGroup = new QButtonGroup;
	optionGroup->setObjectName("optionGroup");
	initOptionGroup();

	featureGroup = new QButtonGroup;
	featureGroup->setObjectName("featureGroup");
	disabledFeatureGroup();

	StartPage *startPage = new StartPage(optionGroup);
	startPage->setObjectName("startPage");
	this->setPage(Page_Start, startPage);
	FeaturesPage *featuresPage = new FeaturesPage(featureGroup,1);
	featuresPage->setObjectName("featuresPage");
	this->setPage(Page_Features, featuresPage);
	//this->setPage(Page_Execute, new ExecutePage();

	this->setStartId(Page_Features);
	this->setModal(false);

	//this->setOption(QWizard::HaveCustomButton1);
	//this->setButtonText(QWizard::CustomButton1,"Execute");
	//connect(this, SIGNAL(customButtonClicked(int)), this, SLOT(executeNextStep(int)));

	this->setOption(QWizard::NoBackButtonOnStartPage,true);
	//setOption(HaveHelpButton, true);
	//setPixmap(QWizard::LogoPixmap, QPixmap(":/images/logo.png"));
	//connect(this, SIGNAL(helpRequested()), this, SLOT(showHelp()));

	
 }
//*****************************************************************************************

void PatternAnalysisWizard::initFeatureGroup(void)
{
	if(!m_table) return;

	featureGroup->setExclusive(false);
	bool zernike_added = false;
	for (int c=0; c<m_table->GetNumberOfColumns(); ++c)
	{
		const char * name = m_table->GetColumnName(c);
		std::string current_name = name;
		if( current_name.find("train")!=std::string::npos || current_name.find("prediction")!=std::string::npos || current_name.find("confidence")!=std::string::npos || current_name.find("time")!=std::string::npos)
			continue;
		if( current_name.find("Zern")!=std::string::npos)
		{
			if(zernike_added == false)
			{
				QCheckBox * check = new QCheckBox("zernike moments");
				check->setObjectName("zernikeMoments");
				check->setChecked(true);
				featureGroup->addButton(check, 0);
				zernike_added = true;
			}
			continue;
		}
		QCheckBox * check = new QCheckBox(QString(name));
		check->setObjectName(QString(name));
		check->setChecked(true);
		featureGroup->addButton(check, c);
	}
}

void PatternAnalysisWizard::disabledFeatureGroup(void)
{
	if(!m_table) return;

	featureGroup->setExclusive(false);
	bool zernike_added = false;
	for (int d=0; d<mod_table->GetNumberOfColumns(); ++d)
	{
	    const char * mod_name = mod_table->GetColumnName(d);
	    std::string model_name = mod_name;
	    if(model_name.find("Zern")!=std::string::npos)
		{
			if(zernike_added == false)
			{
				QCheckBox * check = new QCheckBox("zernike moments");
				check->setObjectName("zernikeMoments");
				check->setCheckable(false);  
				featureGroup->addButton(check);
				zernike_added = true;
			}
			continue;
		}
		QCheckBox * check = new QCheckBox(QString(mod_name));
		check->setObjectName(mod_name);
		check->setCheckable(false);  
		featureGroup->addButton(check);		
	}	
}

void PatternAnalysisWizard::initOptionGroup(void)
{
	optionGroup->setExclusive(true);
	QRadioButton *outlierButton = new QRadioButton(tr("Detect Outliers (using SVM)"));
	outlierButton->setObjectName("outlierButton");
	QRadioButton *classifyButton = new QRadioButton(tr("Classify (using KPLS)"));
	classifyButton->setObjectName("classifyButton");
	QRadioButton *createTrainButton = new QRadioButton(tr("Create Training Model... (using SEGMODEL)"));
	createTrainButton->setObjectName("createTrainButton");
	QRadioButton *appendTrainButton = new QRadioButton(tr("Append Training Model... (using APPENDMODEL)"));
	appendTrainButton->setObjectName("appendTrainButton");
	QRadioButton *activeButton = new QRadioButton(tr("Choose Features for Active Learning..."));
	activeButton->setObjectName("activeButton");
	QRadioButton *activeModelButton = new QRadioButton(tr("Extract Table From Active Model..."));
	activeModelButton->setObjectName("activeModelButton");
	QRadioButton *clusButton = new QRadioButton(tr("Select Features for Clustering..."));
	clusButton->setObjectName("clusButton");

	optionGroup->addButton(outlierButton, 0);
	optionGroup->addButton(classifyButton, 1);
	optionGroup->addButton(createTrainButton, 2);
	optionGroup->addButton(appendTrainButton, 3);
	optionGroup->addButton(activeButton, 4);
	optionGroup->addButton(activeModelButton, 5);
	optionGroup->addButton(clusButton, 6);

	switch(m_module)
		{
		    case _SVM: 
				outlierButton->setChecked(true);
				break;
			case _KPLS:
			    classifyButton->setChecked(true);
		        break;
			case _SEGMODEL:
				createTrainButton->setChecked(true);
				break;
			case _APPENDMODEL:
				appendTrainButton->setChecked(true);
				break;
			case _ACTIVE:
				activeButton->setChecked(true);
				break;
			case _ACTIVEMODEL:
				activeModelButton->setChecked(true);
				break;
			case _CLUS:
				clusButton->setChecked(true);
				break;
		}
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
			switch(id)
			{
			    case 0: 
					runSVM();
					break;
				case 1:
				    runKPLS();
			        break;
				case 2:
					saveModel();
					break;
				case 3:
					appendModel(mod_table, filename);
					break;
				case 4:
					extractTable(true);
					break;
				case 5:
					extractTableFromModel(mod_table);
					break;
				case 6:
					extractTable(true);
					break;

			}
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
FeaturesPage::FeaturesPage(QButtonGroup *fGroup, int caller, QWidget *parent)
	: QWizardPage(parent)
{
	if(caller == 0)
		setTitle(tr("Choose Features"));
	else
        setTitle(tr("Selected Features"));

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addWidget( initFeatureBox(fGroup) );
	setLayout(layout);
}
//selectedFeaturesPage::selectedFeaturesPage(QButtonGroup *fGroup, QWidget *parent)
//	: QWizardPage(parent)
//{
//	setTitle(tr("Selected Features"));
//
//	QVBoxLayout *layout = new QVBoxLayout;
//	layout->addWidget( initFeatureBox(fGroup) );
//	setLayout(layout);
//}
QGroupBox * FeaturesPage::initFeatureBox(QButtonGroup *fGroup)
{
	featureGroup = fGroup;

	QGroupBox *groupBox = new QGroupBox(tr("Features"));
	groupBox->setObjectName("Features");
	QVBoxLayout *featureLayout = new QVBoxLayout;

	QHBoxLayout *selLayout = new QHBoxLayout;
	QLabel *selLabel = new QLabel(tr("Select:"));
	QPushButton *allButton = new QPushButton(tr("All"));
	allButton->setObjectName("allButton");
	QPushButton *nonButton = new QPushButton(tr("None"));
	nonButton->setObjectName("nonButton");
	connect(allButton, SIGNAL(clicked()), this, SLOT(selectAll()));
	connect(nonButton, SIGNAL(clicked()), this, SLOT(selectNone()));
	selLayout->addWidget(selLabel);
	selLayout->addWidget(allButton);
	selLayout->addWidget(nonButton);
	selLayout->addStretch(10);
	featureLayout->addLayout(selLayout);

	QWidget *groupWidget = new QWidget;
	groupWidget->setObjectName("groupWidget");
	QVBoxLayout *vLayout = new QVBoxLayout;

	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		vLayout->addWidget( buttons.at(b) );
	}
	groupWidget->setLayout(vLayout);

	QScrollArea *scrollArea = new QScrollArea;
	scrollArea->setObjectName("scrollArea");
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

//****************************************************************************
// ExecutePage
//****************************************************************************
/*
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

	std::vector< std::vector< double > > features;	//Will contain the normalized features
	features.resize( m_table->GetNumberOfRows() );

	//exract data from the model and get min/max values:
	for(int r=0; r<(int)features.size(); ++r)
	{
		for(int c=0; c<(int)columnsToUse.size(); ++c)
		{
			int col = columnsToUse.at(c);
			double val = m_table->GetValue(r,col).ToDouble();
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

	prob.l = (int)features.size();					//Number of objects
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

	//If need to create a new column do so now:
	vtkAbstractArray * output = m_table->GetColumnByName(columnForPrediction);
	if(output == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForPrediction );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		m_table->AddColumn(column);
	}
	for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  //Set all values to 0
	{
		m_table->SetValueByName(row,columnForPrediction, 0);
	}
	for(int i = 0; i < (int)outliers.size(); ++i)					//Set outliers to 1
	{
		int row = outliers.at(i);
		m_table->SetValueByName(row,columnForPrediction, 1);
	}
	
	emit changedTable();

	/*
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
	*/
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
// RunKPLS - run kpls to classify objects
//****************************************************************************
void PatternAnalysisWizard::runKPLS()
{
#ifdef USE_KPLS
	//Find out which features are checked (which columns to use).
	std::vector<int> columnsToUse;
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		if( buttons.at(b)->isChecked() )
		{
			if(featureGroup->id(buttons.at(b)) != 0)
			{
				columnsToUse.push_back( featureGroup->id(buttons.at(b)) );
			}
			else
			{
				for(int col=0; col<(int)m_table->GetNumberOfColumns(); ++col)
				{
					std::string col_name = m_table->GetColumnName(col);
					if(col_name.find("Zern")!=std::string::npos)
					{
						columnsToUse.push_back(col);
					}
				}
			}
		}
	}
	
	if(columnsToUse.size() <= 0)
		return;

	this->KPLSrun(columnsToUse);

	emit changedTable();
	emit enableModels();
#else
	QMessageBox::information(this, tr("MESSAGE"), tr("FARSIGHT was not compiled with KPLS library"));
#endif
}

void PatternAnalysisWizard::KPLSrun(std::vector<int> columnsToUse){
#ifdef USE_KPLS
	//Setup up the kpls:
	KPLS *kpls = new KPLS();
	kpls->SetLatentVars(5);
	kpls->SetSigma(20);

	int num_rows = m_table->GetNumberOfRows();
	int num_cols = (int)columnsToUse.size();

	MATRIX data = kpls->GetDataPtr(num_rows, num_cols);
	VECTOR ids = kpls->GetIDPtr();
	VECTOR training = kpls->GetTrainingPtr();

	std::set<int> outcomes;
	//extract data from the table:
	for(int r=0; r<num_rows; ++r)
	{
		ids[r] = m_table->GetValue(r,0).ToDouble();
		for(int c=0; c<num_cols; ++c)
		{
			double val = m_table->GetValue(r, columnsToUse.at(c)).ToDouble();
			data[r][c] = val;
		}
		training[r] = m_table->GetValueByName(r,columnForTraining).ToDouble();
		outcomes.insert(int(training[r]));
	}

	if(outcomes.size() <= 1)
	{
		delete kpls;
		return;
	}

	kpls->InitVariables();
	kpls->ScaleData();
	kpls->Train();
	kpls->Classify();

	VECTOR predictions = kpls->GetPredictions();

	//If need to create a new column do so now:
	vtkAbstractArray * output = m_table->GetColumnByName(columnForPrediction);
	if(output == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForPrediction );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		m_table->AddColumn(column);
	}

	for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  
	{
		m_table->SetValueByName(row, columnForPrediction, vtkVariant(predictions[row]));
	}

	delete kpls;
#else
	std::cerr<<"Compiled without the KPLS library"<<std::endl;
#endif
}

PatternAnalysisWizardNoGUI::PatternAnalysisWizardNoGUI(vtkSmartPointer<vtkTable> table, const char * trainColumn, const char * resultColumn){
	this->m_table = table;
	this->columnForTraining = trainColumn;
	this->columnForPrediction = resultColumn;
}

void PatternAnalysisWizardNoGUI::KPLSrun1(std::vector<int> columnsToUse){
#ifdef USE_KPLS
	//Setup up the kpls:
	KPLS *kpls = new KPLS();
	kpls->SetLatentVars(5);
	kpls->SetSigma(20);

	int num_rows = m_table->GetNumberOfRows();
	int num_cols = (int)columnsToUse.size();

	MATRIX data = kpls->GetDataPtr(num_rows, num_cols);
	VECTOR ids = kpls->GetIDPtr();
	VECTOR training = kpls->GetTrainingPtr();

	std::set<int> outcomes;
	//extract data from the table:
	for(int r=0; r<num_rows; ++r)
	{
		ids[r] = m_table->GetValue(r,0).ToDouble();
		for(int c=0; c<num_cols; ++c)
		{
			double val = m_table->GetValue(r, columnsToUse.at(c)).ToDouble();
			data[r][c] = val;
		}
		training[r] = m_table->GetValueByName(r,columnForTraining).ToDouble();
		outcomes.insert(int(training[r]));
	}

	if(outcomes.size() <= 1)
	{
		delete kpls;
		return;
	}

	kpls->InitVariables();
	kpls->ScaleData();
	kpls->Train();
	kpls->Classify();

	VECTOR predictions = kpls->GetPredictions();

	//If need to create a new column do so now:
	vtkAbstractArray * output = m_table->GetColumnByName(columnForPrediction);
	if(output == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForPrediction );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		m_table->AddColumn(column);
	}

	for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  
	{
		m_table->SetValueByName(row, columnForPrediction, vtkVariant(predictions[row]));
	}

	delete kpls;
#else
	std::cerr<<"Compiled without the KPLS library"<<std::endl;
#endif
}

//****************************************************************************
// Save Segmentation Model
//****************************************************************************
void PatternAnalysisWizard::saveModel(void)
{
    //Selected Features
	std::vector<int> columnsToUse;
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
    //for(int b = 0; b<10; ++b)
	{
		//buttons.at(b)->setCheckable(false);
		if( buttons.at(b)->isChecked() )
		{
			if(featureGroup->id(buttons.at(b)) != 0)
			{
				columnsToUse.push_back( featureGroup->id(buttons.at(b)) );
			}
			else
			{
				for(int col=0; col<(int)m_table->GetNumberOfColumns(); ++col)
				{
					std::string col_name = m_table->GetColumnName(col);
					if(col_name.find("Zern")!=std::string::npos)
					{
						columnsToUse.push_back(col);
					}
				}
			}
		}
	}	
	if(columnsToUse.size() <= 0)
		return;


	//Save Dialog
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Training Model As..."),lastPath, tr("TEXT(*.txt)"));
	if(filename == "")
		return;
	lastPath = QFileInfo(filename).absolutePath();
	
	//Load the model into a new table 
    vtkSmartPointer<vtkTable> new_table = vtkSmartPointer<vtkTable>::New();
	new_table->Initialize();
   
	for(int c=0; c<(int)columnsToUse.size(); ++c)
	{
	    vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( m_table->GetColumnName(columnsToUse.at(c)) );
		new_table->AddColumn(column);
	}
	vtkSmartPointer<vtkDoubleArray> column1 = vtkSmartPointer<vtkDoubleArray>::New();
	column1->SetName( "Class" );
	new_table->AddColumn(column1);
	vtkSmartPointer<vtkDoubleArray> column2 = vtkSmartPointer<vtkDoubleArray>::New();
	column2->SetName( "ID" );
	new_table->AddColumn(column2);
   	for(int row = 0; row < (int)m_table->GetNumberOfRows(); ++row)
	{		
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c =0;c<(int)columnsToUse.size();++c)
			model_data1->InsertNextValue(m_table->GetValue(row,columnsToUse.at(c)));
		for(int col=((int)m_table->GetNumberOfColumns())-1; col>=0; --col)
		{	
			std::string current_column = m_table->GetColumnName(col);
			if(current_column.find("prediction") != std::string::npos )
			{
				model_data1->InsertNextValue(m_table->GetValue(row, col));
				break;
			}	
		}
		model_data1->InsertNextValue(m_table->GetValueByName(row, "ID"));
		new_table->InsertNextRow(model_data1);
	}
	//*********************************************************************
	//Save Model
	std::string Filename = filename.toStdString();
	
	//This function writes the features to a text file
	ofstream outFile; 
	outFile.open(Filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return;
	}
	//Write the headers:
	for(int c=0; c<new_table->GetNumberOfColumns(); ++c)
	{
		outFile << new_table->GetColumnName(c) << "\t";
	}
	outFile << "\n";
	//Write out the features:
	std::string current_column;
	for(int row = 0; row < new_table->GetNumberOfRows(); ++row)
	{
		for(int c=0; c < new_table->GetNumberOfColumns(); ++c)
		{
			std::stringstream out;
			current_column = new_table->GetColumnName(c);
			if((current_column.compare("ID") == 0) || (current_column.compare("Class") == 0))
				out << std::fixed << new_table->GetValue(row,c).ToInt();
			else
	            out << std::setprecision(3) << std::fixed << new_table->GetValue(row,c).ToFloat();
	        outFile << out.str() << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
    //******************************************************************************************
}




//****************************************************************************
// Extract the table with the features selected by the user
//****************************************************************************
void PatternAnalysisWizard::extractTable(bool flag)
{	

	extractedTable = true;
    //Selected Features
	std::vector<int> columnsToUse;
	QList<QAbstractButton *> buttons = featureGroup->buttons();
	for(int b = 0; b<buttons.size(); ++b)
	{
		//buttons.at(b)->setCheckable(false);
		
		if( buttons.at(b)->isChecked() )
		{
				columnsToUse.push_back( featureGroup->id(buttons.at(b)) );
				for(int col=0; col<(int)m_table->GetNumberOfColumns(); ++col)
				{

					std::string col_name = m_table->GetColumnName(col);
					if(col_name.find("Zern")!=std::string::npos)
					{
						columnsToUse.push_back(col);
					}
				}
			}
		}
		
	if(columnsToUse.size() <= 0)
		return;


	//Get the new table
    new_table = vtkSmartPointer<vtkTable>::New();
	new_table->Initialize();
   
	for(int c=0; c<(int)columnsToUse.size(); ++c)
	{
		
	    vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( m_table->GetColumnName(columnsToUse.at(c)) );
		new_table->AddColumn(column);
	}

   	for(int row = 0; row < (int)m_table->GetNumberOfRows(); ++row)
	{		
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c =0;c<(int)columnsToUse.size();++c)
			model_data1->InsertNextValue(m_table->GetValue(row,columnsToUse.at(c)));
		new_table->InsertNextRow(model_data1);
	}

	if(flag)
	{
		emit start_training(new_table);
	}

}


//****************************************************************************
// Extract the table with the features present in the model
//****************************************************************************
void PatternAnalysisWizard::extractTableFromModel(vtkSmartPointer<vtkTable> mod_table)
{	

	extractedTable = true;
    //Selected Features
	
	new_table = vtkSmartPointer<vtkTable>::New();
	new_table->Initialize();

	for(int col=0; col<(int)mod_table->GetNumberOfColumns(); ++col)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( mod_table->GetColumnName(col) );
		new_table->AddColumn(column);
	}

	for(int row=0; row<(int)m_table->GetNumberOfRows(); ++row)
	{
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int col=0; col<(int)mod_table->GetNumberOfColumns(); ++col)
		{
			std::string column_name = mod_table->GetColumnName(col);
			model_data1->InsertNextValue(m_table->GetValueByName(row,column_name.c_str()));

		}
		new_table->InsertNextRow(model_data1);
	}
}

//****************************************************************************
// Append Segmentation Model
//****************************************************************************
void PatternAnalysisWizard::appendModel(vtkSmartPointer<vtkTable> mod_table, QString filename)
{
	int r = mod_table->GetValue((int)mod_table->GetNumberOfRows()-1, (int)mod_table->GetNumberOfColumns()-1).ToInt();
   	for(int row = 0; row < (int)m_table->GetNumberOfRows(); ++row)
	{		
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c = 0;c<(int)mod_table->GetNumberOfColumns();++c)
		{
			 std::string mod_column = mod_table->GetColumnName(c);
			 if(mod_column.compare("Class") == 0)
			 {
				 for(int col=((int)m_table->GetNumberOfColumns())-1; col>=0; --col)
				{	
					std::string current_column = m_table->GetColumnName(col);
					if(current_column.find("prediction") != std::string::npos )
					{
						model_data1->InsertNextValue(m_table->GetValue(row, col));
						break;
					}	
				}
			 }

			 else
			 {
		        for(int d=0; d<(int)m_table->GetNumberOfColumns(); ++d)
		        {
			         std::string m_column = m_table->GetColumnName(d);
			         if( m_column.compare(mod_column) == 0)
				     {
					     if(mod_column.compare("ID") == 0)
						 {
							 int id = r + m_table->GetValue(row,d).ToInt();
							 model_data1->InsertNextValue(vtkVariant(id));
						 }
						 else
							 model_data1->InsertNextValue(m_table->GetValue(row,d));
				         break;
				     }
		         }
			 }
		}
		mod_table->InsertNextRow(model_data1);
	}
	//*********************************************************************
	//Save Model
	std::string Filename = filename.toStdString();
	
	//This function writes the features to a text file
	ofstream outFile; 
	outFile.open(Filename.c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return;
	}
	//Write the headers:
	for(int c=0; c<mod_table->GetNumberOfColumns(); ++c)
	{
		outFile << mod_table->GetColumnName(c) << "\t";
	}
	outFile << "\n";
	//Write out the features:
	std::string current_column;
	for(int row = 0; row < mod_table->GetNumberOfRows(); ++row)
	{
		for(int c=0; c < mod_table->GetNumberOfColumns(); ++c)
		{
			std::stringstream out;
			current_column = mod_table->GetColumnName(c);
			if((current_column.compare("ID") == 0) || (current_column.compare("Class") == 0))
				out << std::fixed << mod_table->GetValue(row,c).ToInt();
			else
	            out << std::setprecision(3) << std::fixed << mod_table->GetValue(row,c).ToFloat();
	        outFile << out.str() << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
    //******************************************************************************************

}
