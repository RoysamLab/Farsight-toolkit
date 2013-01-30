
#include "ALforNucEditor.h"

ALforNucEd::ALforNucEd(bool val)
{	
	pixel_class = val;
	confidence_thresh = 0.5;
	prediction_col_name = "prediction_active";
	confidence_col_name = "confidence";
	HeatmapWin = NULL;
	mclr = NULL;
}

ALforNucEd::~ALforNucEd()
{
	if(this->HeatmapWin)
	{
		delete this->HeatmapWin;
	}
	if(mclr)
	{
		delete mclr;
	}	
}

void ALforNucEd::SetLabelView(LabelImageViewQT *view)
{
	labelView = view;
	if(pixel_class)
	{
		std::map<int, ftk::Object::Point> * pixelLocationMap;
		for(int row=0; row<(int)trainingTable->GetNumberOfRows(); ++row)
		{
			int id = trainingTable->GetValue(row,0).ToInt();
			ftk::Object::Point centroid;
			centroid.x = trainingTable->GetValue(row,1).ToInt();
			centroid.y = trainingTable->GetValue(row,2).ToInt();
			centroid.z = 0;
			(*pixelLocationMap)[id] = centroid;
		}
		labelView->SetCenterMapPointer(pixelLocationMap);
	}

}

void ALforNucEd::RunALClassification(bool val)
{	
	from_model = val;
	if(!from_model)
	{
		id_time.clear();
		if(classificationTables.size() > 1)
		{
			labelView->SetCurrentTimeVal(0);
		}

		if(!trainingTable) return;


		//Default confidence threshold is 50 % 
		ClassNameConfidenceDialog *dialog = new ClassNameConfidenceDialog(this);
		if( dialog->exec() )
		{
			classification_name = dialog->getClassName();
			confidence_thresh = dialog->getConfThresh();		
		}
		else
			return;
		delete dialog;


		TrainingDialog *d = new TrainingDialog(trainingTable, "train","active",trainingTable->GetNumberOfRows() ,this);
		d->exec();


		//Clear the Gallery 
		//gallery.clear();	

		// Remove the training examples from the list of ids.
		//Get the list of ids
		for(int i=0;i<trainingTable->GetNumberOfRows(); ++i)
		{
			if(trainingTable->GetValueByName(i,"train_default1").ToDouble()==-1) 
			{
				std::pair<double,double> temp_pair;
				temp_pair.first = trainingTable->GetValue(i,0).ToDouble();
				if(classificationTables.size() == 1)
					temp_pair.second = 0;
				else
					temp_pair.second = trainingTable->GetValueByName(i,"time").ToDouble();
				id_time.push_back(temp_pair);
			}
		}

		if(classificationTables.size() > 1)
			trainingTable->RemoveColumnByName("time");

		// If the user did not hit cancel 
		if(d->result())
		{
			pWizard = new PatternAnalysisWizard( trainingTable, PatternAnalysisWizard::_ACTIVE,"","", this);
			pWizard->setObjectName("pWizard");
			connect(pWizard, SIGNAL(start_training(vtkSmartPointer<vtkTable>)), this, SLOT(Start_Training(vtkSmartPointer<vtkTable>)));
			pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
			pWizard->show();
		}
	}

	else
	{
		if(!trainingTable) return;

		//if(pWizard)
		//{
		//	delete pWizard;
		//}

		QString fileName  = QFileDialog::getOpenFileName(this, "Select training model to open", "",
			tr("TXT Files (*.txt)"));
		if(fileName == "")
			return;		

		vtkSmartPointer<vtkTable> active_model_table = ftk::LoadTable(fileName.toStdString());

		// to generate the Active Learning Matrix
		act_learn_matrix.set_size((int)active_model_table->GetNumberOfColumns() , (int)active_model_table->GetNumberOfRows() - 2);
		for(int row = 2; row<(int)active_model_table->GetNumberOfRows(); ++row)
		{
			for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
			{
				act_learn_matrix.put(col, row-2, active_model_table->GetValue(row,col).ToDouble());
			}
		}

		//to generate the std_deviation and the mean vectors
		std_dev_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
		mean_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
		for(int col=1; col<(int)active_model_table->GetNumberOfColumns(); ++col)
		{
			std_dev_vec.put(col-1, active_model_table->GetValue(0,col).ToDouble());
			mean_vec.put(col-1, active_model_table->GetValue(1,col).ToDouble());
		}

		active_model_table->RemoveRow(0);
		active_model_table->RemoveRow(0);
		active_model_table->RemoveColumn(0);

		ClassNameConfidenceDialog *dialog = new ClassNameConfidenceDialog(this);
		if( dialog->exec() )
		{
			classification_name = dialog->getClassName();
			confidence_thresh = dialog->getConfThresh();		
		}
		else
			return;
		delete dialog;

		pWizard = new PatternAnalysisWizard( trainingTable, active_model_table, "", PatternAnalysisWizard::_ACTIVEMODEL,"","", this);
		pWizard->setWindowTitle(tr("Pattern Analysis Wizard"));
		pWizard->exec();

		if(pWizard->result())
		{	
			// Extracted Table containing features of trained model 
			pawTable = pWizard->getExtractedTable();

			if(mclr)
			{
				delete mclr;
			}
			mclr = new MCLR();
			// Number of features and classes needed in "add_bias" fuction of MCLR
			mclr->Set_Number_Of_Classes((int)active_model_table->GetNumberOfRows());
			mclr->Set_Number_Of_Features((int)active_model_table->GetNumberOfColumns());
			Start_Classification();		
		}

	}
}


void ALforNucEd::Start_Training(vtkSmartPointer<vtkTable> pTable)
{
	//pawTable does not have the id column  nor the train_default column
	pawTable = pTable;
	
	//// Delete the prediction column if it exists
	std::vector<std::string> prediction_names = ftk::GetColumsWithString( "prediction_active" , pawTable);
	if(prediction_names.size()>0)
		pawTable->RemoveColumnByName("prediction_active");

	vnl_vector<double> class_list(pawTable->GetNumberOfRows()); 

	for(int row = 0; (int)row < pawTable->GetNumberOfRows(); ++row)  
	{
		class_list.put(row,vtkVariant(trainingTable->GetValueByName(row,"train_default1")).ToDouble());
	}

	if(mclr)
	{
		delete mclr;
	}
	mclr = new MCLR();
	double sparsity = 1;
	double max_info = -1e9;

	// Normalize the feature matrix
	vnl_matrix<double> Feats = mclr->Normalize_Feature_Matrix(mclr->tableToMatrix(pawTable, id_time));
	mclr->Initialize(Feats,sparsity,class_list,"",pawTable);
	mclr->Get_Training_Model();

	// Get the active query based on information gain
	int active_query = mclr->Active_Query();

	//active_queries = mclr->ALAMO(active_query);		
	active_queries = mclr->Submodular_AL(active_query,mclr->testData);

	if(classificationTables.size() > 1)
		labelView->SetCurrentTimeVal(mclr->id_time_val.at(active_query).second);

	std::vector<std::pair<int,int> > dummy_vector;
	dummy_vector.clear();
	ALDialogPopUP(true, dummy_vector);

}

void ALforNucEd::ALDialogPopUP(bool first_pop, std::vector<std::pair<int,int> > query_labels)
{
	bool user_stop_alDialog_flag = false;
	int atleast_one_chosen = 0;

	if(!first_pop)		
	{
		for(int i=0; i<active_queries.size(); ++i)
		{
			atleast_one_chosen = atleast_one_chosen + query_labels[i].second;
			if(query_labels[i].second == -1)
			{	
				QMessageBox::critical(this, tr("Oops"), tr("Please select a class for all the cells"));
				this->show();
				alDialog =  new ActiveLearningDialog(snapshots, mclr->test_table, mclr->numberOfClasses, active_queries, mclr->top_features);	
				connect(alDialog, SIGNAL(retrain(bool, std::vector<std::pair<int,int> >)), this, SLOT(ALDialogPopUP(bool, std::vector<std::pair<int,int> >)));
				connect(alDialog, SIGNAL(start_classification(bool)), this, SLOT(Start_Classification(bool)));
				alDialog->show();
				i=0;		
				return;
			}
		}

#ifndef _MSC_VER
		if(this->HeatmapWin == NULL)
		{
			this->HeatmapWin = new Heatmap();
		}	
		this->HeatmapWin->setModels(mclr->Rearrange_Table(pawTable));
		this->HeatmapWin->reRunClus();
		this->HeatmapWin->showGraphforNe();
#endif
		// Update the data & refresh the training model and refresh the Training alDialog 		
		mclr->Update_Train_Data(query_labels);

		// Update the gallery
		if(mclr->stop_training && atleast_one_chosen!=0)
		{
			QMessageBox msgBox;
			msgBox.setText("I understand the classification problem.");
			msgBox.setInformativeText("Do you want to stop training and classify ? ");
			msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
			msgBox.setDefaultButton(QMessageBox::Ok);
			int ret = msgBox.exec();

			switch (ret) 
			{
			case QMessageBox::Ok:
				// Save was clicked
				user_stop_alDialog_flag = true;
				break;
			case QMessageBox::Cancel:
				mclr->stop_training = false;
				break;
			default:
				// should never be reached
				break;
			}
		} 

		if(user_stop_alDialog_flag)
		{
			Start_Classification(true);
			return;
		}

		mclr->Get_Training_Model();
		int active_query = mclr->Active_Query();
		//active_queries = mclr->ALAMO(active_query);
		active_queries = mclr->Submodular_AL(active_query,mclr->testData);

	}// END if(!first_pop)		

	
	snapshots.resize(active_queries.size());
	// Collect all the snapshots
	for(int i=0;i<active_queries.size(); ++i)
		//for(int i=0;i<1; ++i)
	{	
		if(classificationTables.size() > 1)
			labelView->SetCurrentTimeVal(mclr->id_time_val.at(active_queries[i]).second);
		if(!pixel_class)
			snapshots[i] =(labelView->getSnapshotforID(mclr->id_time_val.at(active_queries[i]).first));	
		else
			snapshots[i] = labelView->getSnapshotforID_1(mclr->id_time_val.at(active_queries[i]).first);

	}
	
	
	//mclr->test_table is the pawTable obtained above
	alDialog =  new ActiveLearningDialog(snapshots, mclr->test_table, mclr->numberOfClasses, active_queries, mclr->top_features);
	connect(alDialog, SIGNAL(retrain(bool, std::vector<std::pair<int,int> >)), this, SLOT(ALDialogPopUP(bool, std::vector<std::pair<int,int> >)));
	connect(alDialog, SIGNAL(start_classification(bool)), this, SLOT(Start_Classification(bool)));
	alDialog->show();
	
	#ifdef	USE_Clusclus
	//this->HeatmapWin = new Heatmap();
	//this->HeatmapWin->setModels(pawTable, this->selection);
	//this->HeatmapWin->setPriority(mclr->Get_Feature_Order());
	//this->HeatmapWin->runClus();
	//this->HeatmapWin->showGraph();
	#endif

}

void ALforNucEd::Start_Classification(bool create_model)
{
	if(create_model)
	{
		activeModel = mclr->CreateActiveLearningModel(pawTable);
	}

	Perform_Classification();

	//this->HeatmapforActivelearning(table, mclr->numberOfClasses);
}

void ALforNucEd::Perform_Classification()
{
	for(int i=0; i<classificationTables.size() ; ++i)
	{	
		//test_table contains individual tables 
		vtkSmartPointer<vtkTable> test_table  = vtkSmartPointer<vtkTable>::New();
		test_table->Initialize();

		test_table->SetNumberOfRows(classificationTables.at(i)->GetNumberOfRows());
		for(int col=0; col<pawTable->GetNumberOfColumns(); ++col)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName(pawTable->GetColumnName(col));
			test_table->AddColumn(column);	
		}
		for(int row = 0; row < (int)classificationTables.at(i)->GetNumberOfRows(); ++row)
		{		
			vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
			for(int c =0;c<(int)test_table->GetNumberOfColumns();++c)
				model_data1->InsertNextValue(classificationTables.at(i)->GetValueByName(row,test_table->GetColumnName(c)));
			test_table->InsertNextRow(model_data1);
		}	

		////// Final Data  to classify after the active training
		vnl_matrix<double> data_classify;
		if(from_model)
		{
			data_classify =  mclr->Normalize_Feature_Matrix_w(mclr->tableToMatrix_w(test_table), std_dev_vec, mean_vec);
		}
		else
		{
			data_classify =  mclr->Normalize_Feature_Matrix(mclr->tableToMatrix(test_table, mclr->id_time_val));
		}
		data_classify = data_classify.transpose();

		vnl_matrix<double> currprob;
		if(from_model)
		{
			currprob = mclr->Test_Current_Model_w(data_classify, act_learn_matrix);
		}
		else
		{
			currprob = mclr->Test_Current_Model(data_classify);
		}

		if(classification_name != "")
		{
			prediction_col_name = "prediction_active_" + classification_name;
			//prediction_col_name = "prediction_active";
			confidence_col_name = "confidence_" + classification_name;
		}

		//// Add the Prediction Column 
		std::vector< std::string > prediction_column_names = ftk::GetColumsWithString(prediction_col_name.c_str() , classificationTables.at(i) );
		if(prediction_column_names.size() == 0)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName(prediction_col_name.c_str());
			column->SetNumberOfValues( classificationTables.at(i)->GetNumberOfRows() );
			classificationTables.at(i)->AddColumn(column);
		}

		// Add the confidence column
		std::vector< std::string > confidence_column_names = ftk::GetColumsWithString(confidence_col_name.c_str() , classificationTables.at(i) );
		if(confidence_column_names.size() == 0)
		{
			vtkSmartPointer<vtkDoubleArray> column_confidence = vtkSmartPointer<vtkDoubleArray>::New();
			column_confidence->SetName(confidence_col_name.c_str());
			column_confidence->SetNumberOfValues( classificationTables.at(i)->GetNumberOfRows() );
			classificationTables.at(i)->AddColumn(column_confidence);
		}

		for(int row = 0; (int)row < classificationTables.at(i)->GetNumberOfRows(); ++row)  
		{
			vnl_vector<double> curr_col = currprob.get_column(row);
			classificationTables.at(i)->SetValueByName(row,confidence_col_name.c_str(), vtkVariant(curr_col(curr_col.arg_max())));
			if(curr_col(curr_col.arg_max()) > confidence_thresh) 
			{
				classificationTables.at(i)->SetValueByName(row,prediction_col_name.c_str(), vtkVariant(curr_col.arg_max()+1));						
			}
			else
			{
				classificationTables.at(i)->SetValueByName(row,prediction_col_name.c_str(), vtkVariant(0));
			}
		}
		//prediction_names = ftk::GetColumsWithString( prediction_col_name.c_str() , classificationTables.at(i) );
		//selection->clear();			
	}

	emit Classification_Done();

}


//***************************************************************************
//***********************************************************************************
//***********************************************************************************
// A dialog to get the confidence threshold and classification name:
//***********************************************************************************


ClassNameConfidenceDialog::ClassNameConfidenceDialog(QWidget *parent)
: QDialog(parent)
{
	classNameLabel = new QLabel("Enter Classification Name : ");
	class_name = new QLineEdit();
	class_name->setMinimumWidth(30);
	class_name->setFocusPolicy(Qt::StrongFocus);
	classNameLayout = new QHBoxLayout;
	classNameLayout->addWidget(classNameLabel);
	classNameLayout->addWidget(class_name);
	
	confidenceLabel = new QLabel("Specify Confidence Threshold %: ");
	conf_thresh = new QLineEdit();
	conf_thresh->setMinimumWidth(30);
	conf_thresh->setFocusPolicy(Qt::StrongFocus);
	confLayout = new QHBoxLayout;
	confLayout->addWidget(confidenceLabel);
	confLayout->addWidget(conf_thresh);
	
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);

	layout = new QVBoxLayout;
	layout->addLayout(classNameLayout);
	layout->addLayout(confLayout);
	layout->addLayout(bLayout);
	this->setLayout(layout);
	this->setWindowTitle(tr("Confidence Threshold"));

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}

std::string ClassNameConfidenceDialog::getClassName()
{
	std::string className;
	QString input = class_name->displayText();
	className = input.toStdString();
	return className;
}

double ClassNameConfidenceDialog::getConfThresh()
{
	double CT = 0.5;
	QString input = conf_thresh->displayText();
	if(input != "")
	{
		CT = input.toDouble()/100;
	}
	return CT;
}
