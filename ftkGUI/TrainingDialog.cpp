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

#include "TrainingDialog.h"

TrainingDialog::TrainingDialog(vtkSmartPointer<vtkTable> table, const char * trainColumn,std::string mode,int firstTableSize, QWidget *parent)
	: QDialog(parent){

	m_table = table;
	columnForTraining = trainColumn;
	firstTableRows = firstTableSize;

	this->GetTrainingNames( m_table );

	this->setWindowTitle(tr("Training"));
	this->setModal(false);

	QLabel *topLabel = new QLabel(tr("Please enter a comma separated list of Example Object IDs: "));
	QLabel *topLabel1 = new QLabel(tr("Select the classifier to be trained or enter a new classifier: "));


	classesLayout = new QVBoxLayout;
	inputsLayout = new QVBoxLayout;

	classesLayout->addLayout( createFirstExclusiveGroup() );

	addButton = new QPushButton(tr("Add Class"));
	connect(addButton, SIGNAL(clicked()), this, SLOT(addClass()));
	addButton->setDefault(false);
	addButton->setAutoDefault(false);
	delButton = new QPushButton(tr("Remove Class"));
	connect(delButton, SIGNAL(clicked()), this, SLOT(remClass()));
	delButton->setDefault(false);
	delButton->setAutoDefault(false);
	QHBoxLayout *bLayout = new QHBoxLayout;
	bLayout->addWidget(addButton);
	bLayout->addWidget(delButton);
	bLayout->addStretch(20);

	saveButton = new QPushButton(tr("Save Training Set"));
	connect(saveButton, SIGNAL(clicked()), this, SLOT(saveModel(void)));
	saveButton->setDefault(false);
	saveButton->setAutoDefault(false);

	loadButton = new QPushButton(tr("Load Training Set"));
	connect(loadButton, SIGNAL(clicked()), this, SLOT(loadModel()));
	loadButton->setDefault(false);
	loadButton->setAutoDefault(false);

	quitButton = new QPushButton(tr("Cancel"));
	connect(quitButton, SIGNAL(clicked()), this, SLOT(reject()));
	quitButton->setDefault(false);
	quitButton->setAutoDefault(false);
	doneButton = new QPushButton(tr("Update Selected"));
	connect(doneButton, SIGNAL(clicked()), this, SLOT(accept()));
	doneButton->setDefault(false);
	doneButton->setAutoDefault(false);

	QHBoxLayout *endLayout = new QHBoxLayout;
	endLayout->addStretch(25);
	endLayout->addWidget(loadButton);
	endLayout->addWidget(saveButton);
	endLayout->addWidget(quitButton);
	endLayout->addWidget(doneButton);

	QVBoxLayout *masterLayout = new QVBoxLayout;
	masterLayout->addWidget(topLabel1);
	masterLayout->addLayout(classesLayout);
	masterLayout->addWidget(topLabel);
	masterLayout->addLayout(inputsLayout);
	masterLayout->addLayout(bLayout);
	masterLayout->addStretch(20);
	masterLayout->addLayout(endLayout);


	if(mode=="active")
	{
		topLabel->setText("Please enter one example for each class");
		topLabel1->hide();
		radio6->hide();
		className->hide();
		loadButton->setEnabled(FALSE);
		saveButton->setEnabled(FALSE);
	}

	this->setLayout(masterLayout);

	//this->addClass();
	if( class_names.empty() )
		this->differentclassselected(6);
	else
		this->differentclassselected(1);

	use_train_gui = true;

	this->addClass();
}

void TrainingDialog::addClass(void)
{
	if(inputLabels.size() == 10)
		return;

	LabelImageViewQT *labelimview = new LabelImageViewQT;
	QString class_color_name = labelimview->GetColorNameFromTable( inputValues.size()+1 );
	delete labelimview;
	//Create input box
	QLabel *label = new QLabel("Class " + QString::number(inputValues.size() + 1) + " ("  + class_color_name + ")" + ": ");
	inputLabels.push_back( label );

	QLineEdit *inVals = new QLineEdit();
	inVals->setMinimumWidth(200);
	inVals->setFocusPolicy(Qt::StrongFocus);
	inputValues.push_back( inVals );

	QHBoxLayout *layout = new QHBoxLayout;
	layout->addWidget( inputLabels.back() );
	layout->addWidget( inputValues.back() );
	layout->addStretch(1);
	iLayouts.push_back(layout);
	inputsLayout->addLayout( iLayouts.back() );

	for(int i=1; i<inputValues.size(); ++i)
		QWidget::setTabOrder(inputValues.at(i-1), inputValues.at(i));
	QWidget::setTabOrder(inputValues.back(), addButton);
	QWidget::setTabOrder(addButton,delButton);
	QWidget::setTabOrder(delButton,quitButton);
	QWidget::setTabOrder(quitButton,saveButton);
	QWidget::setTabOrder(saveButton, inputValues.front());

	inputValues.back()->setFocus();

}

void TrainingDialog::remClass(void)
{
	if(inputLabels.size() == 0)
		return;

	delete inputLabels.back();
	delete inputValues.back();

	inputLabels.remove(inputLabels.size()-1);
	inputValues.remove(inputValues.size()-1);

	delete iLayouts.back();
	iLayouts.remove(iLayouts.size()-1);

}


void TrainingDialog::saveModel(void){
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Model As..."),lastPath, tr("TEXT(*.xml)"));
	if(filename == "")
		return;
	lastPath = QFileInfo(filename).absolutePath();

	this->accept();

	if(training.size()<2){
		QMessageBox::critical(this, tr("Oops"), tr("Please enter ids for atleast 2 classes to save the model"));
		this->show();
		return;
	}

	//Load the model into a new table 
	vtkSmartPointer<vtkTable> new_table = vtkSmartPointer<vtkTable>::New();
	new_table->Initialize();

	vtkSmartPointer<vtkVariantArray> model_data = vtkSmartPointer<vtkVariantArray>::New();
	model_data = m_table->GetRow(1);
	//Make a copy of the table and delete the prediction columns if present
	vtkSmartPointer<vtkTable> m_table_cpy = vtkSmartPointer<vtkTable>::New();
	for(int i =0;i<model_data->GetNumberOfValues();++i){
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( m_table->GetColumnName(i) );
		m_table_cpy->AddColumn(column);
	}
	for(int row = 0; row < (int)m_table->GetNumberOfRows(); ++row){
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int i =0;i<m_table->GetNumberOfColumns();++i){
			model_data1->InsertNextValue(m_table->GetValue(row,i));
		}
		m_table_cpy->InsertNextRow(model_data1);
	}

	for( int i=0; i<m_table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = m_table->GetColumnName(i);
		if( current_column.find("prediction") != std::string::npos )
			m_table_cpy->RemoveColumnByName( current_column.c_str() );
	}

	//Add columns to the model table
	vtkSmartPointer<vtkVariantArray> model_data2 = vtkSmartPointer<vtkVariantArray>::New();
	model_data2 = m_table_cpy->GetRow(1);
	for(int i =0;i<model_data2->GetNumberOfValues();++i){
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		std::string current_column;
		current_column = m_table_cpy->GetColumnName(i);
		column->SetName( m_table_cpy->GetColumnName(i) );
		new_table->AddColumn(column);
	}

	//Add only the rows which have been selected by the user and do not repeat ids
	for( int j=0; j<(int)training_names.size(); ++j ){
		for(int row = 0; row < (int)m_table_cpy->GetNumberOfRows(); ++row){
			vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
			if(m_table_cpy->GetValueByName(row, training_names.at(j).c_str())!=-1){
				if( j == 1 )
					if(m_table_cpy->GetValueByName(row, training_names.at(j-1).c_str())!=-1)	continue;
				if( j == 2 )
					if(m_table_cpy->GetValueByName(row, training_names.at(j-1).c_str())!=-1 ||
					   m_table_cpy->GetValueByName(row, training_names.at(j-2).c_str())!=-1 )	continue;
				for(int i =0;i<m_table_cpy->GetNumberOfColumns();++i){
					model_data1->InsertNextValue(m_table_cpy->GetValue(row,i));
				}
				new_table->InsertNextRow(model_data1);
			}
		}
	}

	//Write the model into an XML file.	
	bool ok = ftk::SaveTable(filename.toStdString(),new_table);
	if(!ok){
		cerr << "problem writing model to " << filename.toStdString() << endl;
	}
}

void TrainingDialog::loadModel(void){
	QString fileName  = QFileDialog::getOpenFileName(this, "Select training model to open", lastPath,
									tr("TXT Files (*.xml)"));
	if(fileName == "")
		return;
	lastPath = QFileInfo(fileName).absolutePath();

	this->loadModelFromFile(fileName.toStdString());

	//Exit:
	QDialog::accept();
}

void TrainingDialog::loadModelFromFile( std::string file_name ){
	model_table = ftk::LoadTable(file_name);
	if(!model_table) return;
	//Append the data to the current table
	this->GetTrainingNames( model_table );
	this->Append();
	return;
}

//Update the features in this table whose names match (sets doFeat)
//Will creat new rows if necessary:
void TrainingDialog::Append(){
	for( int i=0; i<(int)training_names.size(); ++i ){
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( training_names.at(i).c_str() );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		column->FillComponent( 0,-1 );
		m_table->AddColumn(column);
	}

	int maxrows = m_table->GetNumberOfRows();
	int maxrowid = m_table->GetValueByName(maxrows-1,"ID").ToInt();

	vtkSmartPointer<vtkVariantArray> model_data = vtkSmartPointer<vtkVariantArray>::New();

	for(int row = 0; (int)row < model_table->GetNumberOfRows(); ++row){
		model_data = model_table->GetRow(row);
		m_table->InsertNextRow(model_data); 
		m_table->SetValue(maxrows+row,0,maxrowid+row+1);
	}
		
		//vtkSmartPointer<vtkStringArray> discrim_column = vtkSmartPointer<vtkStringArray>::New();
		//discrim_column->SetName("Original/Model" );
		//discrim_column->SetNumberOfValues(m_table->GetNumberOfRows());
		//m_table->AddColumn(discrim_column);
		//
		//for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  
		//	{	
		//		std::cout<<"I was here"<<std::endl;
		//		if(row<maxrow2)
		//		{
		//			m_table->SetValueByName(row, "Original/Model", "Original");
		//			std::cout<<"I was here"<<std::endl;
		//		}
		//		if(row>=maxrow2)
		//		{
		//			m_table->SetValueByName(row, "Original/Model", "Model");
		//			std::cout<<"I was there"<<std::endl;
		//		}
		//	}

	 if( use_train_gui ) emit changedTable();
	 else return;
}



void TrainingDialog::accept(void){
	//Update the table:
	this->inputToTable();

	//Exit:
	QDialog::accept();
}

void TrainingDialog::inputToTable(void){
	this->parseInputValues();	//inputs to training set
	this->updateTable();		//training set to table
}

void TrainingDialog::tableToInput(void){
	this->parseTableValues();	//table to training set
	this->updateInputs();		//training set to input boxes
}

//Takes the values from the input boxes and puts them in the training vector.
void TrainingDialog::parseInputValues(void){
	training.clear();

	for(int c=0; c<inputValues.size(); ++c){
		std::set<int> ids;
		QString input = inputValues.at(c)->displayText();
		QStringList values = input.split(",");
		for(int i=0; i<values.size(); ++i){
			QString str = values.at(i);
			int v = str.toInt();
			ids.insert( v );
		}
		training[ c+1 ] = ids;
	}
}

//Takes the value from the columnForTraining in the table and puts them into the training vector
void TrainingDialog::parseTableValues(void)
{
	training.clear();

	if(!m_table) return;

	vtkAbstractArray * output = m_table->GetColumnByName(columnForTraining);
	if(output == 0) return;

	//Iterate through table and populate the training vectors
	for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  
	{
		int cls = m_table->GetValueByName(row, columnForTraining).ToInt();
		int id = m_table->GetValue(row,0).ToInt();
		if(cls > 0)
			training[ cls ].insert( id);
	}
}

// Go from training set to table;
void TrainingDialog::updateTable(void){
	std::map< int, std::set<int> >::iterator it;

	//if(training.size() == 0)
	//	return;

	if(!m_table) return;

	if( class_names.size() < 5 ){
		if( radio6->isChecked() ){
			new_class_name.clear();
			QString input_string = className->displayText();
			std::string class_name_prefix;
			class_name_prefix = input_string.toStdString();
			if( !class_name_prefix.empty() ){
				new_class_name = "train_" + class_name_prefix;
				columnForTraining = new_class_name.c_str();
			}
		}
	}

	//If need to create a new column do so now:
	vtkAbstractArray * output = m_table->GetColumnByName(columnForTraining);
	if(output == 0){
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForTraining );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		column->FillComponent(0,-1);
		m_table->AddColumn(column);
	}

	for(int row = 0; (int)row < firstTableRows; ++row){
		int id = m_table->GetValue(row,0).ToInt();
		bool idFound = false;
		for(it=training.begin(); it!=training.end(); ++it){
			int cls = (*it).first;
			if( (*it).second.find(id) != (*it).second.end() ){
				m_table->SetValueByName(row, columnForTraining, vtkVariant( cls ));
				idFound = true;
				break;
			}
		}
		
		std::vector<std::string> time_name = ftk::GetColumsWithString( "time" , m_table);
		int time_col = 0;
		if(time_name.size()>0)
		{
			 time_col = m_table->GetValueByName(row,"time").ToInt();
		}

		if(!idFound || time_col!=0)
			m_table->SetValueByName(row, columnForTraining, vtkVariant( -1 ));
	}
	emit changedTable();
}

//Go from training set to input boxes
void TrainingDialog::updateInputs(void){
	std::map< int, std::set<int> >::iterator it;

	//Make sure number of classes matches:
	int numClassBoxes = (int)inputValues.size();
	int numClasses = (int)training.size();

	if(numClasses > numClassBoxes)	//Need more boxes
	{
		for(int i=numClassBoxes; i<numClasses; ++i){
			this->addClass();
		}
	}
	else if(numClassBoxes > numClasses) //Need to remove boxes
	{
		for(int i=numClasses; i<numClassBoxes; ++i){
			this->remClass();
		}
	}

	//Now fill in the boxes:
	int c=0;
	for(it=training.begin(); it!=training.end(); ++it){
		//int cls = (*it).first;
		std::set<int> ids = (*it).second;
		std::set<int>::iterator s_it = ids.begin();
		QString input;
		for(s_it=ids.begin(); s_it!=ids.end(); ++s_it){
			if(s_it != ids.begin())
				input.append(",");
			input.append(QString::number(*s_it));
		}
		inputValues.at(c)->setText(input);
		++c;
	}
}

void TrainingDialog::GetTrainingNames( vtkSmartPointer<vtkTable> table ){
	//getcoulmn names for m_table and look for strings starting with train_
	training_names.clear();
	class_names.clear();
	for( int i=0; i<table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = table->GetColumnName(i);
		if( current_column.find("train") != std::string::npos ){
			training_names.push_back( current_column );
			std::string::iterator it;
			it=current_column.begin();
			current_column.erase ( current_column.begin(), current_column.begin()+6 );
			class_names.push_back( current_column );
		}
	}
}

QHBoxLayout *TrainingDialog::createFirstExclusiveGroup(){
	QHBoxLayout *vbox = new QHBoxLayout;//(tr("Select Classifier to Train"));
	if( !class_names.empty() ){
		radio1 = new QRadioButton(tr( class_names.at(0).c_str() ));
		connect(radio1, SIGNAL( pressed() ), this, SLOT(class_1_selected()));
		vbox->addWidget(radio1);
		if( class_names.size() > 1 ){
			radio2 = new QRadioButton(tr( class_names.at(1).c_str() ));
			connect(radio2, SIGNAL( pressed() ), this, SLOT(class_2_selected()));
			vbox->addWidget(radio2);
			if( class_names.size() > 2 ){
				radio3 = new QRadioButton(tr( class_names.at(2).c_str() ));
				connect(radio3, SIGNAL( pressed() ), this, SLOT(class_3_selected()));
				vbox->addWidget(radio3);
				if( class_names.size() > 3 ){
					radio4 = new QRadioButton(tr( class_names.at(3).c_str() ));
					connect(radio4, SIGNAL( pressed() ), this, SLOT(class_4_selected()));
					vbox->addWidget(radio4);
					if( class_names.size() > 4 ){
						radio5 = new QRadioButton(tr( class_names.at(4).c_str() ));
						connect(radio5, SIGNAL( pressed() ), this, SLOT(class_5_selected()));
						vbox->addWidget(radio5);
					}
				}
			}
		}
	}
	if( class_names.empty() || class_names.size() < 5 ){
		radio6 = new QRadioButton(tr( "New Classifier:" ));
		connect(radio6, SIGNAL( pressed() ), this, SLOT(class_6_selected()));
		vbox->addWidget(radio6);
		className = new QLineEdit();
		className->setMinimumWidth(100);
		className->setFocusPolicy(Qt::StrongFocus);
		//QHBoxLayout *layout = new QHBoxLayout;
		vbox->addWidget( className );
		vbox->addStretch(1);
	}
	vbox->addStretch(1);

	return vbox;
}

void TrainingDialog::class_1_selected(void){
	this->differentclassselected( 1 );
}

void TrainingDialog::class_2_selected(void){
	this->differentclassselected( 2 );
}

void TrainingDialog::class_3_selected(void){
	this->differentclassselected( 3 );
}

void TrainingDialog::class_4_selected(void){
	this->differentclassselected( 4 );
}

void TrainingDialog::class_5_selected(void){
	this->differentclassselected( 5 );
}

void TrainingDialog::class_6_selected(void){
	this->differentclassselected( 6 );
}

void TrainingDialog::differentclassselected(int selected_class){
	if( !class_names.empty() ){
		if( selected_class==1 && !radio1->isChecked() ){
			selected_class = 1; //Set training mode to this one
			columnForTraining = training_names.at(0).c_str();
			radio1->setChecked(true);
		}
		if( class_names.size() > 1 ){
			if( selected_class==2 && !radio2->isChecked() ){
				selected_class = 2; //Set training mode to this one
				columnForTraining = training_names.at(1).c_str();
				radio2->setChecked(true);
			}
		}
		if( class_names.size() > 2 ){
				if( selected_class==3 && !radio3->isChecked() ){
					selected_class = 3; //Set training mode to this one
					columnForTraining = training_names.at(2).c_str();
					radio3->setChecked(true);
				}
		}
		if( class_names.size() > 3 ){
				if( selected_class==4 && !radio4->isChecked() ){
					selected_class = 4; //Set training mode to this one
					columnForTraining = training_names.at(3).c_str();
					radio4->setChecked(true);
				}
		}
		if( class_names.size() > 4 ){
				if( selected_class==5 && !radio5->isChecked() ){
					selected_class = 5; //Set training mode to this one
					columnForTraining = training_names.at(4).c_str();
					radio5->setChecked(true);
				}
		}
	}
	if( class_names.empty() || class_names.size() < 5 ){
		if( selected_class==6 && !radio6->isChecked() ){
			selected_class = 6; //Set training mode to this one
			std::ostringstream StrStream;
			if( class_names.empty() )
					StrStream << 1;
				else if( class_names.size()==1 )
					StrStream << 2;
				else if( class_names.size()==2 )
					StrStream << 3;
				else if( class_names.size()==3 )
					StrStream << 4;
				else if( class_names.size()==4 )
					StrStream << 5;
			default_training_name.clear();
			default_training_name = "train_default";
			default_training_name.append( StrStream.str() );
			columnForTraining = default_training_name.c_str();
			radio6->setChecked(true);
		}
	}
	this->tableToInput();
	return;
}
TrainingDialogNoGUI::TrainingDialogNoGUI(vtkSmartPointer<vtkTable> table){
	//This method has been created for batch processing mode and should be avoided when you expect a trainingdialog box
	m_table = table;
}

void TrainingDialogNoGUI::loadModelFromFile1( std::string file_name ){
	model_table = ftk::LoadTable(file_name);
	if(!model_table) return;
	//Append the data to the current table
	this->GetTrainingNames1( model_table );
	this->Append1();
	return;
}

//Update the features in this table whose names match (sets doFeat)
//Will creat new rows if necessary:
void TrainingDialogNoGUI::Append1(){
	for( int i=0; i<(int)training_names.size(); ++i ){
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( training_names.at(i).c_str() );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		column->FillComponent( 0,-1 );
		m_table->AddColumn(column);
	}

	int maxrows = m_table->GetNumberOfRows();
	int maxrowid = m_table->GetValueByName(maxrows-1,"ID").ToInt();

	vtkSmartPointer<vtkVariantArray> model_data = vtkSmartPointer<vtkVariantArray>::New();

	for(int row = 0; (int)row < model_table->GetNumberOfRows(); ++row){
		model_data = model_table->GetRow(row);
		m_table->InsertNextRow(model_data); 
		m_table->SetValue(maxrows+row,0,maxrowid+row+1);
	}

	return;
}

void TrainingDialogNoGUI::GetTrainingNames1( vtkSmartPointer<vtkTable> table ){
	//getcoulmn names for m_table and look for strings starting with train_
	training_names.clear();
	class_names.clear();
	for( int i=0; i<table->GetNumberOfColumns(); ++i ){
		std::string current_column;
		current_column = table->GetColumnName(i);
		if( current_column.find("train") != std::string::npos ){
			training_names.push_back( current_column );
			std::string::iterator it;
			it=current_column.begin();
			current_column.erase ( current_column.begin(), current_column.begin()+6 );
			class_names.push_back( current_column );
		}
	}
}
