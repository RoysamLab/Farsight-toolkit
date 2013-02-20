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
#include "Active_Learning_Dialog.h"
#include "math.h"

//Constructors:
Active_Learning_Dialog::Active_Learning_Dialog(std::vector<QImage> snapshot, vtkSmartPointer<vtkTable> table,int num_classes,std::vector<int> active_query,std::vector<int>top_feats,QWidget *parent)
: QDialog(parent)
{
	rejectFlag = false;
	queries = active_query;
	num_of_classes = num_classes;
	query_label.resize(snapshot.size());
	buttongroup.resize(snapshot.size());

	this->setWindowTitle(tr("ALAMO Window: Specify Class"));
	this->setModal(true);
	//Master Layout
	QGridLayout * layout = new QGridLayout;


	for(int i=0;i<snapshot.size();++i)
	{
	std::vector<QHBoxLayout *> rows  = Sample_Details(snapshot[i],table,num_classes,active_query[i],top_feats,i);
	layout->addLayout(rows[0],2*i,0,0);
	layout->addLayout(rows[1],2*i+1,0,0);
	}
	
	//Sample_Details(snapshot,table,num_classes,active_query,top_feats);
	////QLabel *channelLabel = new QLabel("Please ensure all the relevant channels which might affect classification are ON ", this);	
	////Master Layout
	//layout->addLayout(topRow,2,0,0);
	//layout->addLayout(botRow,3,0,0);


	//Done Button
	QPushButton *doneButton = new QPushButton("Done");
	connect(doneButton, SIGNAL(clicked()), this, SLOT(finished()));
	doneButton->setDefault(false);
	doneButton->setAutoDefault(false);


	QPushButton *nextButton = new QPushButton("Next");
	connect(nextButton, SIGNAL(clicked()), this, SLOT(accept()));
	nextButton->setDefault(false);
	nextButton->setAutoDefault(true);


	//Cancel Button
	QPushButton *cancelButton = new QPushButton("Cancel");
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(rejectAction()));
	cancelButton->setDefault(false);
	cancelButton->setAutoDefault(false);

	//Top-row of the window 
	QHBoxLayout *finalRow = new QHBoxLayout;
	finalRow->addWidget(nextButton,0,0);
	finalRow->addWidget(doneButton,0,0);
	finalRow->addWidget(cancelButton,0,0);

	layout->addLayout(finalRow,2*snapshot.size()+2,0,0);
	this->setLayout(layout);
	
	
	for(int i =0 ; i<query_label.size() ; ++i)
	{
		query_label[i].second = -1;// Used to check if no radiobutton was selected for any sample
	}

	//QLabel *channelLabel = new QLabel("Please ensure all the relevant channels which might affect classification are ON ", this);	
	this->resize(600,600);
}



//Constructors:
Active_Learning_Dialog::Active_Learning_Dialog(std::string validate,std::vector<QImage> snapshot, vtkSmartPointer<vtkTable> table,int classval,std::vector<int> rowvals,int num_classes,QWidget *parent)
: QDialog(parent)
{
	this->setWindowTitle(tr("Validation Window: Validate Classifier Performance"));
	this->setModal(true);

	rejectFlag = false;
	
	queries = rowvals;
	num_of_classes = num_classes;

	query_label.resize(rowvals.size());
	buttongroup.resize(rowvals.size());

	//Master Layout
	QGridLayout * layout = new QGridLayout;
	
	QHBoxLayout *topRow = new QHBoxLayout;
	std::vector<int>::iterator it; 
	
	//set font    
	QFont font;
	font.setPointSize(15);
	font.setBold(true);
	
	QString x = "All these samples have been classified as class  "+QString::number(classval);
	QLabel *ClassLabel = new QLabel(x, this);
	ClassLabel->setFont(font);
	layout->addWidget(ClassLabel,0,0);


	//Shows 5 samples at a time
	for(int i=0;i<snapshot.size();++i)
	{
		std::vector<QHBoxLayout *> rows  = Validation_Sample_Details(snapshot[i],table,classval,rowvals[i],i,num_classes);
		layout->addLayout(rows[0],2*i+1,0,0);
		layout->addLayout(rows[1],2*i+2,0,0);
	}

	//Done Button
	QPushButton *doneButton = new QPushButton("Done");
	connect(doneButton, SIGNAL(clicked()), this, SLOT(finished()));
	doneButton->setDefault(false);
	doneButton->setAutoDefault(false);


	//Cancel Button
	QPushButton *cancelButton = new QPushButton("Cancel");
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(rejectAction()));
	cancelButton->setDefault(false);
	cancelButton->setAutoDefault(false);


	//QPushButton *nextButton = new QPushButton("Next");
	//connect(nextButton, SIGNAL(clicked()), this, SLOT(accept()));
	//nextButton->setDefault(false);
	//nextButton->setAutoDefault(true);

	//Top-row of the window 
	QHBoxLayout *finalRow = new QHBoxLayout;
	//finalRow->addWidget(nextButton,0,0);
	finalRow->addWidget(doneButton,0,0);
	finalRow->addWidget(cancelButton,0,0);

	layout->addLayout(finalRow,14,0,0);
	this->setLayout(layout);

	for(int i =0 ; i<query_label.size() ; ++i)
	{
		query_label[i].second = -1;// Used to check if no radiobutton was selected for any sample
	}

	//QLabel *channelLabel = new QLabel("Please ensure all the relevant channels which might affect classification are ON ", this);	
	this->resize(600,600);

}


std::vector<QHBoxLayout *> Active_Learning_Dialog::Sample_Details(QImage snapshot, vtkSmartPointer<vtkTable> table,int num_classes,int active_query,std::vector<int>top_feats,int group)
{

	QLabel *imageLabel = new QLabel(this);
	imageLabel->setPixmap(QPixmap::fromImage(snapshot));
	imageLabel->resize(imageLabel->pixmap()->size());
	
	finish = true;

	
	//Top-row of the window 
	QHBoxLayout *topRow = new QHBoxLayout;
	topRow->addWidget(imageLabel,0,0);
	
	QLabel *enterClassLabel = new QLabel("Select Class: ", this);
	topRow->addWidget(enterClassLabel,0,0);
	
	buttongroup[group] = new QButtonGroup(this);
	buttongroup[group]->setExclusive(true);

	// radiobutons for each class
	for(int i =0 ; i< num_classes ; ++i)
	{
		QRadioButton  *class_button = new QRadioButton(QString::number(i+1), this);
		topRow->addWidget(class_button,0,0);
		button_vector.push_back(class_button);
		buttongroup[group]->addButton(class_button);
		connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class()));
	}
			
	QRadioButton *class_button = new QRadioButton("I am not sure" , this);
	connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class()));
	button_vector.push_back(class_button);
	buttongroup[group]->addButton(class_button);
	topRow->addWidget(class_button,0,0);

	
	//Remove the train column 
	table->RemoveColumnByName("train");
	int numb_cols = MIN(5,table->GetNumberOfColumns());
 
	QStringList sl;

	for(int i=0;i<numb_cols;++i)
	{
		sl<<QString(table->GetColumnName(top_feats[i]));
	}
	
	// Table displaying 5 top features
	QTableWidget *tableWidget = new QTableWidget(1, numb_cols, this);
	tableWidget->setHorizontalHeaderLabels(sl);
	tableWidget->setRowCount(1);
    tableWidget->setColumnCount(numb_cols);

	for(int i=0;i<numb_cols;++i)
	{
	  QTableWidgetItem *newItem = new QTableWidgetItem(QString::number(table->GetValue(active_query,top_feats[i]).ToDouble()));
	  newItem->setFlags(newItem->flags() & ~Qt::ItemIsEditable);
	  tableWidget->setItem(0, i, newItem);
	}

	//classes_chosen.resize(num_classes+1);

	//for(int i =0 ; i<=num_classes ; ++i)
	//{
	//	classes_chosen[i] = -1;// Used to check if no radiobutton was selected for any sample
	//}

	//Bottom-row of the window 
	QHBoxLayout *botRow = new QHBoxLayout;
	botRow->addWidget(tableWidget,0,0);
	//botRow->addWidget(channelLabel,1,0);
	std::vector<QHBoxLayout *> rows;
	rows.push_back(topRow);
	rows.push_back(botRow);

	return rows;
}


std::vector<QHBoxLayout *> Active_Learning_Dialog::Validation_Sample_Details(QImage snapshot, vtkSmartPointer<vtkTable> table,int classval,int PIA_query,int group,int num_classes)
{	

	QLabel *imageLabel = new QLabel(this);
	imageLabel->setPixmap(QPixmap::fromImage(snapshot));
	imageLabel->resize(imageLabel->pixmap()->size());
	
	finish = true;
	
	//Top-row of the window 
	QHBoxLayout *topRow = new QHBoxLayout;
	topRow->addWidget(imageLabel,0,0);
	
	QLabel *enterClassLabel = new QLabel("Please indicate the correct class: ", this);
	topRow->addWidget(enterClassLabel,0,0);
	
	buttongroup[group] = new QButtonGroup(this);
	buttongroup[group]->setExclusive(true);

	// radiobutons for each class
	for(int i =0 ; i< num_classes ; ++i)
	{
		QRadioButton  *class_button = new QRadioButton(QString::number(i+1), this);
		topRow->addWidget(class_button,0,0);
		button_vector.push_back(class_button);
		buttongroup[group]->addButton(class_button);
		connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class_Validate()));
	}
			
	//QRadioButton *class_button = new QRadioButton("I am not sure" , this);
	//connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class()));
	//button_vector.push_back(class_button);
	//buttongroup[group]->addButton(class_button);
	//topRow->addWidget(class_button,0,0);


	QStringList sl;
	for(int i=0;i<table->GetNumberOfColumns();++i)
	{
		sl<<QString(table->GetColumnName(i));
	}

	// Table displaying features of the sample
	QTableWidget *tableWidget = new QTableWidget(1, table->GetNumberOfColumns(), this);
	tableWidget->setHorizontalHeaderLabels(sl);
	tableWidget->setRowCount(1);
    tableWidget->setColumnCount(table->GetNumberOfColumns());

	for(int i=0;i<table->GetNumberOfColumns();++i)
	{
	  QTableWidgetItem *newItem = new QTableWidgetItem(QString::number(table->GetValue(PIA_query,i).ToDouble()));
	  newItem->setFlags(newItem->flags() & ~Qt::ItemIsEditable);
	  tableWidget->setItem(0, i, newItem);
	}

	//Bottom-row of the window 
	QHBoxLayout *botRow = new QHBoxLayout;
	botRow->addWidget(tableWidget,0,0);
	//botRow->addWidget(channelLabel,1,0);
	std::vector<QHBoxLayout *> rows;
	rows.push_back(topRow);
	rows.push_back(botRow);

	return rows;
}


void Active_Learning_Dialog::finished()
{
	finish = false;
	this->accept();
}


void Active_Learning_Dialog::rejectAction()
{
	this->rejectFlag = true;
	this->reject();
}



void Active_Learning_Dialog::Set_Class()
{
	for(int i =0; i< button_vector.size() ; ++i)
	{
		QRadioButton * test_button = button_vector.at(i);
		if(test_button->isChecked())
		{
			if(test_button->text().toStdString()=="I am not sure")
			{
				query_label.at(i/(num_of_classes+1)).first = queries[i/(num_of_classes+1)];
				query_label.at(i/(num_of_classes+1)).second = 0;
			}
			else
			{	
				query_label.at(i/(num_of_classes+1)).first = queries[i/(num_of_classes+1)];
				query_label.at(i/(num_of_classes+1)).second = test_button->text().toInt();
			}
		}
	}
}

void Active_Learning_Dialog::Set_Class_Validate()
{
	for(int i =0; i< button_vector.size() ; ++i)
	{
		QRadioButton * test_button = button_vector.at(i);
		if(test_button->isChecked())
		{
				query_label.at(i/(num_of_classes)).first = queries[i/(num_of_classes)];
				query_label.at(i/(num_of_classes)).second = test_button->text().toInt();
		}
	}
}



void Active_Learning_Dialog::User_Validation_Response()
{
	for(int i =0; i< button_vector.size() ; ++i)
	{
		QRadioButton * test_button = button_vector.at(i);
		if(test_button->isChecked())
		{
			query_label.at(i/2).first =  -1;//dummy
			query_label.at(i/2).second = i%2;
		}
	}
}

