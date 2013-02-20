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
#include "GenericALDialog.h"

//Constructors:
GenericALDialog::GenericALDialog(vtkSmartPointer<vtkTable> table,int num_classes,int active_query,std::vector<int>top_feats,QWidget *parent)
: QDialog(parent)
{

	//temp_pair.first = snapshot;
		
	//QLabel *imageLabel = new QLabel(this);
	//imageLabel->setPixmap(QPixmap::fromImage(snapshot));
	//imageLabel->resize(imageLabel->pixmap()->size());
	
	

	this->setWindowTitle(tr("Active Learning Window: Specify Class"));
	this->setModal(false);


	finish = true;

	//Master Layout
	QGridLayout * layout = new QGridLayout;
	
	//// textbox
	//lineEdit = new QLineEdit();
	//lineEdit->setMaxLength(15);

	//Top-row of the window 
	QHBoxLayout *topRow = new QHBoxLayout;
	//topRow->addWidget(imageLabel,1,0);
	
	//Done Button
	doneButton = new QPushButton("Done");

	connect(doneButton, SIGNAL(clicked()), this, SLOT(finished()));
	doneButton->setDisabled(true);
	doneButton->setDefault(false);
	doneButton->setAutoDefault(false);

	nextButton = new QPushButton("Next");

	connect(nextButton, SIGNAL(clicked()), this, SLOT(accept()));
	nextButton->setDisabled(true);
	nextButton->setDefault(false);
	nextButton->setAutoDefault(true);
	
	QLabel *enterClassLabel = new QLabel("Select Class: ", this);
	topRow->addWidget(enterClassLabel,0,0);
	


	// radiobutons for each class
	for(int i =1 ; i<= num_classes ; ++i)
	{
		QRadioButton *class_button = new QRadioButton(QString::number(i), this);
		topRow->addWidget(class_button,0,0);
		button_vector.push_back(class_button);
		connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class()));
		// Class 1 is always the default selection
	}
		
		QRadioButton *class_button = new QRadioButton("I am not sure" , this);
		connect(class_button, SIGNAL(clicked()),this,SLOT(Set_Class()));
		button_vector.push_back(class_button);
		topRow->addWidget(class_button,0,0);
		topRow->addWidget(nextButton,0,0);
		topRow->addWidget(doneButton,0,0);
		topRow->addStretch(20);


	QLabel *channelLabel = new QLabel("Please ensure all the relevant channels which might affect classification are ON ", this);


	//Remove the train column 
	table->RemoveColumnByName("train");
	int numb_cols = MIN(5,table->GetNumberOfColumns());


	QStringList sl;

	for(int i=0;i<numb_cols;++i)
	{
		sl<<QString(table->GetColumnName(top_feats[i]));
	}
	
	// Table displaying 10 top features
	QTableWidget *tableWidget = new QTableWidget(1, numb_cols, this);
	tableWidget->setHorizontalHeaderLabels(sl);


	for(int i=0;i<numb_cols;++i)
	{
	  QTableWidgetItem *newItem = new QTableWidgetItem(QString::number(table->GetValue(active_query,top_feats[i]).ToDouble()));
	  newItem->setFlags(newItem->flags() & ~Qt::ItemIsEditable);
	  tableWidget->setItem(0, i, newItem);
	}
	
	class_selected = -1; // Used to check if no radiobutton was selected
	
	//Bottom-row of the window 
	QVBoxLayout *botRow = new QVBoxLayout;
	botRow->addWidget(tableWidget,1,0);
	botRow->addWidget(channelLabel,1,0);
	
	layout->addLayout(topRow,0,0,0);
	layout->addLayout(botRow,1,0,0);

	this->setLayout(layout);
}

void GenericALDialog::finished()
{
	finish = false;
	this->accept();
}

void GenericALDialog::Set_Class()
{
	nextButton->setDisabled(false);
	doneButton->setDisabled(false);
	for(int i =0; i< button_vector.size() ; ++i)
	{
		QRadioButton * test_button = button_vector.at(i);
		if(test_button->isChecked())
		{
			//std::cout<<test_button->text().toStdString() << " was selected" <<std::endl;
			if(test_button->text().toStdString()=="I am not sure")
				class_selected = 0;
			else
				class_selected = test_button->text().toInt();
		}
	}

	//temp_pair.second = class_selected;	

}


