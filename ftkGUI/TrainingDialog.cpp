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

TrainingDialog::TrainingDialog(vtkSmartPointer<vtkTable> table, const char * trainColumn, QWidget *parent)
	: QDialog(parent)
{
	this->setWindowTitle(tr("Training"));

	QLabel *topLabel = new QLabel(tr("Please enter comma separated list of IDs: "));

	inputsLayout = new QVBoxLayout;

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

	saveButton = new QPushButton(tr("Save Model"));
	connect(saveButton, SIGNAL(clicked()), this, SLOT(saveModel()));
	saveButton->setDefault(false);
	saveButton->setAutoDefault(false);
	quitButton = new QPushButton(tr("Cancel"));
	connect(quitButton, SIGNAL(clicked()), this, SLOT(reject()));
	quitButton->setDefault(false);
	quitButton->setAutoDefault(false);
	doneButton = new QPushButton(tr("Done"));
	connect(doneButton, SIGNAL(clicked()), this, SLOT(accept()));
	doneButton->setDefault(false);
	doneButton->setAutoDefault(false);

	QHBoxLayout *endLayout = new QHBoxLayout;
	endLayout->addStretch(20);
	endLayout->addWidget(saveButton);
	endLayout->addWidget(quitButton);
	endLayout->addWidget(doneButton);

	QVBoxLayout *masterLayout = new QVBoxLayout;
	masterLayout->addWidget(topLabel);
	masterLayout->addLayout(inputsLayout);
	masterLayout->addLayout(bLayout);
	masterLayout->addStretch(20);
	masterLayout->addLayout(endLayout);

	this->setLayout(masterLayout);

	this->addClass();

	m_table = table;
	columnForTraining = trainColumn;
}

void TrainingDialog::addClass(void)
{
	if(inputLabels.size() == 10)
		return;

	//Create input box
	QLabel *label = new QLabel("Class " + QString::number(inputValues.size() + 1) + ": ");
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

void TrainingDialog::saveModel(void)
{
	QMessageBox::critical(this, tr("Oops"), tr("This function not yet implemented"));
}

void TrainingDialog::accept(void)
{
	//Update the table:
	this->parseInputValues();
	this->updateTable();

	//Exit:
	QDialog::accept();
}

void TrainingDialog::parseInputValues(void)
{
	training.clear();

	for(int c=0; c<inputValues.size(); ++c)
	{
		std::set<int> ids;
		QString input = inputValues.at(c)->displayText();
		QStringList values = input.split(",");
		for(int i=0; i<values.size(); ++i)
		{
			QString str = values.at(i);
			int v = str.toInt();
			ids.insert( v );
		}
		training.push_back(ids);
	}
}

void TrainingDialog::updateTable(void)
{
	if(training.size() == 0)
		return;

	//If need to create a new column do so now:
	vtkAbstractArray * output = m_table->GetColumnByName(columnForTraining);
	if(output == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( columnForTraining );
		column->SetNumberOfValues( m_table->GetNumberOfRows() );
		column->FillComponent(0,-1);
		m_table->AddColumn(column);
	}

	for(int row = 0; (int)row < m_table->GetNumberOfRows(); ++row)  
	{
		int id = m_table->GetValue(row,0).ToInt();
		for(int c=0; c<training.size(); ++c)
		{
			if( training.at(c).find(id) != training.at(c).end() )
			{
				m_table->SetValueByName(row, columnForTraining, vtkVariant( c + 1 ));
				break;
			}
		}
	}

	emit changedTable();
}