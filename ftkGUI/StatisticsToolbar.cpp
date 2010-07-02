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
#include "StatisticsToolbar.h"

StatisticsToolbar::StatisticsToolbar(QWidget *parent)
	:QWidget(parent)
{
	

	this->inputDataTable = vtkSmartPointer<vtkTable>::New();
	this->StatisticsTable = new QTableView();
	this->StatisticsModel = new QStandardItemModel(this);
	this->StatisticsSelectionModel = new QItemSelectionModel(this->StatisticsModel, this);	
    this->statisticsDockWidget = new QDockWidget(tr("Statistics Toolbar"));
	
	

}


StatisticsToolbar::~StatisticsToolbar()
{
}



void StatisticsToolbar::setTable(vtkSmartPointer<vtkTable> dataTable)
{
	SetUpHeaders(dataTable);
	int rows = dataTable->GetNumberOfRows();
	std::vector<double> columnStatistics;
	
	for(int i = 0; i< dataTable->GetNumberOfColumns(); ++i)
	{	//compute statistics
		vtkAbstractArray *Column = dataTable->GetColumn(i);
		columnStatistics = ComputeStatistics(Column, rows);	

		
	
		for(int row = 0; row < (int)rowHeaders.size(); ++row)
		{
			this->StatisticsModel->setData(this->StatisticsModel->index(row, i), columnStatistics.at(row));
		}
	}
	
	this->StatisticsTable->setModel(StatisticsModel);
	statisticsDockWidget->setWidget(StatisticsTable);
	this->resize(500, 250);
	
}



std::vector<double> StatisticsToolbar::ComputeStatistics(vtkAbstractArray *Column, int rows)
{
	
	double average = this->Average(Column, rows);
	double stDeviation = this->StDeviation(Column, average, rows);
	QList<double> dataList = this->SortColumn(Column, rows);
    double min = dataList.at(0);
	double max = dataList.back();
	double median = floor(dataList.at(dataList.size()/2));
	double mode = Mode(dataList);

	std::cout<< average<< "\t" << stDeviation << "\t" << min << "\t" << max << "\t" << median << "\t" << mode<< endl;
	
	
	//put stats into a vector

	std::vector<double> row;
	row.push_back(average);
	row.push_back(stDeviation);
	row.push_back(min);
	row.push_back(max);
	row.push_back(median);
	row.push_back(mode);
	
	return row;
}



double StatisticsToolbar::Average(vtkAbstractArray *Column, int rows )
{
	double sum = 0;
	int j;
	for(j = 0; j< rows; j++)
	{
		sum = (Column->GetVariantValue(j).ToDouble()) + sum;
	}
	double average = (sum/rows);
	return average;
}



double StatisticsToolbar::StDeviation(vtkAbstractArray *Column, double average, int rows)
{
	double difference = 0;
	for(int j = 0; j< rows-1; ++j)
	{
		difference = (((Column->GetVariantValue(j).ToDouble() - average)*(Column->GetVariantValue(j).ToDouble() - average)) + difference); 
	}
	double stDeviation =sqrt( difference / average );
	return stDeviation;
}



QList<double> StatisticsToolbar::SortColumn(vtkAbstractArray *Column, int rows)
{
	this->myList = new QList<double>;
	this->myList->push_front(Column->GetVariantValue(0).ToDouble());
	for(int j = 1; j<rows; j++)
	{
		
		double newValue = Column->GetVariantValue(j).ToDouble();
		if(!(newValue >= 0))
			newValue = 0;
				
		if(newValue >= (this->myList->back()))//greater than last item
		{
			this->myList->push_back(newValue);

			
		}
		else if(newValue <(this->myList->front()))//less than first item
		{
			this->myList->push_front(newValue);
		}	
		else 
		{
			this->myList->push_front(newValue);//insert at beginning
			double nextValue = this->myList->at(1);
			int n = 0;//current index of new item in list
			while(newValue > nextValue) //while new number is larger than next in list..
			{
				this->myList->replace(n, nextValue);
                
				this->myList->replace((n+1), newValue);
				n = n+1;
				if(n+1 == 3)
					break;
				else
				    nextValue = this->myList->at(n+1);
			}
			//int size = this->myList->size();
		}
		
		
	}
	return *myList;

}

double StatisticsToolbar::Mode(QList<double> dataList)
{
	int maxCount = 1;
	double modeValue = -1;

	for(int i = 0; i< dataList.size(); ++i)
	{
		int count = dataList.count(dataList.at(i));
		if (count > maxCount)
		{
			maxCount = count;
		    modeValue = dataList.at(i);
		}
	}
	return modeValue;


}



void StatisticsToolbar::SetUpHeaders(vtkSmartPointer<vtkTable> dataTable)
{
	

	for (int col = 0; col< dataTable->GetNumberOfColumns(); col++)
	{
		this->colHeaders.push_back(dataTable->GetColumnName(col));
	}
	/*this->colHeaders.push_back("First");
	this->colHeaders.push_back("Second");
	this->colHeaders.push_back("Third");*/
	this->rowHeaders.push_back("Average");
	this->rowHeaders.push_back("Standard Deviation");
	this->rowHeaders.push_back("Minimum");
	this->rowHeaders.push_back("Maximum");
	this->rowHeaders.push_back("Median");
	this->rowHeaders.push_back("Mode");
	

    unsigned int numColHeaders = this->colHeaders.size();
	unsigned int numRowHeaders = this->rowHeaders.size();
	this->StatisticsModel->setColumnCount(numColHeaders);
	this->StatisticsModel->setRowCount(numRowHeaders);

	for(unsigned int i = 0; i < numColHeaders; i++)
	{
		this->StatisticsModel->setHeaderData(i,Qt::Horizontal,this->colHeaders.at(i));
		for(unsigned int j = 0; j < numRowHeaders; j++)
		{
			this->StatisticsModel->setHeaderData(j,Qt::Vertical, this->rowHeaders.at(j));
	    }
	}
	
}



	