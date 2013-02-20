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
	this->dataTable = vtkSmartPointer<vtkTable>::New();
	this->selectionTable = vtkSmartPointer<vtkTable>::New();
	this->StatisticsTable = new QTableView();
	this->StatisticsModel = new QStandardItemModel(this);
	this->StatisticsSelectionModel = new QItemSelectionModel(this->StatisticsModel, this);	
	this->statisticsDockWidget = new QDockWidget(tr("Statistics Toolbar"));
	this->Selection = new ObjectSelection();
	std::vector<int> selectedRowNumbers;
	
	

}


StatisticsToolbar::~StatisticsToolbar()
{
}



void StatisticsToolbar::setTable(vtkSmartPointer<vtkTable> inputDataTable, ObjectSelection * Selection)
{ 
	dataTable = NULL;
	SetUpHeaders(inputDataTable);
	this->Selection = Selection;
	std::vector<double> selectedData;
	std::vector<int> IDList = this->GetSelectedIDs();
	this->selectedRowNumbers.clear();
	//bool found;	
    dataTable = inputDataTable;

	if(IDList.size() > 0)
	{
		
		for(unsigned int i = 0; i < IDList.size(); i++)
		{
			bool found = false;
			unsigned int j = 0;
			while ((!found) && (j< inputDataTable->GetNumberOfRows()) )
			{
				if(IDList[i] == inputDataTable->GetValue(j,0).ToInt())
				{
				this->selectedRowNumbers.push_back(j);//inserted row number corresponding to selected row of input table
				found = true;
				}
				else
				{
					j++;
				}
			}
		}
	}
	
	
	

	int rows = dataTable->GetNumberOfRows();
	std::cout << rows << endl;
	std::vector<double> columnStatistics;
	
	for(int i = 1; i< dataTable->GetNumberOfColumns(); ++i)
	{	//compute statistics
		vtkAbstractArray *Column = dataTable->GetColumn(i);
		columnStatistics = ComputeStatistics(Column, rows);	
		
		for(int row = 0; row < (int)rowHeaders.size(); ++row)
		{
			this->StatisticsModel->setData(this->StatisticsModel->index(row, i-1), columnStatistics.at(row));
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
	unsigned int j=0;
	double average;
	
	if(this->selectedRowNumbers.size() > 0)
	{
		for(unsigned int i=0; i< this->selectedRowNumbers.size(); i++)
		{
			j = this->selectedRowNumbers.at(i);//set row #
			sum = (Column->GetVariantValue(j).ToDouble()) + sum;
		}
		average = (sum/(selectedRowNumbers.size()));
	}

	else
	{
	

	for(j = 0; (int)j< rows; j++)
	{
		sum = (Column->GetVariantValue(j).ToDouble()) + sum;
	}
	average = (sum/rows);
	}

	return average;
}



double StatisticsToolbar::StDeviation(vtkAbstractArray *Column, double average, int rows)
{
	double difference = 0;
    double stDeviation = 0;
	unsigned int j = 0;

	if(this->selectedRowNumbers.size() > 0)
	{
		for(unsigned int i=0; i< this->selectedRowNumbers.size(); i++)
		{
			j = this->selectedRowNumbers.at(i);//set row #
			difference = (((Column->GetVariantValue(j).ToDouble() -average) * (Column->GetVariantValue(j).ToDouble() - average)) + difference);
			
		}
		stDeviation = sqrt( difference / (selectedRowNumbers.size()-1) );
	}
	else
	{
		for(j = 0; (int)j< rows; ++j)
		{
			difference = (((Column->GetVariantValue(j).ToDouble() - average)*(Column->GetVariantValue(j).ToDouble() - average)) + difference); 
			
		}
		stDeviation = sqrt( difference / (rows-1) );
	}
	
	return stDeviation;
}


QList<double> StatisticsToolbar::SortColumn(vtkAbstractArray *Column, int rows)
{
	this->myList = new QList<double>;
	double newValue;
	unsigned int j = 0;
	
	if(this->selectedRowNumbers.size() > 0)
	{
		j = this->selectedRowNumbers.at(0);//set row #
		newValue = Column->GetVariantValue(j).ToDouble();
		this->myList->push_back(newValue);

		for(unsigned int i=1; i< this->selectedRowNumbers.size(); i++)
		{
			j = this->selectedRowNumbers.at(i);//set row #
			newValue = Column->GetVariantValue(j).ToDouble();
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
			else //in middle?
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

	}

	else//nothing selected, compute on whole table
	{
		this->myList->push_front(Column->GetVariantValue(0).ToDouble());
		for(j = 1; (int)j<rows; j++)
		{
		
			newValue = Column->GetVariantValue(j).ToDouble();
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
	

	for (int col = 1; col< dataTable->GetNumberOfColumns(); col++)//start at 1 so id column is not included
	{
		this->colHeaders.push_back(dataTable->GetColumnName(col));
	}
	
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


std::vector<int> StatisticsToolbar::GetSelectedIDs()
{   
	std::vector<int> SelectedIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{	
		SelectedIDs.push_back(*it);
	}
	
	return SelectedIDs;
	
	
}

//void StatisticsToolbar::makeTable(vtkSmartPointer<vtkTable> inputDataTable)
//{
//	std::vector<double> selectedData;
//	std::vector<int> IDList = this->GetSelectedIDs();
//	this->selectedRowNumbers.clear();
//	bool found;	
//
//	if(IDList.size() > 0)
//	{
//		
//		for(unsigned int i = 0; i < IDList.size(); i++)
//		{
//			bool found = false;
//			unsigned int j = 0;
//			while ((!found) && (j< inputDataTable->GetNumberOfRows()) )
//			{
//				if(IDList[i] == inputDataTable->GetValue(j,0).ToInt())
//				{
//				this->selectedRowNumbers.push_back(j);//inserted row number corresponding to selected row of input table
//				found = true;
//				}
//				else
//				{
//					j++;
//				}
//			}
//		}
//	}
//	
//	this->dataTable = inputDataTable;	
//	
//}










	
