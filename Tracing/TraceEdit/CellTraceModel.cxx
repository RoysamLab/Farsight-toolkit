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

#include <vector>
#include <string>
#include "TraceLine.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
CellTraceModel::CellTraceModel()
{	
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->Cells.clear();
}
CellTraceModel::CellTraceModel(std::vector<CellTrace*> Cells)
{	
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->setCells(Cells);
}
CellTraceModel::~CellTraceModel()
{	
}
void CellTraceModel::setCells(std::vector<CellTrace*> Cells)
{
	this->Cells.clear();
	this->Cells = Cells;
	this->SyncModel();
}
void CellTraceModel::SetupHeaders()
{
	this->headers.clear();
	this->headers.push_back("Root Trace");
	this->headers.push_back("Segments");
	this->headers.push_back("Stems");
	this->headers.push_back("Branch Pt");
	this->headers.push_back("Leaf Nodes");
	this->headers.push_back("Min Leaf Level");
	this->headers.push_back("Min Leaf Path Length");
	this->headers.push_back("Max Leaf Level");
	this->headers.push_back("Max Leaf Path Length");
	this->headers.push_back("Ave Leaf Level");
	this->headers.push_back("Ave Leaf Path Length");
	this->headers.push_back("Total Euclidian Length");
	this->headers.push_back("Total Path Length");
	this->headers.push_back("Average Segment Path Length");
	this->headers.push_back("Total Volume");
	this->headers.push_back("Soma X");
	this->headers.push_back("Soma Y");
	this->headers.push_back("Soma Z");
	this->headers.push_back("Trace File");
	int numHeaders = (int)this->headers.size();
	vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( this->headers.at(i).toStdString().c_str() );
		this->DataTable->AddColumn(column);
    }
}
void CellTraceModel::SyncModel()
{	
	this->DataTable->Initialize();	
	this->Selection->clear();
	this->SetupHeaders();
	for (int i = 0; i < (int) this->Cells.size(); i ++)
	{
		this->DataTable->InsertNextRow(this->Cells[i]->DataRow());
	}
}

vtkSmartPointer<vtkTable> CellTraceModel::getDataTable()
{
	return this->DataTable;
}
ObjectSelection * CellTraceModel::GetObjectSelection()
{
	return this->Selection;
}
void CellTraceModel::SelectByRootTrace(std::vector<TraceLine*> roots)
{
	this->Selection->clear();
	std::set<long int> ID;
	unsigned int i = 0;
	for (i = 0; i < roots.size(); i++)
	{
		ID.insert((long)roots.at(i)->GetId());
		//std::cout<< "root" << roots.at(i)->GetId()<< " at" <<i<< std::endl;
	}
	this->Selection->select(ID);
}
std::set<long int> CellTraceModel::GetSelecectedIDs()
{
	std::set<long int> allSelectedIDs, nextIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		int id = (int) *it;
		bool found = false;
		unsigned int j = 0;
		while(!found&&(j<this->Cells.size()))
		{
			if (id == this->Cells.at(j)->rootID())
			{
				nextIDs = this->Cells.at(j)->TraceIDsInCell();
				std::set<long>::iterator k;
				for (k = nextIDs.begin(); k != nextIDs.end(); k++)
				{
					allSelectedIDs.insert(*k);
				}
				found = true;
			}else
			{
				j++;
			}
		}//end while !found
	}//end for selected
	return allSelectedIDs;
}
std::vector<CellTrace*> CellTraceModel::GetSelecectedCells()
{
	std::vector<CellTrace*> selectedCell;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		int id = (int) *it;
		bool found = false;
		unsigned int j = 0;
		while(!found&&(j<this->Cells.size()))
		{
			if (id == this->Cells.at(j)->rootID())
			{
				selectedCell.push_back(this->Cells.at(j));
				found = true;
			}else
			{
				j++;
			}
		}//end while !found
	}//end for selected
	return selectedCell;
}
int CellTraceModel::getCellCount()
{
	return this->Cells.size();
}
CellTrace * CellTraceModel::GetCellAt( int i)
{
	CellTrace* currentCell;
	if (i < this->Cells.size())
	{
		currentCell = this->Cells.at(i);
	}
	else
	{
		currentCell = this->Cells.back();
	}
	this->Selection->select(currentCell->rootID());
	return currentCell;
}