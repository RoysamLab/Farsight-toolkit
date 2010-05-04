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
	this->headers.push_back("Leaf Nodes");
	this->headers.push_back("Min Terminal Level");
	this->headers.push_back("Max Terminal Level");
	this->headers.push_back("Total Euclidian Length");
	this->headers.push_back("Total Path Length");
	this->headers.push_back("Total Volume");
	this->headers.push_back("ave Terminal Path Length");
	int numHeaders = (int)this->headers.size();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkDoubleArray>::New();
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