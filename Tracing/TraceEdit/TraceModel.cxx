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

#include <QtGui/QStandardItemModel>
#include <QtGui/QItemSelectionModel>

#include "TraceLine.h"
#include "TraceModel.h"

TraceModel::TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders)
{  
	this->Model = new QStandardItemModel(0, 0, this);
	this->SelectionModel = new QItemSelectionModel(this->Model, this);
//standard headers	
	this->headers.push_back("ID");
	this->headers.push_back("Trace Size");
	this->headers.push_back("Type");
	this->headers.push_back("Parent");
	if (FeatureHeaders.size() >=1)
	{
		for (int i = 0; i< (int)FeatureHeaders.size(); i++)
		{
			this->headers.push_back(FeatureHeaders[i].c_str());
		}
	}
	this->NumFeatures = this->headers.size();
	this->SetTraces(trace_lines);
}
TraceModel::TraceModel(std::vector<TraceLine*> trace_lines)
{	
	this->Model = new QStandardItemModel(0, 0, this);
	this->SelectionModel = new QItemSelectionModel(this->Model, this);
//standard headers	
	this->headers.push_back("ID");
	this->headers.push_back("Trace Size");
	this->headers.push_back("Type");
	this->headers.push_back("Parent");	
	this->NumFeatures = this->headers.size();
	this->SetTraces(trace_lines);
}
void TraceModel::SetTraces(std::vector<TraceLine*> trace_lines)
{
	this->TraceLines.clear();
	this->TraceLines = trace_lines;
	this->SyncModel();
}
void TraceModel::SetupHeaders()
{	
	int numHeaders = this->headers.size();
	this->Model->setColumnCount(numHeaders);
	for(int i=0; i < numHeaders; ++i)
    {
		this->Model->setHeaderData(i, Qt::Horizontal, 
			this->headers.at(i));
    }
}

void TraceModel::SyncModel()
{
	if(this->GetTraces().size() == 0)
	{
		return;
	}
	this->Model->setColumnCount(0);
	this->Model->setRowCount(0);
	this->SetupHeaders();
	std::vector< std::vector< double > > data;
	for (int i = 0; i < (int)this->TraceLines.size(); ++i)
	{
		std::vector<double> row;
		row.push_back(this->TraceLines.at(i)->GetId());
		row.push_back(this->TraceLines.at(i)->GetSize());
		row.push_back((int)this->TraceLines.at(i)->GetType());
		row.push_back(this->TraceLines.at(i)->GetParentID());
		for (int j = 0; j < (int)this->TraceLines.at(i)->Features.size(); ++j)
		{
			row.push_back(this->TraceLines.at(i)->Features.at(j));
		}
		data.push_back(row);
	}//end for traces.size  
	for (int row=0; row<(int)data.size(); ++row)
    {
    this->Model->insertRow(row);
    for(int col=0; col < this->Model->columnCount(); ++col)
      {
      this->Model->setData(this->Model->index(row, col), data.at(row).at(col));
      }
    }
	this->MapTracesToRows();
  //let the views know that the model changed
  emit modelChanged();
}
void TraceModel::MapTracesToRows()
{
	if (this->GetTraces().size() == 0)
	{
		return;
	}
	QModelIndex index;
	int ID;
	this->IDToRowMap.clear();
	for (int row = 0; row < this->Model->rowCount(); ++row)
	{
		index = this->Model->index(row, TraceModel::IDColumn);
		ID = this->Model->data(index).toInt();
		this->IDToRowMap.insert(ID, row);
	}
}

void TraceModel::SelectByIDs(int ID)
{
	int row = this->IDToRowMap.value(ID);
	QItemSelection selection;
	selection.clear();
	QModelIndex index1 = this->Model->index(row, 0, QModelIndex());
	QModelIndex index2 = this->Model->index(row, (this->Model->columnCount())-1, QModelIndex());
	selection.select(index1, index2);
	this->SelectionModel->select(selection, QItemSelectionModel::Select);
	emit selectionChanged();
}
std::vector<int> TraceModel::GetSelecectedIDs()
{
	std::vector<int> SelectedIDs;
	QModelIndexList selected = this->SelectionModel->selectedRows();
	for (int i = 0; i < selected.size(); ++i)
	{
		int row = selected.at(i).row();
		int id = this->Model->data(this->Model->index(row, TraceModel::IDColumn)).toInt();
		SelectedIDs.push_back(id);
	}
	return SelectedIDs;
}
std::vector<TraceLine*> TraceModel::GetSelectedTraces()
{
	std::vector<TraceLine*> selectedTrace;
	QModelIndexList selected = this->SelectionModel->selectedRows();
	for (int i = 0; i < selected.size(); ++i)
	{
		int row = selected.at(i).row();
		selectedTrace.push_back( this->GetTraces().at(row));
	}
	return selectedTrace;
}