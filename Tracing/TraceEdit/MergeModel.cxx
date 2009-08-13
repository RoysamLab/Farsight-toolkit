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

#include <QtGui/QStandardItemModel>
#include <QtGui/QItemSelectionModel>

#include "TraceGap.h"
#include "TraceLine.h"
#include "MergeModel.h"

////////////////////////////////////////////////////////////////////////////////
MergeModel::MergeModel(std::vector<TraceGap*> gaps)
{
  this->Model = new QStandardItemModel(0,0);
  this->SelectionModel = new QItemSelectionModel(this->Model);
  this->SetTraceGaps(gaps);
}

////////////////////////////////////////////////////////////////////////////////
MergeModel::~MergeModel()
{
  delete this->Model;
  this->Model = NULL;
  delete this->SelectionModel;
  this->SelectionModel = NULL;
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SetupHeaders()
{
  std::vector<QString> headers;
  headers.push_back("ID");
  headers.push_back("Trace 1");
  headers.push_back("Trace 2");
  headers.push_back("Gap");
  headers.push_back("Angle");
  headers.push_back("D");
  headers.push_back("Length");
  headers.push_back("Smoothness");
  headers.push_back("Cost");
  int numHeaders = headers.size();
  this->Model->setColumnCount(numHeaders);
  for(int i=0; i<(int)headers.size(); ++i)
    {
    this->Model->setHeaderData(i, Qt::Horizontal, headers.at(i));
    }
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SyncModel()
{
  if(this->GetTraceGaps().size() == 0)
    {
    return;
    }

  //clear the model
  this->Model->setColumnCount(0);
  this->Model->setRowCount(0);
  this->SetupHeaders();
  //and then repopulate it with data from the trace gaps
  std::vector< std::vector< double > > data;
  std::vector<TraceGap*> Gaps = this->GetTraceGaps();
  for(unsigned int i = 0;i < this->GetTraceGaps().size(); i++)
    {
    std::vector<double> row;
	row.push_back(Gaps[i]->compID);
    row.push_back(Gaps[i]->Trace1->GetId());
    row.push_back(Gaps[i]->Trace2->GetId());
    row.push_back(Gaps[i]->dist);
    row.push_back(Gaps[i]->angle);
    row.push_back(Gaps[i]->maxdist);
    row.push_back(Gaps[i]->length);
    row.push_back(Gaps[i]->smoothness);
    row.push_back(Gaps[i]->cost);
    data.push_back(row);
    }

  for (int row=0; row<(int)this->GetTraceGaps().size(); ++row)
    {
    //create a new row
    this->Model->insertRow(row);
    //insert the data for a gap in this row
    for(int col=0; col < this->Model->columnCount(); ++col)
      {
      this->Model->setData(this->Model->index(row, col), data.at(row).at(col));
      }
    }
  //let the views know that the model changed
  emit modelChanged();

  //and refresh the mapping between map ids and rows
  this->MapGapIDsToRows();
  this->MapTracesToRows();

}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::MapGapIDsToRows()
{
  if(this->GetTraceGaps().size() == 0)
    {
		return;
    }
  QModelIndex index;
  int ID;
  this->IDToRowMap.clear();
  for (int row = 0; row < this->Model->rowCount(); row++)
    {
    index = this->Model->index(row, MergeModel::IDColumn);
    ID = this->Model->data(index).toInt();
    this->IDToRowMap.insert(ID,row);
    }
}
void MergeModel::MapTracesToRows()
{
	if(this->GetTraceGaps().size() == 0)
    {
		return;
    }
	QModelIndex index1, index2;
	int T1, T2;
	this->Trace1ToRow.clear();
	this->Trace2ToRow.clear();
	for (int row = 0; row < this->Model->rowCount(); row++)
    {
//map trace 1
		index1 = this->Model->index(row, this->Trace1Col);
		T1 = this->Model->data(index1).toInt();
		this->Trace1ToRow.insert(T1,row);
//map trace 2		
		index2 = this->Model->index(row, this->Trace2Col);
		T2 = this->Model->data(index2).toInt();
		this->Trace2ToRow.insert(T2,row);
    }
}
void MergeModel::SelectbyTraceID(int id)
{
	int row1 = -1 , row2= -1;
	row1 = this->Trace1ToRow.value(id);
	row2 = this->Trace2ToRow.value(id);
	QItemSelection selection;
	selection.clear();
	if (row1 > 0)
	{
		QModelIndex index1 = this->Model->index(row1, 0, QModelIndex());
		QModelIndex index2 = this->Model->index(row1,(this->Model->columnCount())-1, QModelIndex());
		selection.select(index1,index2);
		this->SelectionModel->select(selection, QItemSelectionModel::Toggle);
		std::cout << "\tRow 1\t" << row1;
	}
	if (row2 > 0)
	{	
		QModelIndex index1 = this->Model->index(row2, 0, QModelIndex());
		QModelIndex index2 = this->Model->index(row2,(this->Model->columnCount())-1, QModelIndex());
		selection.select(index1,index2);
		this->SelectionModel->select(selection, QItemSelectionModel::Toggle);
		std::cout << "\tRow 2\t" << row2 << std::endl; 
	}	
	//emit modelChanged();	
}
////////////////////////////////////////////////////////////////////////////////
void MergeModel::SetTraceGaps(std::vector<TraceGap *> gaps)
{
  this->TraceGaps = gaps;
  this->SyncModel();
}

////////////////////////////////////////////////////////////////////////////////
std::vector<TraceGap *> MergeModel::GetTraceGaps()
{
  return this->TraceGaps;
}

////////////////////////////////////////////////////////////////////////////////
QStandardItemModel* MergeModel::GetModel()
{
  return this->Model;
}

////////////////////////////////////////////////////////////////////////////////
QItemSelectionModel* MergeModel::GetSelectionModel()
{
  return this->SelectionModel;
}

/**
 * Call this slot when deleting gaps.
 * Note that actual changes to the trace data are performed by the TraceObject
 * class.
 **/
////////////////////////////////////////////////////////////////////////////////
void MergeModel::deleteTrigger(void)
{
  //Extract a list of IDs
  QModelIndexList selectedIndices = this->SelectionModel->selectedRows();
  std::vector<int> ids(0);
  for (int selectedIndex = 0;
       selectedIndex < selectedIndices.size();
       ++selectedIndex)
    {
    int row = selectedIndices.at(selectedIndex).row();
    QList<QStandardItem *> items = this->Model->takeRow(row);
    for(int i = 0; i < items.size(); ++i)
      {
      delete items.at(i);
      }
    }
  this->MapGapIDsToRows();
  emit modelChanged();
}

/**
 * Call this slot when merging gaps.
 * Note that actual changes to the trace data are performed by the TraceObject
 * class.
 **/

 //how is this different?  the gap still goes bye-bye.
 //I guess it needs to emit a different signal so View3D knows what to do.
 //View3D already knows what to do, it's the one driving the action.
 //but it does need to know when the model's selection changes.
////////////////////////////////////////////////////////////////////////////////
void MergeModel::mergeTrigger(void)
{
}

////////////////////////////////////////////////////////////////////////////////
std::vector<int> MergeModel::GetSelectedGapIDs()
{
  std::vector<int> selectedGapIDs;
	QModelIndexList selectedIndices = this->SelectionModel->selectedRows();
	for (int selectedIndex = 0;
       selectedIndex < selectedIndices.size();
       ++selectedIndex) 
	  {
		int row = selectedIndices.at(selectedIndex).row();
		int id = this->Model->data(
      this->Model->index(row, MergeModel::IDColumn) ).toInt();
		selectedGapIDs.push_back(id);
	  }
  return selectedGapIDs;
}

