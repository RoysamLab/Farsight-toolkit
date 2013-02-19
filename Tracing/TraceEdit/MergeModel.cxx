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
	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->Selection = new ObjectSelection();
	this->SetTraceGaps(gaps);
}

////////////////////////////////////////////////////////////////////////////////
MergeModel::~MergeModel()
{
  //explicit deletion isn't necessary because Qt objects take care of killing
  //their own children
}

////////////////////////////////////////////////////////////////////////////////
void MergeModel::SetupHeaders()
{
	std::vector<std::string> headers;
	headers.push_back("Gap ID");
	headers.push_back("Trace 1 ID");
	headers.push_back("Trace 2 ID");
	headers.push_back("Gap Distance");
	headers.push_back("Gap Angle");
	headers.push_back("Max Distance");
	headers.push_back("Path Length");
	headers.push_back("Tortuosity");
	headers.push_back("Merging Cost");
	//int numHeaders = headers.size();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i<(int)headers.size(); ++i)
	{		
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( headers.at(i).c_str() );
		this->DataTable->AddColumn(column);
	}
}
ObjectSelection * MergeModel::GetObjectSelection()
{
	return this->Selection;
}
////////////////////////////////////////////////////////////////////////////////
void MergeModel::SyncModel()
{
  if(this->GetTraceGaps().size() == 0)
    {
    return;
    }

  //clear the model 
  this->DataTable->Initialize();
  this->Selection->clear();
  this->SetupHeaders();
  //and then repopulate it with data from the trace gaps
  std::vector<TraceGap*> Gaps = this->GetTraceGaps();
  for(unsigned int i = 0;i < this->GetTraceGaps().size(); i++)
    {
		vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
		DataRow->InsertNextValue(Gaps[i]->compID);
		DataRow->InsertNextValue(Gaps[i]->Trace1->GetId());
		DataRow->InsertNextValue(Gaps[i]->Trace2->GetId());
		DataRow->InsertNextValue(Gaps[i]->dist);
		DataRow->InsertNextValue(Gaps[i]->angle);
		DataRow->InsertNextValue(Gaps[i]->maxdist);
		DataRow->InsertNextValue(Gaps[i]->length);
		DataRow->InsertNextValue(Gaps[i]->smoothness);
		DataRow->InsertNextValue(Gaps[i]->cost);
		this->DataTable->InsertNextRow(DataRow);
    }/*
  this->MapGapIDsToRows();
  this->MapTracesToRows();*/

}
vtkSmartPointer<vtkTable> MergeModel::getDataTable()
{
	return this->DataTable;
}
////////////////////////////////////////////////////////////////////////////////
void MergeModel::SelectbyTraceID(int id)
{
	unsigned int i = 0;
	int id1 =-1, id2= -1;
	bool foundT1 = false, foundT2 = false;
	while ((i < this->TraceGaps.size()) && !(foundT1 && foundT2))
	{
		if (!foundT1)
		{
			if (id == this->TraceGaps.at(i)->Trace1->GetId())
			{
				id1 = this->TraceGaps.at(i)->compID;
				foundT1 = true;
			}
		}//end search t1
		if (!foundT2)
		{
			if (id == this->TraceGaps.at(i)->Trace2->GetId())
			{
				id2 = this->TraceGaps.at(i)->compID;
				foundT2 = true;
			}
		}//end search t2
		i++;
	}
	if (foundT1 && foundT2)
	{
		if (this->Selection->isSelected((long) id1)
			!=this->Selection->isSelected((long) id2))
		{
			this->Selection->remove((long) id1);
			this->Selection->remove((long) id2);
		}
		else
		{
			this->Selection->toggle((long) id1);
			this->Selection->toggle((long) id2);
		}
	}
	else if(foundT1)
	{
		this->Selection->toggle((long) id1);
	}
	else if(foundT2)
	{
		this->Selection->toggle((long) id2);
	}
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
std::vector<int> MergeModel::GetSelectedGapIDs()
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

