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
#include "TraceModel.h"

TraceModel::TraceModel(std::vector<TraceLine*> trace_lines, std::vector<std::string> FeatureHeaders)
{  
	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->Selection = new ObjectSelection();
//standard headers	
	this->stdHeaders();
	if (FeatureHeaders.size() >=1)
	{
		for (int i = 0; i< (int)FeatureHeaders.size(); i++)
		{
			this->headers.push_back(FeatureHeaders[i].c_str());
		}
	}
	this->NumFeatures = this->headers.size();
	this->SetupHeaders();
	this->SetTraces(trace_lines);
}
TraceModel::TraceModel(std::vector<TraceLine*> trace_lines)
{	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
//standard headers	
	this->stdHeaders();
	this->NumFeatures = this->headers.size();
	this->SetTraces(trace_lines);
}
void TraceModel::stdHeaders()
{
	this->headers.push_back("ID");
	this->headers.push_back("# of Bits");
	this->headers.push_back("Path Length");
	this->headers.push_back("Euclidian Length");
	this->headers.push_back("Tortuosity");
	this->headers.push_back("Radius");
	this->headers.push_back("Volume");
	this->headers.push_back("Type");
	this->headers.push_back("Parent");	
	this->headers.push_back("Root ID");
	this->headers.push_back("Level");
	this->headers.push_back("Path To Root");
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
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( this->headers.at(i).toStdString().c_str() );
		this->DataTable->AddColumn(column);
    }
}
ObjectSelection * TraceModel::GetObjectSelection()
{
	return this->Selection;
}
void TraceModel::SyncModel()
{
	if(this->GetTraces().size() == 0)
	{
		return;
	}
	this->DataTable->Initialize();	
	this->Selection->clear();
	this->SetupHeaders();
	
	for (int i = 0; i < (int)this->TraceLines.size(); ++i)
	{
		vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
		
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetId());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetSize());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetLength());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetEuclidianLength());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetFragmentationSmoothness());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetRadii());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetVolume());
		DataRow->InsertNextValue((int)this->TraceLines.at(i)->GetType());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetParentID());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetRootID());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetLevel());
		DataRow->InsertNextValue(this->TraceLines.at(i)->GetPathLength());
		for (int j = 0; j < (int)this->TraceLines.at(i)->Features.size(); ++j)
		{
			DataRow->InsertNextValue(this->TraceLines.at(i)->Features.at(j));
		}
		this->DataTable->InsertNextRow(DataRow);
	}//end for traces.size  
	//this->MapTracesToRows();
}
vtkSmartPointer<vtkTable> TraceModel::getDataTable()
{
	return this->DataTable;
}
void TraceModel::SelectByIDs(int ID)
{
	this->Selection->toggle((long) ID);
	emit selectionChanged();
}
void TraceModel::SelectByIDs(std::set<long int> ID)
{
	this->Selection->add(ID);
	emit selectionChanged();
}
void TraceModel::SelectByIDs(std::vector<int> IDs)
{
	this->blockSignals(1);
	std::set<long int> ID;
	for (unsigned int i = 0; i < IDs.size(); i++)
	{
		ID.insert((long)IDs.at(i));
	}
	this->Selection->add(ID);
	this->blockSignals(0);
	emit selectionChanged();
}
std::vector<int> TraceModel::GetSelecectedIDs()
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
std::vector<TraceLine*> TraceModel::GetSelectedTraces()
{
	std::vector<TraceLine*> selectedTrace;
	std::vector<int> IDList = this->GetSelecectedIDs();
	for ( unsigned int i = 0; i< IDList.size(); i++)
	{
		bool found = false; 
		unsigned int j = 0;
		while ((!found )&&(j < this->TraceLines.size()))
		{
			if (this->TraceLines[j]->GetId()==IDList[i])
			{
				selectedTrace.push_back(this->TraceLines[j]);
				found= true;
			}
			else
			{
				j++;
			}
		}//end search for trace
	}//finished with id search
	return selectedTrace;
}
std::vector<TraceLine*> TraceModel::getRoots()
{
	std::vector<TraceLine*> roots;
	std::vector<int> IDList;
	std::vector<TraceLine*> selectedTrace = this->GetSelectedTraces();
	for (unsigned int i = 0; i < selectedTrace.size(); i++)
	{
		IDList.push_back(selectedTrace.at(i)->GetRootID());
	}
	for ( unsigned int i = 0; i< IDList.size(); i++)
	{
		bool found = false; 
		unsigned int j = 0;
		while ((!found )&&(j < this->TraceLines.size()))
		{
			if (this->TraceLines[j]->GetId()==IDList[i])
			{
				roots.push_back(this->TraceLines[j]);
				found= true;
			}
			else
			{
				j++;
			}
		}//end search for trace
	}//finished with id search
	return roots;
}
