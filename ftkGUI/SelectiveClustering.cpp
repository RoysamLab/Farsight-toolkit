/*=========================================================================

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

#include "SelectiveClustering.h"

SelectiveClustering::SelectiveClustering()
{
	/*! 
	* Initalize selective clustering
	*/
	this->ClusterMap.clear();
	this->ObjectTable = vtkSmartPointer<vtkTable>::New();
}

vtkIdType SelectiveClustering::AddCluster(std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* Create a new cluster from a selection
	*/
	vtkIdType newKey = this->ClusterMap.size();
	this->ClusterMap[newKey] = ClusterSelectionSet;
	return newKey;
}

bool SelectiveClustering::AddCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* Adds a unique new cluster
	* replaces cluster if one with the key currently exists
	*/
	this->iter = this->ClusterMap.find(key);
	if (this->iter == this->ClusterMap.end())
	{
		//Remove old item
		this->ClusterMap.erase(iter);
	}
	this->ClusterMap[key] = ClusterSelectionSet;
	return true;
}

bool SelectiveClustering::RemoveCluster(vtkIdType key)
{
	/*! 
	* Find and removes cluster and selections
	*/
	this->iter = this->ClusterMap.find(key);
	if (this->iter != this->ClusterMap.end())
	{
		//found and removed
		this->ClusterMap.erase(this->iter);
		return true;
	}
	//not found in map
	return false;
}

void SelectiveClustering::ClearClusters()
{
	/*! 
	* All Clusters must go
	*/
	this->ClusterMap.clear();
}

void SelectiveClustering::AddSelectionToCluster(vtkIdType key, vtkIdType ID)
{
	/*!
	* Modify Cluster specifed by the key to add new selections
	*/
	this->iter = this->ClusterMap.find(key);
	(*this->iter).second.insert(ID);
}

void SelectiveClustering::AddSelectionToCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* Insert range of selections into existing cluster
	*/
	this->iter = this->ClusterMap.find(key);
	(*this->iter).second.insert(ClusterSelectionSet.begin(), ClusterSelectionSet.end());
}

void SelectiveClustering::RemoveSelectionFromCluster(vtkIdType key, vtkIdType ID)
{
	/*! 
	* Modifies selection in cluster fount by key
	*/
	this->iter = this->ClusterMap.find(key);
	(*this->iter).second.erase(ID);
}

void SelectiveClustering::RemoveSelectionFromCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* finds and removes selection from specified cluster
	*/
	this->iter = this->ClusterMap.find(key);
	std::set< vtkIdType > Selection = (*this->iter).second;
	std::set< vtkIdType >::iterator LocalIterr = ClusterSelectionSet.begin();
	for (; LocalIterr != ClusterSelectionSet.end(); LocalIterr++)
	{
		Selection.erase(*LocalIterr);
	}
	(*this->iter).second = Selection;
}

std::set< vtkIdType > SelectiveClustering::SelectionFromCluster(vtkIdType key)
{
	/*! 
	* returns the selection set from specified cluster
	*/
	this->iter = this->ClusterMap.find(key);
	return (*iter).second;
}

std::set< vtkIdType > SelectiveClustering::GetAllSelections()
{
	/*! 
	* returns one selection set for all clusters
	*/
	std::set< vtkIdType > selections;
	this->iter = this->ClusterMap.begin();
	for (; this->iter != this->ClusterMap.end(); this->iter++)
	{
		std::set< vtkIdType > tempSet = (*this->iter).second;
		selections.insert(tempSet.begin(), tempSet.end());
	}
	return selections;
}

bool SelectiveClustering::SetObjectTable(vtkSmartPointer<vtkTable>InputObjectTable)
{
	/*! 
	* Set table of objects to create selections and clusters from
	*/
	vtkIdType rows = InputObjectTable->GetNumberOfRows();
	if ( rows == 0)
	{
		return false;
	}
	this->ObjectTable->Initialize();
	this->ObjectTable = InputObjectTable;
	this->NumberOfObjects = rows;
	return true;
}

vtkSmartPointer<vtkTable> SelectiveClustering::GetTableOfAllSelected()
{
	/*! 
	* Return a table containing all objects referenced in the clusters 
	*/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	selectedTable->Initialize();
	std::set< vtkIdType > selectedIDs = this->GetAllSelections();
	if ( selectedIDs.size() != 0)
	{
		return selectedTable; //should it return null?
	}

	vtkIdType NumRows = this->ObjectTable->GetNumberOfRows();
	for (vtkIdType row = 0; row <= NumRows; row++)
	{
		/*vtkIdType rowObjId = this->ObjectTable->GetValue(row, 1);
		std::set< vtkIdType >::iterator FoundAt = selectedIDs.find(rowObjId);
		if (FoundAt != selectedIDs.end())
		{
			vtkVariantArray * RowCopy = this->ObjectTable->GetRow(row);
			selectedTable->InsertNextRow(RowCopy);
		}*/
	}

	return selectedTable;
}