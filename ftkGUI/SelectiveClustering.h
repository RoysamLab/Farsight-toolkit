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
#ifndef SELECTIVECLUSTERING_H
#define SELECTIVECLUSTERING_H

#include <set>
#include <vector>
#include <string>
#include <map>
#include <QObject>
#include "vtkTable.h"
#include "vtkVariant.h"
#include "ftkUtils.h"
#include "vtkSmartPointer.h"

class SelectiveClustering: public QObject
{
	Q_OBJECT;

public:

	SelectiveClustering();

	vtkIdType AddCluster(std::set<vtkIdType> ClusterSelectionSet);
	bool AddCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);

	bool RemoveCluster(vtkIdType key);
	void ClearClusters();
	
	void AddSelectionToCluster(vtkIdType key, vtkIdType ID);
	void AddSelectionToCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);

	void RemoveSelectionFromCluster(vtkIdType key, vtkIdType ID);
	void RemoveSelectionFromCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);

	vtkIdType ClusterSelectionSize(vtkIdType key);
	vtkIdType NumberOfClusters();

	std::set< vtkIdType > SelectionFromCluster(vtkIdType key);
	std::set< vtkIdType > GetAllSelections();

	// Object Table Functions
	bool SetObjectTable(vtkSmartPointer<vtkTable> InputObjectTable);
	vtkSmartPointer<vtkTable> GetTableOfAllSelected();
	vtkSmartPointer<vtkTable> GetTableOfSelectedFromCluster(vtkIdType key);

signals:
	void SelectionChanged();
	void ClusterChanged();
	void DataChanged();
	
private:
	// Cluster map and iterator 
	std::map<vtkIdType, std::set< vtkIdType > > ClusterMap;
	std::map<vtkIdType, std::set< vtkIdType > >::iterator iter;

	// Object table to cluster
	vtkSmartPointer<vtkTable> ObjectTable;
	vtkIdType NumberOfObjects;
};
#endif