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
#include <QtGui/QStatusBar>
#include <QtGui/QMenuBar>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QFormLayout>

#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QStringList>

#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QCheckBox>
#include <QListWidget>

#include <QtGui/QCloseEvent>
#include <QtGui/QDialog>

#include "ftkGUI/ObjectSelection.h"

#include "vtkTable.h"
#include "vtkVariant.h"
#include "ftkUtils.h"
#include "vtkSmartPointer.h"
class ClusterManager;

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
	vtkIdType GetNumberOfSelections();
	vtkIdType GetNumberOfObjects() { return NumberOfObjects; }
	std::set< vtkIdType > GetClusterIDs();
	QStringList GetClusterIDsList();

	std::set< vtkIdType > SelectionFromCluster(vtkIdType key);
	std::set< vtkIdType > GetAllSelections();

	// Object Table Functions
	bool SetObjectTable(vtkSmartPointer<vtkTable> InputObjectTable);
	vtkSmartPointer<vtkTable> GetTableOfAllSelected();
	vtkSmartPointer<vtkTable> GetTableOfSelectedFromCluster(vtkIdType key);
	std::set< vtkIdType > cluster_operator_ADD(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_SUBTRACT(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_AND(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_XOR(vtkIdType key1,vtkIdType key2);

signals:
	void SelectionChanged();
	void ClusterChanged();
	void DataChanged();
	
private:
	// Cluster map and iterator 
	std::map<vtkIdType, std::set< vtkIdType > > ClusterMap;
	std::map<vtkIdType, std::set< vtkIdType > >::iterator iter;

	//Map ID to row in table
	std::map< vtkIdType, vtkIdType> TableIDMap;
	//Table ID map Maps Object ID to Row
	std::map< vtkIdType, vtkIdType>::iterator TableIDIter;

	// Object table to cluster
	vtkSmartPointer<vtkTable> ObjectTable;
	vtkIdType NumberOfObjects;

	//Private Table manip functions
	void CopySelectedIntoTable( std::set< vtkIdType > selectedIDs, 
		vtkSmartPointer<vtkTable> selectedTable);

};

class ClusterManager : public QDialog
{
	Q_OBJECT
public:
	ClusterManager();
	void setClusteringModel(SelectiveClustering * newClusterModel);
	void setObjectSelection(ObjectSelection * ObjSelection);

public slots:

	void SelectionToClusterModification();
	void ChangeInClusters();
	void ChangeInObjectSelection();

private:

	//Selection Classes 
	SelectiveClustering * ClusterModel;
	ObjectSelection * LegacyObjectSelection;

	//Private Functions
	std::set< vtkIdType > ObjectSelectionToIDSet();

	//QT Gui Layouts
	QHBoxLayout * MainLayout;
	QListWidget * ClusterListView;

	QLabel * NumClusters;
	QLabel * NumObjects;
	QLabel * NumSelected;

	QPushButton * AddClusterButton;
};
#endif