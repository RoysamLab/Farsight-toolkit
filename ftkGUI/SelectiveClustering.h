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
#include <sstream>
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
#include "QvtkTableView.h"

#include "vtkTable.h"
#include "vtkVariant.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"
#include "vtkAnnotationLink.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkDataRepresentation.h"
#include "vtkCallbackCommand.h"
// PieChart
#include "vtkChartPie.h"
#include "vtkPlotPie.h"
#include "vtkPlot.h"
#include "vtkColorSeries.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include <vtkQtTableView.h>
#include <vtkQtTableModelAdapter.h>

#include "ftkUtils.h"
#include "vtkSmartPointer.h"
#include "SelectionUtilities.h"
class ClusterManager;

class SelectiveClustering: public QObject
{
	Q_OBJECT;

public:

	SelectiveClustering();
	~SelectiveClustering();

	// Add and remove clusters
	vtkIdType AddCluster(std::set<vtkIdType> ClusterSelectionSet);
	bool AddCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);
	bool RemoveCluster(vtkIdType key);
	bool RemoveCluster(vtkSmartPointer<vtkIdTypeArray> SelectedClusters);
	void ClearClusters();
	bool ValidKey(vtkIdType key);
	
	//Modify Clusters
	void AddSelectionToCluster(vtkIdType key, vtkIdType ID);
	void AddSelectionToCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);

	void RemoveSelectionFromCluster(vtkIdType key, vtkIdType ID);
	void RemoveSelectionFromCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);
	void UpdateSelectionFromCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet);

	//Cluster information 
	vtkIdType ClusterSelectionSize(vtkIdType key);
	vtkIdType NumberOfClusters();
	vtkIdType GetNumberOfSelections();
	vtkIdType GetNumberOfObjects() { return NumberOfObjects; }
	std::set< vtkIdType > GetClusterIDs();
	QStringList GetClusterIDsList();

	vtkSmartPointer<vtkTable> ClusterFeatureTable();

	std::set< vtkIdType > SelectionFromCluster(vtkIdType key);
	std::vector< long int> SelectionIDsFromCluster(vtkIdType key);

	std::set< vtkIdType > GetAllSelections();
	void CopySelectedIntoTable( std::set< vtkIdType > selectedIDs, 
		vtkSmartPointer<vtkTable> selectedTable);

	//Cluster Table

	vtkSmartPointer<vtkTable> GetClusterTable();

	// Object Table Functions
	bool SetObjectTable(vtkSmartPointer<vtkTable> InputObjectTable);
	void update();
	vtkSmartPointer<vtkTable> GetTableOfAllSelected();
	vtkSmartPointer<vtkTable> GetTableOfSelectedFromCluster(vtkIdType key);
	std::set< vtkIdType > cluster_operator_ADD(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_SUBTRACT(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_AND(vtkIdType key1,vtkIdType key2);
	std::set< vtkIdType > cluster_operator_XOR(vtkIdType key1,vtkIdType key2);
	
	vtkSmartPointer<vtkAnnotationLink> ClusterAnnotationLink;
	
	vtkSmartPointer<vtkAnnotationLink> ObjectAnnotationLink;

	std::map< vtkIdType, vtkIdType> GetObjectTableIDMap();
	std::map< vtkIdType, vtkIdType> GetClusterTableIDMap();
	void emitSelectionFinished();

signals:
	void selectionFinished();
	void SelectionChanged();
	void ClusterChanged();
	void DataChanged();
	
private:
	// Cluster map and iterator 
	std::map<vtkIdType, std::set< vtkIdType > > ClusterMap;
	std::map<vtkIdType, std::set< vtkIdType > >::iterator iter;

	//Map ID to row in table
	std::map< vtkIdType, vtkIdType> ObjectTableIDMap;
	//Table ID map Maps Object ID to Row ID
	std::map< vtkIdType, vtkIdType>::iterator TableIDIter;

	// Object table to cluster
	vtkSmartPointer<vtkTable> ObjectTable;
	vtkIdType NumberOfObjects;


	// Cluster Table 
	vtkSmartPointer<vtkTable> ClusterTable;
	//map is < OBJ ID, Table Row>
	std::map< vtkIdType, vtkIdType> ClusterTableIDMap;
	void CreateClusterTableHeaders();
	void AddRowToClusterTable(vtkIdType Key, vtkVariant ClusterSize, vtkVariant ClusterName);
	void RemoveRowFromClusterTable(vtkIdType Key);

	vtkSmartPointer<vtkVariantArray> CondenseClusterToFeatureRow(vtkIdType Key);

};

class ClusterManager : public QDialog
{
	Q_OBJECT
public:
	ClusterManager();
	~ClusterManager();
	void setManagerTitle(std::string newTitle);
	void setClusteringModel(SelectiveClustering * newClusterModel);
	void setObjectSelection(ObjectSelection * ObjSelection);
	vtkSmartPointer<vtkIdType> TableID;

public slots:

	void SelectionToClusterModification();

	void ClearClusters();
	void RemoveSelectedClusters();
	void ShowClusterFeatures();

	void ShowClusterObjectTables();
	void CloseClusterObjectTables();

	void ChangeInClusters();
	void ChangeInData();
	

	void ChangeInObjectSelection();
	void ShowDistribution();
	void RunOperatorOnSelectedClusters();

private:

	//Selection Classes 
	SelectiveClustering * ClusterModel;
	ObjectSelection * LegacyObjectSelection;


	//QT Gui Layouts
	QVBoxLayout* MainLayout;
	QHBoxLayout* HLayout ;
	QHBoxLayout* HOperatorDisplayLayout;
	QVBoxLayout* VLayout;
	QListWidget * ClusterListView;

	QvtkTableView * QVTKClusterTableView;
	QvtkTableDialog * ClusterFeatureDialog;
	std::vector<QvtkTableDialog *> ClusterObjectTables;
	
	vtkSmartPointer<vtkIdTypeArray> GetClusterTableSelections();

	QLabel * NumClusters;
	QLabel * NumObjects;
	QLabel * NumSelected;
	
	QComboBox *OperatorList, *Operand1, *Operand2, *ActionType;
	QPushButton * AddClusterButton;
	QPushButton * RunOperatorButton;
	QPushButton * ClearClusterButton;
	QPushButton * RemoveClusterButton;
	QPushButton * ClusterFeaturesButton;
	QPushButton * ShowObjectTables;
	QPushButton * HideObjectTables;
	QPushButton * ShowDistributionButton;

	bool AnnotationLinkSetUp;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
};
#endif
