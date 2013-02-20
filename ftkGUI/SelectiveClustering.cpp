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

	this->ObjectTableIDMap.clear();
	this->ObjectTable = vtkSmartPointer<vtkTable>::New();

	this->ClusterTable = vtkSmartPointer<vtkTable>::New();
	this->ClusterTableIDMap.clear();
	this->CreateClusterTableHeaders();

	this->ClusterAnnotationLink = vtkSmartPointer<vtkAnnotationLink>::New();

	this->ObjectAnnotationLink = vtkSmartPointer<vtkAnnotationLink>::New();
}

SelectiveClustering::~SelectiveClustering()
{
}

vtkIdType SelectiveClustering::AddCluster(std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* Create a new cluster from a selection
	*/
	vtkIdType newKey = this->ClusterMap.size();
	this->iter = this->ClusterMap.find(newKey);
	while ((this->ClusterMap.size() > 0) &&(this->iter != this->ClusterMap.end()))
	{
		newKey++;
		this->iter = this->ClusterMap.find(newKey);
	}
	this->ClusterMap[newKey] = ClusterSelectionSet;
	this->AddRowToClusterTable(newKey, vtkVariant(ClusterSelectionSet.size()), NULL);

	emit ClusterChanged();
	return newKey;
}

void SelectiveClustering::emitSelectionFinished()
{
	emit selectionFinished();
}

bool SelectiveClustering::AddCluster(vtkIdType key, std::set<vtkIdType> ClusterSelectionSet)
{
	/*! 
	* Adds a unique new cluster
	* replaces cluster if one with the key currently exists
	*/
	this->iter = this->ClusterMap.find(key);
	if ((this->ClusterMap.size() > 0) && (this->iter != this->ClusterMap.end()))
	{
		//Remove old item
		this->ClusterMap.erase(iter);
	}
	this->ClusterMap[key] = ClusterSelectionSet;
	this->AddRowToClusterTable(key, vtkVariant(ClusterSelectionSet.size()), NULL);
	emit ClusterChanged();
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
		SelectionUtilities::RemoveRowAndReMapTable(key, this->ClusterTable, this->ClusterTableIDMap);
		//this->RemoveRowFromClusterTable(key);
		emit ClusterChanged();
		return true;
	}
	//not found in map
	return false;
}

bool SelectiveClustering::RemoveCluster(vtkSmartPointer<vtkIdTypeArray> SelectedClusters)
{
	/*! 
	* Find and removes multiple clusters and selections
	*/
	bool removed = false;
	for (vtkIdType count = 0; count < SelectedClusters->GetSize(); count++)
	{
		vtkIdType key = SelectedClusters->GetValue(count);
		//std::cout<< "\n" << key;
		this->iter = this->ClusterMap.find(key);
		if (this->iter != this->ClusterMap.end())
		{
			//found and removed
			//std::cout<< " removed ";
			this->ClusterMap.erase(this->iter);
			SelectionUtilities::RemoveRowAndReMapTable(key, this->ClusterTable, this->ClusterTableIDMap);
			//this->RemoveRowFromClusterTable(key);
			removed = true;
		}
		//
	}
	if (removed)
	{
		emit ClusterChanged();
	}
	return removed;
}

void SelectiveClustering::ClearClusters()
{
	/*! 
	* All Clusters must go
	*/
	this->ClusterMap.clear();
	this->CreateClusterTableHeaders();
	this->ClusterTableIDMap.clear();
	emit ClusterChanged();
}


bool SelectiveClustering::ValidKey(vtkIdType key)
{
	/*!
	* 
	*/
	return (this->ClusterMap.end() != this->ClusterMap.find(key));
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

vtkIdType SelectiveClustering::ClusterSelectionSize(vtkIdType key)
{
	/*! 
	* How Many objects in the cluster
	*/
	this->iter = this->ClusterMap.find(key);
	std::set< vtkIdType > Selection = (*this->iter).second;
	return (vtkIdType) Selection.size();
}

vtkIdType SelectiveClustering::NumberOfClusters()
{
	/*! 
	* Count of clusters created
	*/
	return (vtkIdType) this->ClusterMap.size();
}

vtkIdType SelectiveClustering::GetNumberOfSelections()
{
	/*! 
	* Count of total number of selections
	*/
	vtkIdType totalCount = 0;
	this->iter = this->ClusterMap.begin();
	for (; this->iter != this->ClusterMap.end(); this->iter++)
	{
		std::set< vtkIdType > tempSet = (*this->iter).second;
		totalCount += tempSet.size();
	}
	return totalCount;
}

std::set< vtkIdType > SelectiveClustering::GetClusterIDs()
{
	/*! 

	* returns set containing the ID of each cluster
	*/
	std::set< vtkIdType > ClusterIDs;
	this->iter = this->ClusterMap.begin();
	for (; this->iter != this->ClusterMap.end(); this->iter++)
	{
		ClusterIDs.insert((*this->iter).first);
	}
	return ClusterIDs;
}
QStringList SelectiveClustering::GetClusterIDsList()
{
	/*! 
	* Returns QSTring List of Cluster Ids for display
	*/
	QStringList clusterList;
	this->iter = this->ClusterMap.begin();
	for (; this->iter != this->ClusterMap.end(); this->iter++)
	{
		clusterList << QString::number((int)(*this->iter).first);
	}
	return clusterList;
}
std::set< vtkIdType > SelectiveClustering::SelectionFromCluster(vtkIdType key)
{
	/*! 
	* returns the selection set from specified cluster
	*/
	this->iter = this->ClusterMap.find(key);
	return (*iter).second;
}

std::vector< long int> SelectiveClustering::SelectionIDsFromCluster(vtkIdType key)
{
	this->iter = this->ClusterMap.find(key);
	std::set< vtkIdType > clusterSet = (*iter).second;
	std::vector< long int> ids;
	std::set< vtkIdType >::iterator iter;
	for( iter = clusterSet.begin(); iter != clusterSet.end(); iter++)
	{
		if(ObjectTableIDMap.find( *iter) != ObjectTableIDMap.end())
		{
			ids.push_back((long int)ObjectTableIDMap.find( *iter)->second);
		}
	}
	return ids;
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

vtkSmartPointer<vtkTable> SelectiveClustering::GetClusterTable()
{
	return this->ClusterTable;
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
	this->ObjectTable->DeepCopy( InputObjectTable);
	this->update();
	return true;
}

void SelectiveClustering::update()
{
	this->NumberOfObjects = this->ObjectTable->GetNumberOfRows();
	for ( vtkIdType currow = 0; currow < this->NumberOfObjects; currow++ )
	{
		vtkIdType rowObjId = this->ObjectTable->GetValue( currow, 0 ).ToTypeInt64();
		this->ObjectTableIDMap[ rowObjId ] = currow;
	}
	this->ClusterMap.clear();
	this->CreateClusterTableHeaders();
	this->ClusterTableIDMap.clear();
	emit DataChanged();
}

vtkSmartPointer<vtkTable> SelectiveClustering::GetTableOfAllSelected()
{
	/*! 
	* Return a table containing all objects referenced in the clusters 
	*/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	selectedTable->Initialize();
	std::set< vtkIdType > selectedIDs = this->GetAllSelections();
	//std::cout<< " Total obj Sel " << selectedIDs.size() << std::endl;
	if ( selectedIDs.size() == 0)
	{
		return selectedTable; //should it return null?
	}
	this->CopySelectedIntoTable(selectedIDs, selectedTable);

	return selectedTable;
}

vtkSmartPointer<vtkTable> SelectiveClustering::GetTableOfSelectedFromCluster(vtkIdType key)
{
	/*! 
	* Return a table containing all objects from specified cluster
	*/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	selectedTable->Initialize();
	std::set< vtkIdType > selectedIDs = this->SelectionFromCluster(key);
	if ( selectedIDs.size() == 0)
	{
		return selectedTable; //should it return null?
	}
	this->CopySelectedIntoTable(selectedIDs, selectedTable);

	return selectedTable;
}

void SelectiveClustering::CopySelectedIntoTable(std::set<vtkIdType> selectedIDs, vtkSmartPointer<vtkTable> selectedTable)
{
	/*! 
	* Populate a table containing selected objects 
	*/
	//Add Header Columns
	for(vtkIdType NumberOfColumns = 0; NumberOfColumns < this->ObjectTable->GetNumberOfColumns(); NumberOfColumns++ )
	{
		vtkSmartPointer<vtkVariantArray> col = vtkSmartPointer<vtkVariantArray>::New();
		col->SetName(this->ObjectTable->GetColumnName(NumberOfColumns));
		selectedTable->AddColumn(col);
	}
	// 
	std::set< vtkIdType >::iterator selectionIter = selectedIDs.begin();
	for (; selectionIter != selectedIDs.end(); selectionIter++)
	{
		vtkIdType curSelection = *selectionIter;
		this->TableIDIter = this->ObjectTableIDMap.find(curSelection);
		if (this->TableIDIter != this->ObjectTableIDMap.end())
		{
			vtkVariantArray * RowCopy = this->ObjectTable->GetRow((*this->TableIDIter).second);
			selectedTable->InsertNextRow(RowCopy);
		}
	}
	////Code to compare timing map vs search

	//vtkIdType NumRows = this->ObjectTable->GetNumberOfRows();
	//for (vtkIdType row = 0; row <= this->NumberOfObjects; row++)
	//{
	//	vtkIdType rowObjId = this->ObjectTable->GetValue(row, 0).ToTypeInt64();
	//	//std::cout << "Searching for obj: " << rowObjId << std::endl;
	//	std::set< vtkIdType >::iterator FoundAt = selectedIDs.find(rowObjId);
	//	if (FoundAt != selectedIDs.end())
	//	{
	//		//std::cout << "found obj: " << rowObjId << std::endl;
	//		vtkVariantArray * RowCopy = this->ObjectTable->GetRow(row);
	//		selectedTable->InsertNextRow(RowCopy);
	//	}
	//}
}

void SelectiveClustering::CreateClusterTableHeaders()
{
	/*!
	* Creates an empty table then sets col headers to id size and name
	*/
	this->ClusterTable->Initialize();

	vtkSmartPointer<vtkVariantArray> IDcol = vtkSmartPointer<vtkVariantArray>::New();
	IDcol->SetName("Cluster ID");
	this->ClusterTable->AddColumn(IDcol);

	vtkSmartPointer<vtkVariantArray> Countcol = vtkSmartPointer<vtkVariantArray>::New();
	Countcol->SetName("Number Of Objects");
	this->ClusterTable->AddColumn(Countcol);

	vtkSmartPointer<vtkVariantArray> Namecol = vtkSmartPointer<vtkVariantArray>::New();
	Namecol->SetName("Cluster Name");
	this->ClusterTable->AddColumn(Namecol);
}

void SelectiveClustering::AddRowToClusterTable(vtkIdType Key, vtkVariant ClusterSize, vtkVariant ClusterName)
{
	/*!
	* Builds table of Cluster ID size and name if any
	*/
	vtkSmartPointer<vtkVariantArray> NewRow = vtkSmartPointer<vtkVariantArray>::New();
	NewRow->InsertNextValue(Key); 
	NewRow->InsertNextValue(ClusterSize); 
	NewRow->InsertNextValue(ClusterName);
	vtkIdType RowAdded = this->ClusterTable->InsertNextRow(NewRow);
	this->ClusterTableIDMap[Key] = RowAdded;
}

void SelectiveClustering::RemoveRowFromClusterTable(vtkIdType Key)
{
	/*!
	* removes entry from cluster table 
	*/
	this->TableIDIter = this->ClusterTableIDMap.find(Key);
	if (this->TableIDIter != this->ClusterTableIDMap.end())
	{
		this->ClusterTable->RemoveRow((*this->TableIDIter).second);
	}
	for (vtkIdType row = 0; row != this->ClusterTable->GetNumberOfRows(); row++)
	{
		vtkIdType value = this->ClusterTable->GetValue(row,0).ToTypeInt64();
		this->ClusterTableIDMap[value] = row;
	}
}

vtkSmartPointer<vtkVariantArray> SelectiveClustering::CondenseClusterToFeatureRow(vtkIdType Key)
{
	/*!
	* Returns a row from the statistics of the
	* features of the objects in each cluster 
	*/
	vtkSmartPointer<vtkTable> clustersTable = this->GetTableOfSelectedFromCluster(Key);

	vtkSmartPointer<vtkVariantArray> rowForCluster = vtkSmartPointer<vtkVariantArray>::New();
	rowForCluster->Initialize();
	rowForCluster->InsertNextValue(Key);

	for(vtkIdType col= 1; col != clustersTable->GetNumberOfColumns(); col ++)
	{
		vtkIdType rowCount = clustersTable->GetNumberOfRows();
		//
		double sum = 0; 
		vtkAbstractArray * columnData = clustersTable->GetColumn(col);
		for (vtkIdType row = 0; row != rowCount; row++)
		{
			sum += columnData->GetVariantValue(row).ToDouble();
		}
		sum = sum / (double)rowCount;
		rowForCluster->InsertNextValue(sum);
	}
	return rowForCluster;
}

vtkSmartPointer<vtkTable> SelectiveClustering::ClusterFeatureTable()
{
	/*!
	* 
	*/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	selectedTable->Initialize();	
	
	for(vtkIdType NumberOfColumns = 0; NumberOfColumns < this->ObjectTable->GetNumberOfColumns(); NumberOfColumns++ )
	{
		vtkSmartPointer<vtkVariantArray> col = vtkSmartPointer<vtkVariantArray>::New();
		col->SetName(this->ObjectTable->GetColumnName(NumberOfColumns));
		selectedTable->AddColumn(col);
	}

	this->iter = this->ClusterMap.begin();
	for (; this->iter != this->ClusterMap.end(); this->iter++)
	{
		vtkSmartPointer<vtkVariantArray> NextCluster = this->CondenseClusterToFeatureRow((*this->iter).first);
		selectedTable->InsertNextRow(NextCluster);
	}
	//selectedTable->Dump(16);
	return selectedTable;

}

std::set< vtkIdType > SelectiveClustering::cluster_operator_ADD(vtkIdType key1, vtkIdType key2)
{
	/*!
	*  	adds two clusters - Union
	*   @param1 : vtkIdType key1 argument, 
	*	@param2	: vtkIdType argument
	*	@return the vtktable with the result
	**/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	std::set< vtkIdType > selectedIDs = this->SelectionFromCluster(key1);
	std::set< vtkIdType > tempSet = this->SelectionFromCluster(key2);
	selectedIDs.insert(tempSet.begin(), tempSet.end());
	return selectedIDs;
}


std::set< vtkIdType > SelectiveClustering::cluster_operator_SUBTRACT(vtkIdType key1, vtkIdType key2)
{
	/*!
	*  subtracts cluster with key2 from cluster with key1
	*   @param1 : vtkIdType key1 argument, 
	*	@param2	: vtkIdType argument
	*	@return the vtktable with the result
	**/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	std::set< vtkIdType > selectedIDs = this->SelectionFromCluster(key1);
	std::set< vtkIdType > tempSet = this->SelectionFromCluster(key2);
	std::set< vtkIdType > result ;
	std::set_difference(selectedIDs.begin(),selectedIDs.end(),tempSet.begin(),tempSet.end(),std::inserter(result, result.end()));
	return result;
}

std::set< vtkIdType > SelectiveClustering::cluster_operator_AND(vtkIdType key1, vtkIdType key2)
{
	/*!
	*  	subtracts two clusters
	*   @param1 : vtkIdType key1 argument
	*	@param2	: vtkIdType argument
	*	@return the vtktable with the result
	**/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	std::set< vtkIdType > selectedIDs = this->SelectionFromCluster(key1);
	std::set< vtkIdType > tempSet = this->SelectionFromCluster(key2);
	std::set< vtkIdType > result ;
	std::set_intersection(selectedIDs.begin(),selectedIDs.end(),tempSet.begin(),tempSet.end(),std::inserter(result, result.end()));
	return result;
}

std::set< vtkIdType > SelectiveClustering::cluster_operator_XOR(vtkIdType key1, vtkIdType key2)
{
	/*!
	*  	subtracts two clusters
	*   @param1 : vtkIdType key1 argument
	*	@param2	: vtkIdType argument
	*	@return the vtktable with the result
	**/
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	std::set< vtkIdType > selectedIDs = this->SelectionFromCluster(key1);
	std::set< vtkIdType > tempSet = this->SelectionFromCluster(key2);
	std::set< vtkIdType > result ;
	std::set_symmetric_difference(selectedIDs.begin(),selectedIDs.end(),tempSet.begin(),tempSet.end(),std::inserter(result, result.end()));
	return result;
}

std::map< vtkIdType, vtkIdType> SelectiveClustering::GetObjectTableIDMap()
{
	/*! 
	* Find Object Row by ID
	*/
	return this->ObjectTableIDMap;
}


std::map< vtkIdType, vtkIdType> SelectiveClustering::GetClusterTableIDMap()
{
	/*! 
	* Find Cluster Row by Cluster ID
	*/
	return this->ClusterTableIDMap;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

ClusterManager::ClusterManager()
{
	/*! 
	* Create Gui 
	*/
	this->ClusterModel = new SelectiveClustering();
	this->setModal(false);

	this->NumObjects = new QLabel(" None ");
	this->NumClusters = new QLabel(" None ");
	this->NumSelected = new QLabel(" None ");

	this->AddClusterButton = new QPushButton(" Add Cluster ");
	this->RunOperatorButton = new QPushButton(" Run Operator ");

	connect(this->AddClusterButton, SIGNAL(clicked()), this, SLOT(SelectionToClusterModification()));
	connect(this->RunOperatorButton, SIGNAL(clicked()), this, SLOT(RunOperatorOnSelectedClusters()));

	this->ClearClusterButton = new QPushButton(" Clear Cluster ");
	connect(this->ClearClusterButton, SIGNAL(clicked()), this, SLOT(ClearClusters()));

	this->RemoveClusterButton = new QPushButton(" Remove Cluster ");
	connect(this->RemoveClusterButton, SIGNAL(clicked()), this, SLOT(RemoveSelectedClusters()));

	this->ClusterFeaturesButton = new QPushButton(" Show Cluster ");
	connect(this->ClusterFeaturesButton, SIGNAL(clicked()), this, SLOT(ShowClusterFeatures()));

	this->ShowObjectTables = new QPushButton(" Show Cluster Tables");
	connect(this->ShowObjectTables, SIGNAL(clicked()), this, SLOT(ShowClusterObjectTables()));

	this->HideObjectTables = new QPushButton(" Close Cluster Tables");
	connect(this->HideObjectTables, SIGNAL(clicked()), this, SLOT(CloseClusterObjectTables()));
	
	this->ShowDistributionButton = new QPushButton(" Show Distribution ");
	connect(this->ShowDistributionButton, SIGNAL(clicked()), this, SLOT(ShowDistribution()));
	
	//this->ClusterListView = new QListWidget(this);

	this->QVTKClusterTableView = new QvtkTableView();
	this->ClusterFeatureDialog = NULL;
	this->ClusterObjectTables.clear();
	
	this->OperatorList = new QComboBox(this);
	OperatorList->addItem("ADD");
	OperatorList->addItem("SUBRACT");
	OperatorList->addItem("AND");
	OperatorList->addItem("XOR");
	
	this->ActionType = new QComboBox(this);	
	ActionType->addItem("Display");
	ActionType->addItem("Override I");
	ActionType->addItem("Override II");


	this->Operand1 = new QComboBox(this);
	this->Operand2  = new QComboBox(this);
	
	QFormLayout *InfoLayout = new QFormLayout();
	InfoLayout->addRow("Number of Objects: ", this->NumObjects);
	InfoLayout->addRow("Number of Clusters: ", this->NumClusters);
	InfoLayout->addRow("Number Selected: ", this->NumSelected);
	InfoLayout->addRow(this->ClearClusterButton);
	InfoLayout->addRow(this->AddClusterButton);
	InfoLayout->addRow(this->RemoveClusterButton);

	HOperatorDisplayLayout= new QHBoxLayout();
	this->HOperatorDisplayLayout->addWidget(Operand1);
	this->HOperatorDisplayLayout->addWidget(OperatorList);
	this->HOperatorDisplayLayout->addWidget(ActionType);
	this->HOperatorDisplayLayout->addWidget(Operand2);
	this->HOperatorDisplayLayout->addWidget(this->RunOperatorButton);
	//InfoLayout->addRow(this->HOperatorDisplayLayout);

	QHBoxLayout* ButtonLayout = new QHBoxLayout();
	//ButtonLayout->addWidget(this->AddClusterButton);
	//ButtonLayout->addWidget(this->RunOperatorButton);
	//ButtonLayout->addWidget(this->ClearClusterButton);
	ButtonLayout->addWidget(this->ClusterFeaturesButton);
	ButtonLayout->addWidget(this->ShowObjectTables);
	ButtonLayout->addWidget(this->HideObjectTables);
	ButtonLayout->addWidget(this->ShowDistributionButton);

	//InfoLayout->addRow("Operator: ", this->OperatorList);

	this->MainLayout = new QVBoxLayout();
	this->HLayout = new QHBoxLayout();

	//this->MainLayout->addWidget(this->ClusterListView);
	//this->HLayout->addWidget(this->ClusterTableView->GetWidget());
	this->HLayout->addWidget(this->QVTKClusterTableView);
	this->HLayout->addLayout(InfoLayout);
	
	this->MainLayout->addLayout(ButtonLayout);
	this->MainLayout->addLayout(HLayout);
	//this->MainLayout->addWidget(this->QVTKClusterTableView);
	this->MainLayout->addLayout(this->HOperatorDisplayLayout);

	//this->MainLayout->addWidget(this->ClusterListView,0,0);
	//this->MainLayout->addWidget(this->ClusterTableView->GetWidget(),0,1);
	//this->MainLayout->addWidget(Operand1,1,0);
	//this->MainLayout->addWidget(OperatorList,1,1);
	//this->MainLayout->addWidget(Operand2,1,2);
	//this->MainLayout->addLayout(InfoLayout,2,0);


	this->setLayout(this->MainLayout);
	this->setWindowTitle(tr("Cluster Manager"));
	AnnotationLinkSetUp = false;
}

ClusterManager::~ClusterManager()
{
	delete this->QVTKClusterTableView;
	for (int i =0 ; this->ClusterObjectTables.size() > i ; i++)
	{
	delete this->ClusterObjectTables[i];
	this->ClusterObjectTables[i] = 0;
	}
	this->ClusterObjectTables.clear();
}

void ClusterManager::setManagerTitle(std::string newTitle)
{
	this->setWindowTitle(newTitle.c_str());
}

void ClusterManager::setClusteringModel(SelectiveClustering * newClusterModel)
{
	/*! 
	* Create link to SelectiveClustering 
	*/
	if(this->ClusterModel)
	{
		delete this->ClusterModel;
	}
	this->ClusterModel = newClusterModel;
	connect(this->ClusterModel, SIGNAL(ClusterChanged()), this, SLOT(ChangeInClusters()));
	connect(this->ClusterModel, SIGNAL(DataChanged()), this, SLOT(ChangeInData()));
}

void ClusterManager::setObjectSelection(ObjectSelection *ObjSelection)
{
	/*! 
	* allows for use of Object Selection
	*/
	this->LegacyObjectSelection = ObjSelection;
	connect(this->LegacyObjectSelection, SIGNAL(changed()), this, SLOT(ChangeInObjectSelection()));

}

void ClusterManager::SelectionToClusterModification()
{
	/*! 
	* Currently test of object Selection to Clusters
	*/
	std::set<long int> curSel = this->LegacyObjectSelection->getSelections();
	std::set< vtkIdType > sel = SelectionUtilities::ObjectSelectionToIDSet(curSel);
	this->ClusterModel->AddCluster(sel);
}

void ClusterManager::ClearClusters()
{
	/*! 
	* Clears from the cluster model signals/slots auto update displays
	*/
	this->ClusterModel->ClearClusters();
}

void ClusterManager::RemoveSelectedClusters()
{
	/*!
	* Gets a IDArray of selected clusters 
	* deletes Clusters fom the model
	* Objects remain
	*/
	vtkSmartPointer<vtkIdTypeArray> SelectedClusters = this->GetClusterTableSelections();
	this->ClusterModel->RemoveCluster(SelectedClusters);
	//this->ClusterModel->GetClusterTable()->Dump(16);
}

void ClusterManager::ShowClusterFeatures()
{
	/*!
	* Shows window with all clusters and average feature values
	*/
	if (!this->ClusterFeatureDialog)
	{
		this->ClusterFeatureDialog = new QvtkTableDialog();
		this->ClusterFeatureDialog->setTitle("Table of feature averages");
	}
	if (this->ClusterFeatureDialog->closed)
	{
		this->ClusterFeatureDialog = new QvtkTableDialog();
		this->ClusterFeatureDialog->setTitle("Table of feature averages");
	}
	this->ClusterFeatureDialog->UpdateView(this->ClusterModel->ClusterFeatureTable(),  this->ClusterModel->ClusterAnnotationLink);
}

void ClusterManager::ShowClusterObjectTables()
{
	/*!
	* Displays a table for each cluster
	* each row is an object and its features
	*/
	if (this->ClusterObjectTables.size() > 0)
	{
		this->CloseClusterObjectTables();
	}
	vtkSmartPointer<vtkIdTypeArray> SelectedClusters = this->GetClusterTableSelections();
	if (SelectedClusters->GetSize() > 0)
	{
		for (vtkIdType count = 0; count < SelectedClusters->GetSize(); count++)
		{
			vtkIdType key = SelectedClusters->GetValue(count);
			if (!this->ClusterModel->ValidKey(key))
			{
				continue;
			}
			QvtkTableDialog* newTable = new QvtkTableDialog();
			std::stringstream temp;
			temp<<"Cluster: " << key;
			newTable->setTitle(temp.str());
			newTable->UpdateView(this->ClusterModel->GetTableOfSelectedFromCluster(key), 
				this->ClusterModel->ObjectAnnotationLink);
			this->ClusterObjectTables.push_back(newTable);
		}
	}
	else
	{
		vtkIdType numClus = this->ClusterModel->NumberOfClusters();
		std::set< vtkIdType > clusterIDs =  this->ClusterModel->GetClusterIDs();
		std::set< vtkIdType >::iterator iter = clusterIDs.begin();
		for (; iter != clusterIDs.end(); iter++)
		{
			QvtkTableDialog* newTable = new QvtkTableDialog();
			std::stringstream temp;
			temp<<"Cluster: " << (*iter);
			newTable->setTitle(temp.str());
			newTable->UpdateView(this->ClusterModel->GetTableOfSelectedFromCluster(*iter), 
				this->ClusterModel->ObjectAnnotationLink);
			this->ClusterObjectTables.push_back(newTable);
		}
	}
	//Table Of all objects
	QvtkTableDialog* newTable = new QvtkTableDialog();
	newTable->setTitle("All Selected Objects");
	newTable->UpdateView(this->ClusterModel->GetTableOfAllSelected(), 
		this->ClusterModel->ObjectAnnotationLink);
	this->ClusterObjectTables.push_back(newTable);

}

void ClusterManager::CloseClusterObjectTables()
{
	/*!
	* Closes and deletes tables views of 
	* the objects in each Cluster 
	*/
	for (int tablesOpen =0 ; this->ClusterObjectTables.size() > tablesOpen ; tablesOpen++)
	{
		if (!this->ClusterObjectTables[tablesOpen]->closed)
		{
			this->ClusterObjectTables[tablesOpen]->close();
		}
		delete this->ClusterObjectTables[tablesOpen];
		this->ClusterObjectTables[tablesOpen] = 0;
	}
	this->ClusterObjectTables.clear();
}

vtkSmartPointer<vtkIdTypeArray> ClusterManager::GetClusterTableSelections()
{
	/*!
	* maps selected vtkTable rows to selected vtkIDS
	*/
	
	return this->QVTKClusterTableView->getSelectedObjects();
}
void ClusterManager::ChangeInClusters()
{
	/*! 
	* updates information when clusters change
	*/
	int numClust = (int) this->ClusterModel->NumberOfClusters();
	this->NumClusters->setNum(numClust);
	this->NumObjects->setNum((int) this->ClusterModel->GetNumberOfObjects());
	this->NumSelected->setNum((int) this->ClusterModel->GetNumberOfSelections());

	//display of clusters in a list
	/*QStringList ClusterList = this->ClusterModel->GetClusterIDsList();
	this->ClusterListView->clear();
	this->ClusterListView->addItems( ClusterList);*/
	this->QVTKClusterTableView->SetInputLink(this->ClusterModel->GetClusterTable(), this->ClusterModel->ClusterAnnotationLink);
	if (this->ClusterFeatureDialog)
	{
		this->ClusterFeatureDialog->UpdateView(this->ClusterModel->ClusterFeatureTable(),  this->ClusterModel->ClusterAnnotationLink);
	}
	//delete the tables
	if (this->ClusterObjectTables.size() > 0)
	{
		this->ShowClusterObjectTables();
	}

	QStringList ClusterList = this->ClusterModel->GetClusterIDsList();
	Operand1->clear();
	Operand2->clear();
	Operand1->addItems(ClusterList);
	Operand2->addItems(ClusterList);

}

void ClusterManager::ChangeInData()
{
	/*! 
	* updates views when sorce data changes
	*/
	int numClust = (int) this->ClusterModel->NumberOfClusters();
	this->NumClusters->setNum(numClust);
	this->NumObjects->setNum((int) this->ClusterModel->GetNumberOfObjects());
	this->NumSelected->setNum((int) this->ClusterModel->GetNumberOfSelections());

	this->QVTKClusterTableView->SetInputLink(this->ClusterModel->GetClusterTable(), this->ClusterModel->ClusterAnnotationLink);
	if (this->ClusterFeatureDialog)
	{
		this->ClusterFeatureDialog->UpdateView(this->ClusterModel->ClusterFeatureTable(),  this->ClusterModel->ClusterAnnotationLink);
	}
	if (this->ClusterObjectTables.size() > 0)
	{
		this->CloseClusterObjectTables();
	}
}

void ClusterManager::ChangeInObjectSelection()
{
	/*! 
	* Signal Slot interface for objectSelection
	*/
}

void ClusterManager::RunOperatorOnSelectedClusters()
{
	/*! 
	* Runs the Operator on the selected Clusters
	*/
	int index = this->OperatorList->currentIndex();
	std::set< vtkIdType > clusterIDs;
	vtkSmartPointer<vtkTable> selectedTable = vtkSmartPointer<vtkTable>::New();
	selectedTable->Initialize();
	
	//create vtkVariant and convert it to vtkIDType .ToTypeInt64(); :::::::::this is to convert the string to vtkIDType
	vtkVariant a ;
	vtkVariant b ;
	a = this->Operand1->currentText().toInt();
	b = this->Operand2->currentText().toInt();
	
	std::stringstream temp;
	
	switch(index)
	{
		case 0: //ADD
			clusterIDs =  this->ClusterModel->cluster_operator_ADD(a.ToTypeInt64(),b.ToTypeInt64());
			temp << a.ToTypeInt64() << " + " << b.ToTypeInt64();
			break;
		case 1: //SUBRACT
			clusterIDs =  this->ClusterModel->cluster_operator_SUBTRACT(a.ToTypeInt64(),b.ToTypeInt64());
			temp << a.ToTypeInt64() << " - " << b.ToTypeInt64();
			break;
		case 2: //AND
			clusterIDs =  this->ClusterModel->cluster_operator_AND(a.ToTypeInt64(),b.ToTypeInt64());
			temp << a.ToTypeInt64() << " And " << b.ToTypeInt64();
			break;
		case 3: //XOR
			clusterIDs =  this->ClusterModel->cluster_operator_XOR(a.ToTypeInt64(),b.ToTypeInt64());
			temp << a.ToTypeInt64() << " XOR " << b.ToTypeInt64();
			break;
		default: 
			std::cerr << "Incorrect Operator = " << index << std::endl;
			break;
	}
	
	//Display the Data in the console
	this->ClusterModel->CopySelectedIntoTable(clusterIDs, selectedTable);
	QvtkTableDialog* newTable = new QvtkTableDialog();
	//selectedTable->Dump(16);
	std::set< vtkIdType > tempClusterIds;
	index  = this->ActionType->currentIndex();
	switch(index)
	{
		case 0: 	

			
			newTable->setTitle(temp.str());
			newTable->UpdateView(selectedTable, this->ClusterModel->ObjectAnnotationLink);
			this->ClusterObjectTables.push_back(newTable);

			break;
		case 1: //Override I
			tempClusterIds = this->ClusterModel->SelectionFromCluster(a.ToTypeInt64());
			this->ClusterModel->RemoveSelectionFromCluster(a.ToTypeInt64(),tempClusterIds);
			this->ClusterModel->AddSelectionToCluster(a.ToTypeInt64(),clusterIDs);

			std::cout<< "Value of Cluster " << a.ToString() << std::endl;
			tempClusterIds = this->ClusterModel->SelectionFromCluster(a.ToTypeInt64());
			selectedTable->Initialize();
			this->ClusterModel->CopySelectedIntoTable(tempClusterIds, selectedTable);
			selectedTable->Dump(16);	

			break;
		case 2: //Override II
			tempClusterIds = this->ClusterModel->SelectionFromCluster(b.ToTypeInt64());
			this->ClusterModel->RemoveSelectionFromCluster(b.ToTypeInt64(),tempClusterIds);
			this->ClusterModel->AddSelectionToCluster(b.ToTypeInt64(),clusterIDs);


			std::cout<< "Value of Cluster " << b.ToString() << std::endl;
			tempClusterIds = this->ClusterModel->SelectionFromCluster(b.ToTypeInt64());
			selectedTable->Initialize();
			this->ClusterModel->CopySelectedIntoTable(tempClusterIds, selectedTable);
			selectedTable->Dump(16);	

			break;
		default: 
			std::cerr << "Incorrect Action Type = " << index << std::endl;
			break;
	}
}
void ClusterManager::ShowDistribution()
{
	
//======================================================================================================================
	cout<<"I have reached into the ShowDistribution Function \n";
	this->ClusterModel->GetClusterTable()->Dump(16);
	
	vtkSmartPointer<vtkTable> distributionTable = vtkSmartPointer<vtkTable>::New();
	vtkSmartPointer<vtkChartPie> PieChart = vtkSmartPointer<vtkChartPie>::New();
	vtkSmartPointer<vtkColorSeries> colorSeries = vtkSmartPointer<vtkColorSeries>::New();
	vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();

	vtkSmartPointer<vtkIntArray> DistributionArray = vtkSmartPointer<vtkIntArray>::New();
	DistributionArray->SetName("Number of Object");
	vtkSmartPointer<vtkStringArray> ClusterIdArray = vtkSmartPointer<vtkStringArray>::New();
	ClusterIdArray->SetName("Cluster ID's");
	
	vtkTable *  clusterTable= this->ClusterModel->GetClusterTable();
	cout<< "this is clusterTable: \n";
	clusterTable->Dump();
	
	vtkIdType rowCount = this->ClusterModel->GetClusterTable()->GetNumberOfRows();
	for (vtkIdType row = 0; row != rowCount; row++)
	{
		DistributionArray->InsertNextValue(clusterTable->GetValue(row, 1).ToInt());
		ClusterIdArray->InsertNextValue(clusterTable->GetValue(row,0).ToString());
	}

	distributionTable->AddColumn(ClusterIdArray);
	distributionTable->AddColumn(DistributionArray);
	cout<< "This is the distributionTable to draw the PieChart: \n";
	distributionTable->Dump(16);
	
	colorSeries->SetColorScheme(vtkColorSeries::WARM);
	vtkPlotPie * pie = vtkPlotPie::SafeDownCast(PieChart->AddPlot(0));
	pie->SetColorSeries(colorSeries);

	#if VTK_MAJOR_VERSION <= 5
	  pie->SetInput(distributionTable);
	#else
	  pie->SetInputData(distributionTable);
	#endif
	 pie->SetInputArray(0,"Number of Object");
	 pie->SetLabels(ClusterIdArray);
	
	PieChart->SetTitle("Cluster Distribution");
	PieChart->SetShowLegend(true);

	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
	view->GetRenderWindow()->SetSize(600,350);
	view->GetScene()->AddItem(PieChart);
	view->GetRenderWindow()->SetMultiSamples(0);
	view->GetInteractor()->Initialize();
	view->GetInteractor()->Start();
	
}
