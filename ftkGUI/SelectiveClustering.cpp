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
	this->ClusterVtkViewUpdater = vtkSmartPointer<vtkViewUpdater>::New();

	this->ObjectAnnotationLink = vtkSmartPointer<vtkAnnotationLink>::New();
	this->ObjectVtkViewUpdater = vtkSmartPointer<vtkViewUpdater>::New();
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
		this->RemoveRowFromClusterTable(key);
		emit ClusterChanged();
		return true;
	}
	//not found in map
	return false;
}

bool SelectiveClustering::RemoveCluster(vtkIdTypeArray *SelectedClusters)
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
			this->RemoveRowFromClusterTable(key);
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
	this->ObjectTable = InputObjectTable;
	this->NumberOfObjects = rows;
	for ( vtkIdType currow = 0; currow <= this->NumberOfObjects; currow++ )
	{
		vtkIdType rowObjId = this->ObjectTable->GetValue( currow, 0 ).ToTypeInt64();
		this->ObjectTableIDMap[ rowObjId ] = currow;
	}
	emit DataChanged();
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
	vtkVariantArray * NewRow = vtkVariantArray::New();
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

	this->ClusterFeaturesButton = new QPushButton(" show Cluster ");
	connect(this->ClusterFeaturesButton, SIGNAL(clicked()), this, SLOT(ShowClusterFeatures()));
	
	//this->ClusterListView = new QListWidget(this);
	this->ClusterTableView = vtkSmartPointer<vtkQtTableView>::New();
	
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
	InfoLayout->addRow(this->RemoveClusterButton);

	HOperatorDisplayLayout= new QHBoxLayout();
	this->HOperatorDisplayLayout->addWidget(Operand1);
	this->HOperatorDisplayLayout->addWidget(OperatorList);
	this->HOperatorDisplayLayout->addWidget(ActionType);
	this->HOperatorDisplayLayout->addWidget(Operand2);
	this->HOperatorDisplayLayout->addWidget(this->RunOperatorButton);
	//InfoLayout->addRow(this->HOperatorDisplayLayout);

	QHBoxLayout* ButtonLayout = new QHBoxLayout();
	ButtonLayout->addWidget(this->AddClusterButton);
	//ButtonLayout->addWidget(this->RunOperatorButton);
	ButtonLayout->addWidget(this->ClearClusterButton);
	ButtonLayout->addWidget(this->ClusterFeaturesButton);

	//InfoLayout->addRow("Operator: ", this->OperatorList);

	this->MainLayout = new QVBoxLayout();
	this->HLayout = new QHBoxLayout();

	//this->MainLayout->addWidget(this->ClusterListView);
	this->HLayout->addWidget(this->ClusterTableView->GetWidget());
	this->HLayout->addLayout(InfoLayout);
	
	this->MainLayout->addLayout(ButtonLayout);
	this->MainLayout->addLayout(HLayout);
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

void ClusterManager::setClusteringModel(SelectiveClustering * newClusterModel)
{
	/*! 
	* Create link to SelectiveClustering 
	*/
	this->ClusterModel = newClusterModel;
	connect(this->ClusterModel, SIGNAL(ClusterChanged()), this, SLOT(ChangeInClusters()));
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
	std::set< vtkIdType > sel = this->ObjectSelectionToIDSet();
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
	* 
	*/
	vtkIdTypeArray * SelectedClusters = this->GetClusterTableSelections();
	this->ClusterModel->RemoveCluster(SelectedClusters);
	//this->ClusterModel->GetClusterTable()->Dump(16);
}

void ClusterManager::ShowClusterFeatures()
{
	this->ClusterModel->ClusterFeatureTable();
}

vtkIdTypeArray * ClusterManager::GetClusterTableSelections()
{
	/*!
	* maps selected vtkTable rows to selected vtkIDS
	*/
	vtkIdTypeArray * SelectedClusters = vtkIdTypeArray::New() ;
	vtkIdTypeArray * SelectedRows = vtkIdTypeArray::New() ;

	this->ClusterTableView->GetSelectedItems(SelectedRows);
	vtkSmartPointer<vtkTable> table = this->ClusterModel->GetClusterTable();
	//table->Dump(16);
	for (vtkIdType count = 0; count < SelectedRows->GetSize(); count++)
	{
		vtkIdType row = SelectedRows->GetValue(count);
		vtkVariant value = table->GetValue(row, 0);
		SelectedClusters->InsertNextValue(value.ToTypeInt64());
	}
	return SelectedClusters;
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

	this->ClusterTableView->RemoveAllRepresentations();
	this->ClusterTableView->Update();
	vtkQtTableModelAdapter adapt(this->ClusterModel->GetClusterTable());
	this->ClusterTableView->SetRepresentationFromInput(adapt.GetVTKDataObject());
	if(!AnnotationLinkSetUp)
	{
		this->ClusterTableView->GetRepresentation()->SetAnnotationLink(this->ClusterModel->ClusterAnnotationLink);
		this->ClusterTableView->GetRepresentation()->SetSelectionType(vtkSelectionNode::PEDIGREEIDS);
		this->ClusterTableView->GetRepresentation()->SetSelectionArrayName("Cluster ID");
		/*this->ClusterModel->ClusterVtkViewUpdater->AddView(this->ClusterTableView);
		this->ClusterModel->ClusterVtkViewUpdater->AddAnnotationLink(this->ClusterModel->ClusterAnnotationLink);*/

		this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New(); 
		this->selectionCallback->SetClientData(this);
		this->selectionCallback->SetCallback ( SelectionCallbackFunction);
		vtkAnnotationLink *link = this->ClusterTableView->GetRepresentation()->GetAnnotationLink();
		link->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);
		AnnotationLinkSetUp = true;
	}
	this->ClusterTableView->Update();

	//display of clusters in a list
	/*QStringList ClusterList = this->ClusterModel->GetClusterIDsList();
	this->ClusterListView->clear();
	this->ClusterListView->addItems( ClusterList);*/

	QStringList ClusterList = this->ClusterModel->GetClusterIDsList();
	Operand1->clear();
	Operand2->clear();
	Operand1->addItems(ClusterList);
	Operand2->addItems(ClusterList);

}

void ClusterManager::ChangeInObjectSelection()
{
	/*! 
	* Signal Slot interface for objectSelection
	*/
}

std::set< vtkIdType > ClusterManager::ObjectSelectionToIDSet()
{
	/*! 
	* convert objectSelection into form selective
	* clustering can use for operations
	*/
	std::set<long int> curSel = this->LegacyObjectSelection->getSelections();
	std::set<long int>::iterator iter = curSel.begin();
	std::set< vtkIdType > Selection;
	for (; iter != curSel.end(); iter++)
	{
		vtkIdType id = (vtkIdType) (*iter);
		Selection.insert(id);
	}
	return Selection;
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
	
	switch(index)
	{
		case 0: //ADD
			clusterIDs =  this->ClusterModel->cluster_operator_ADD(a.ToTypeInt64(),b.ToTypeInt64());
			break;
		case 1: //SUBRACT
			clusterIDs =  this->ClusterModel->cluster_operator_SUBTRACT(a.ToTypeInt64(),b.ToTypeInt64());
			break;
		case 2: //AND
			clusterIDs =  this->ClusterModel->cluster_operator_AND(a.ToTypeInt64(),b.ToTypeInt64());
			break;
		case 3: //XOR
			clusterIDs =  this->ClusterModel->cluster_operator_XOR(a.ToTypeInt64(),b.ToTypeInt64());
			break;
		default: 
			std::cerr << "Incorrect Operator = " << index << std::endl;
			break;
	}
	
	//Display the Data in the console
	this->ClusterModel->CopySelectedIntoTable(clusterIDs, selectedTable);
	selectedTable->Dump(16);
	std::set< vtkIdType > tempClusterIds;
	index  = this->ActionType->currentIndex();
	switch(index)
	{
		case 0: //Display
			//tempClusterIds = this->ClusterModel->SelectionFromCluster(a.ToTypeInt64());
			//this->ClusterModel->RemoveSelectionFromCluster(a.ToTypeInt64(),tempClusterIds);
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
void ClusterManager::SelectionCallbackFunction(vtkObject *caller, unsigned long eventId, void *clientData, void *callData)
{
	/*! 
	* 
	*/ 
	//std::cout<< "in callback function \n";
	vtkAnnotationLink * annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection * selection = annotationLink->GetCurrentSelection();
	ClusterManager * ClusMan = (ClusterManager*)clientData;
	vtkSelection * TableRowSelection = vtkSelection::New();
	vtkSelectionNode* vertices = NULL;
	//vtkSelectionNode* edges = NULL;

	if( selection->GetNode(0))
	{
		if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(0);
		}/*
		else if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
		{
			edges = selection->GetNode(0);
		}*/
	}

	if( selection->GetNode(1))
	{
		if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(1);
		}/*
		else if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
		{
			edges = selection->GetNode(1);
		}*/
	}
	
	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
		vtkIdType numTuples = vertexList->GetNumberOfTuples();
		if( vertexList != NULL && numTuples > 0)
		{
			//std::cout<< "number of selections: " << numTuples << "\nselected:\n";
			for( vtkIdType i = 0; i < numTuples; i++)
			{
				vtkIdType value = vertexList->GetValue(i);
				//std::cout<< value << "\n";
				//selectedRows->InsertTuple(
			}
			//TableRowSelection->GetNode(0)->SetSelectionList(vertexList);
		}
		//ClusMan->ClusterTableView->GetRepresentation()->UpdateSelection(selection);
	}//end verticies != null
	ClusMan->ClusterTableView->GetRepresentation()->GetAnnotationLink()->SetCurrentSelection(selection);
	ClusMan->ClusterTableView->Update();

}