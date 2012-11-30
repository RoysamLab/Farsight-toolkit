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
#include "CellTraceModel.h"
#include "CellTrace.h"
CellTraceModel::CellTraceModel()
{	
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->ColumnSelection = new ObjectSelection();
	
	this->CellClusterSelection = new SelectiveClustering();
	//this->CellClusterManager = new ClusterManager();
	//this->CellClusterManager->setClusteringModel(this->CellClusterSelection );
	//this->CellClusterManager->setObjectSelection(this->Selection);
	//this->CellClusterManager->setManagerTitle("Cell Cluster Manager");

	this->Cells.clear();
	this->graphVisualize = new GraphWindow();
	this->AdditionalHeaders.clear();
}
CellTraceModel::CellTraceModel(std::map< int ,CellTrace*> Cells)
{	
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->ColumnSelection = new ObjectSelection();
	this->graphVisualize = NULL;
	this->setCells(Cells);
}
CellTraceModel::~CellTraceModel()
{	
	delete this->Selection;
  this->Selection = NULL;
	delete this->ColumnSelection;
  this->ColumnSelection = NULL;

  if(this->graphVisualize != NULL)
    {
    delete this->graphVisualize;
    this->graphVisualize = NULL;
    }
  //delete this->CellClusterManager;
  delete this->CellClusterSelection;
}
void CellTraceModel::setCells(std::map< int ,CellTrace*> Cells)
{
	this->Cells.clear();
	this->Cells = Cells;
	//this->CellIDLookupMAP.clear();
	this->SyncModel();
}
void CellTraceModel::SetupHeaders()
{
	this->headers.clear();
	this->headers.push_back("Root Trace");

	this->headers.push_back("Soma X Pos");
	this->headers.push_back("Soma Y Pos");
	this->headers.push_back("Soma Z Pos");

	this->headers.push_back("Width X");
	this->headers.push_back("Height Y");
	this->headers.push_back("Depth Z");

	this->headers.push_back("Soma Radii");
	this->headers.push_back("Soma Surface Area");
	this->headers.push_back("Soma Volume");

	this->headers.push_back("Skewness X");
	this->headers.push_back("Skewness Y");
	this->headers.push_back("Skewness Z");
	this->headers.push_back("Euclidean Skewness");

	this->headers.push_back("Segments");
	this->headers.push_back("Stems");
	this->headers.push_back("Branching Stems");
	this->headers.push_back("Branch Pt");
	this->headers.push_back("Bifurcations");
	this->headers.push_back("Leaf Nodes");

	this->headers.push_back("Min Diameter");
	this->headers.push_back("Ave Diameter");
	this->headers.push_back("Max Diameter");

	this->headers.push_back("Min Diameter Power");
	this->headers.push_back("Ave Diameter Power");
	this->headers.push_back("Max Diameter Power");

	this->headers.push_back("Tot Volume");
	this->headers.push_back("Min Segment Volume");
	this->headers.push_back("Ave Segment Volume");
	this->headers.push_back("Max Segment Volume");
	this->headers.push_back("Tot Surface Area");
	this->headers.push_back("Min Segment Surface Area");
	this->headers.push_back("Ave Segment Surface Area");
	this->headers.push_back("Max Segment Surface Area");
	//this->headers.push_back("Tot Segment Section Area");
	this->headers.push_back("Min Segment Section Area");
	this->headers.push_back("Ave Segment Section Area");
	this->headers.push_back("Max Segment Section Area");

	this->headers.push_back("Min Burk Taper");
	this->headers.push_back("Ave Burk Taper");
	this->headers.push_back("Max Burk Taper");

	this->headers.push_back("Min Hillman Taper");
	this->headers.push_back("Ave Hillman Taper");
	this->headers.push_back("Max Hillman Taper");

	this->headers.push_back("Tot Euclidean Length");
	this->headers.push_back("Ave Segment Euclidean Length");
	this->headers.push_back("Tot Path Length");
	this->headers.push_back("Min Segment Path Length");
	this->headers.push_back("Ave Segment Path Length");
	this->headers.push_back("Max Segment Path Length");

	this->headers.push_back("Min Stem Distance");
	this->headers.push_back("Ave Stem Distance");
	this->headers.push_back("Max Stem Distance");

	this->headers.push_back("Min Contraction");
	this->headers.push_back("Ave Contraction");
	this->headers.push_back("Max Contraction");

	this->headers.push_back("Tot Fragmentation");
	this->headers.push_back("Min Fragmentation");
	this->headers.push_back("Ave Fragmentation");
	this->headers.push_back("Max Fragmentation");

	this->headers.push_back("Min Daughter Ratio");
	this->headers.push_back("Ave Daughter Ratio");
	this->headers.push_back("Max Daughter Ratio");

	this->headers.push_back("Min Parent Daughter Ratio");
	this->headers.push_back("Ave Parent Daughter Ratio");
	this->headers.push_back("Max Parent Daughter Ratio");

	this->headers.push_back("Min Daughter Length Ratio");
	this->headers.push_back("Ave Daughter Length Ratio");
	this->headers.push_back("Max Daughter Length Ratio");

	this->headers.push_back("Min Partition Asymmetry");
	this->headers.push_back("Ave Partition Asymmetry");
	this->headers.push_back("Max Partition Asymmetry");

	this->headers.push_back("Min Rall Power");
	this->headers.push_back("Ave Rall Power");
	this->headers.push_back("Max Rall Power");

	this->headers.push_back("Min Pk");
	this->headers.push_back("Ave Pk");
	this->headers.push_back("Max Pk");

	this->headers.push_back("Min Pk Classic");
	this->headers.push_back("Ave Pk Classic");
	this->headers.push_back("Max Pk Classic");

	this->headers.push_back("Min Pk 2");
	this->headers.push_back("Ave Pk 2");
	this->headers.push_back("Max Pk 2");

	this->headers.push_back("Min Azimuth");
	this->headers.push_back("Ave Azimuth");
	this->headers.push_back("Max Azimuth");
	this->headers.push_back("Min Elevation");
	this->headers.push_back("Ave Elevation");
	this->headers.push_back("Max Elevation");

	this->headers.push_back("Min Local Bifurcation Amp");
	this->headers.push_back("Ave Local Bifurcation Amp");
	this->headers.push_back("Max Local Bifurcation Amp");
	this->headers.push_back("Min Local Bifurcation Tilt");
	this->headers.push_back("Ave Local Bifurcation Tilt");
	this->headers.push_back("Max Local Bifurcation Tilt");
	this->headers.push_back("Min Local Bifurcation Torque");
	this->headers.push_back("Ave Local Bifurcation Torque");
	this->headers.push_back("Max Local Bifurcation Torque");
	
	this->headers.push_back("Min Remote Bifurcation Amp");
	this->headers.push_back("Ave Remote Bifurcation Amp");
	this->headers.push_back("Max Remote Bifurcation Amp");
	this->headers.push_back("Min Remote Bifurcation Tilt");
	this->headers.push_back("Ave Remote Bifurcation Tilt");
	this->headers.push_back("Max Remote Bifurcation Tilt");
	this->headers.push_back("Min Remote Bifurcation Torque");
	this->headers.push_back("Ave Remote Bifurcation Torque");
	this->headers.push_back("Max Remote Bifurcation Torque");

	this->headers.push_back("Min Leaf Level");
	this->headers.push_back("Min Leaf Path Length");
	this->headers.push_back("Ave Leaf Level");
	this->headers.push_back("Ave Leaf Path Length");
	this->headers.push_back("Max Leaf Level");
	this->headers.push_back("Max Leaf Path Length");

	this->headers.push_back("Tot Terminal Compartments");
	this->headers.push_back("Min Terminal Compartments");
	this->headers.push_back("Ave Terminal Compartments");
	this->headers.push_back("Max Terminal Compartments");

	this->headers.push_back("Diameter Threshold Min");
	this->headers.push_back("Ave Diameter Threshold");
	this->headers.push_back("Diameter Threshold Max");
	this->headers.push_back("Terminal Parent Diameter Min");
	this->headers.push_back("Ave Terminal Parent Diameter");
	this->headers.push_back("Terminal Parent Diameter Max");

	this->headers.push_back("Min Terminal Hillman Thresh"); //extend for all branches
	this->headers.push_back("Ave Terminal Hillman Thresh");
	this->headers.push_back("Max Terminal Hillman Thresh");

	this->headers.push_back("Min Branch Pt to Soma Euclidean Distance");
	this->headers.push_back("Ave Branch Pt to Soma Euclidean Distance");
	this->headers.push_back("Max Branch Pt to Soma Euclidean Distance");

	this->headers.push_back("Min Tip to Soma Euclidean Distance");
	this->headers.push_back("Ave Tip to Soma Euclidean Distance");
	this->headers.push_back("Max Tip to Soma Euclidean Distance");

	this->headers.push_back("Overall Magnitude of Tips");
	this->headers.push_back("Overall Azimuth of Tips");
	this->headers.push_back("Overall Elevation of Tips");


	this->headers.push_back("Trace File");

	/*this->headers.push_back("Prediction");
	this->headers.push_back("Confidence");*/
	this->headers.push_back("Distance to Device");
	//this->headers.push_back("Label");

	int size = this->AdditionalHeaders.size();
	for (int k = 0; k < size; k++)
	{	
		this->headers.push_back(QString(this->AdditionalHeaders[k].c_str()));
	}
	
	int numHeaders = (int)this->headers.size();
	//std::cout<<numHeaders << "\t features computed\n";
	
	vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( this->headers.at(i).toStdString().c_str() );
		this->DataTable->AddColumn(column);
    }
}
void CellTraceModel::SyncModel()
{	
	int size = this->AdditionalHeaders.size();
	this->DataTable->Initialize();	
	this->Selection->clear();
	//this->CellIDLookupMAP.clear();
	this->SetupHeaders();
	for (CellIDLookupIter = this->Cells.begin();  CellIDLookupIter != this->Cells.end(); CellIDLookupIter ++)
	{
		if (size == 0)
		{
			this->DataTable->InsertNextRow((*CellIDLookupIter).second->DataRow());
			//this->CellIDLookupMAP[ this->Cells[i]->rootID() ] = this->Cells[i];
		}
		else
		{
			this->DataTable->InsertNextRow((*CellIDLookupIter).second->GetExtendedDataRow(this->AdditionalHeaders));
			//this->CellIDLookupMAP[ this->Cells[i]->rootID() ] = this->Cells[i];
		}
	}
	this->CellClusterSelection->SetObjectTable(this->DataTable);
	//this->CellClusterManager->setVisible(true);
}


vtkSmartPointer<vtkTable> CellTraceModel::getDataTable()
{
	return this->DataTable;
}

void CellTraceModel::setDataTable(vtkSmartPointer<vtkTable> table)
{
	this->DataTable = table;
	this->Selection->clear();
}

vtkSmartPointer<vtkTable> CellTraceModel::getCellBoundsTable()
{
	vtkSmartPointer<vtkTable> CellBoundsTable = vtkSmartPointer<vtkTable>::New();
	CellBoundsTable->Initialize();
	std::vector<QString> BoundsHeaders;
	BoundsHeaders.clear();
	BoundsHeaders.push_back("Cell ID");
	BoundsHeaders.push_back("Soma X Pos");
	BoundsHeaders.push_back("Soma Y Pos");
	BoundsHeaders.push_back("Soma Z Pos");
	BoundsHeaders.push_back("Min X");
	BoundsHeaders.push_back("Max X");
	BoundsHeaders.push_back("Min Y");
	BoundsHeaders.push_back("Max Y");
	BoundsHeaders.push_back("Min Z");
	BoundsHeaders.push_back("Max Z");
	BoundsHeaders.push_back("Skewness X");
	BoundsHeaders.push_back("Skewness Y");
	BoundsHeaders.push_back("Skewness Z");

	int numHeaders = (int)BoundsHeaders.size();
	vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( BoundsHeaders.at(i).toStdString().c_str() );
		CellBoundsTable->AddColumn(column);
    }
	this->CellIDLookupIter = this->Cells.begin();
	for(; CellIDLookupIter != this->Cells.end(); CellIDLookupIter++)
	{
		CellTrace* current = (*CellIDLookupIter).second;
		CellBoundsTable->InsertNextRow(current->BoundsRow());
	}
	return CellBoundsTable;
}

ObjectSelection * CellTraceModel::GetObjectSelection()
{
	return this->Selection;
}
ObjectSelection * CellTraceModel::GetObjectSelectionColumn()
{
	return this->ColumnSelection;
}
void CellTraceModel::SelectByRootTrace(std::vector<TraceLine*> roots)
{
	//this->Selection->clear();
	std::set<long int> ID;
	for (int i = 0; i < roots.size(); i++)
	{
		ID.insert((long)roots.at(i)->GetId());
		//std::cout<< "Root ID selected: " << roots.at(i)->GetId() << std::endl; //print out the root traces selected
	}
	this->Selection->select(ID);
}

void CellTraceModel::SelectByIDs(std::vector<int> IDs)
{
	//this->Selection->clear();
	std::set<long int> ID;
	std::vector<int>::iterator IDs_iterator;
	
	//std::cout << "IDs.size(): " << IDs.size() << std::endl;
	for (IDs_iterator = IDs.begin(); IDs_iterator != IDs.end(); IDs_iterator++)
	{
		ID.insert(*IDs_iterator);
	}
	Selection->select(ID);
}

std::vector<CellTrace*> CellTraceModel::getCells(std::vector<long> IDs)
{
	std::vector<CellTrace*> cells;
	std::vector<long>::iterator IDs_iterator;
	
	for (IDs_iterator = IDs.begin(); IDs_iterator != IDs.end(); IDs_iterator++)
		cells.push_back(Cells[*IDs_iterator]);
	
	return cells;
}

std::set<long int> CellTraceModel::GetSelectedContinuousIDs()
{
	std::set<long int> ids;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		long int id = (long int) *it;

		bool found = false;
		long int j = 0;

		std::map< int ,CellTrace*>::iterator mapIter = this->Cells.begin();
		while(!found && mapIter != this->Cells.end())
		{
			if (id == mapIter->first)
			{
				ids.insert(j);
				found = true;
			}else
			{
				j++;
			}
			mapIter++;
		}//end while !found

	}//end for selected
	return ids;
}

std::set<long int> CellTraceModel::GetSelectedIDs()
{
	std::set<long int> allSelectedIDs, nextIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		this->CellIDLookupIter = this->Cells.find((int) *it);
		if (this->CellIDLookupIter != this->Cells.end())
		{
			nextIDs = ((*this->CellIDLookupIter).second)->TraceIDsInCell();
			std::set<long>::iterator k;
			for (k = nextIDs.begin(); k != nextIDs.end(); k++)
			{
				allSelectedIDs.insert(*k);
			}
		}
	}//end for selected
	return allSelectedIDs;
}

std::vector<TraceLine*> CellTraceModel::GetSelectedTraces()
{
	/*!
	*
	*/
	std::vector<TraceLine*> TracesReturned, tempset;
	TracesReturned.clear();
	tempset.clear();
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		this->CellIDLookupIter = this->Cells.find((int) *it);
		if (this->CellIDLookupIter != this->Cells.end())
		{
			tempset = ((*this->CellIDLookupIter).second)->getSegments();
			TracesReturned.insert(TracesReturned.end(), tempset.begin(), tempset.end());
		}
	}//end for selected
	return TracesReturned;
}

std::vector<CellTrace*> CellTraceModel::GetSelectedCells()
{
	std::vector<CellTrace*> selectedCell;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		this->CellIDLookupIter = this->Cells.find((int) *it);
		if (this->CellIDLookupIter != this->Cells.end())
		{
			selectedCell.push_back((*this->CellIDLookupIter).second);
		}
	}//end for selected
	return selectedCell;
}
unsigned int CellTraceModel::getCellCount()
{
	return (unsigned int) this->Cells.size();
}
CellTrace * CellTraceModel::GetCell( int i)
{
	CellTrace* currentCell;
	CellIDLookupIter = this->Cells.find(i);
	if (CellIDLookupIter != this->Cells.end())
	{
		currentCell = (*CellIDLookupIter).second;
	}
	else
	{
		currentCell = (*this->Cells.rbegin()).second;
	}
	this->Selection->select(currentCell->rootID());
	return currentCell;
}
CellTrace * CellTraceModel::GetCellNoSelection( int i)
{
	CellTrace* currentCell;
	CellIDLookupIter = this->Cells.find(i);
	if (CellIDLookupIter != this->Cells.end())
	{
		currentCell = (*CellIDLookupIter).second;
	}
	else
	{
		currentCell = (*this->Cells.rbegin()).second;
	}
	return currentCell;
}

std::map< int ,CellTrace*>::iterator CellTraceModel::GetCelliterator()
{
	/*!
	*
	*/
	return this->Cells.begin();
}
std::map< int ,CellTrace*>::iterator CellTraceModel::GetCelliteratorEnd()
{
	/*!
	*
	*/
	return this->Cells.end();
}

void CellTraceModel::WriteCellCoordsToFile(const char *fileName)
{
	fstream outputTxt;
	outputTxt.open(fileName, fstream::out);
	this->CellIDLookupIter = this->Cells.begin();
	for(; CellIDLookupIter != this->Cells.end(); CellIDLookupIter++)
	{
		CellTrace* current = (*CellIDLookupIter).second;
		outputTxt << current->BasicFeatureString() << "\n";
	}
	outputTxt.close();
	std::cout << "file written\n";
}
void CellTraceModel::createCellToCellGraph()
{
	std::map< unsigned int, std::vector<double> > centroidMap;
	this->CellIDLookupIter = this->Cells.begin();
	for(; CellIDLookupIter != this->Cells.end(); CellIDLookupIter++)
	{
		CellTrace* current = (*CellIDLookupIter).second;
		unsigned int id = current->rootID();
		std::vector<double> cellCoord;
		cellCoord.push_back(current->somaX);
		cellCoord.push_back(current->somaY);
		cellCoord.push_back(current->somaZ);
		centroidMap[id] = cellCoord;
	}//end for cell size
	kNearestObjects<3>* KNObj = new kNearestObjects<3>(centroidMap);
	//KNObj->setFeatureTable(this->getCellBoundsTable());
	std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	kNeighborIDs = KNObj->k_nearest_neighbors_All(4, 0, 0);	//guess at parameters, dist, no classes 
	//////////////////////////////////////////////////////
	/*std::string full_string;
	std::stringstream ss1;
	ss1 << 4;
	full_string = "D(k=" + ss1.str() + ",class=all)" ;
	DataTable->RemoveColumnByName(full_string.c_str());
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName(full_string.c_str());
	column->SetNumberOfValues((int)DataTable->GetNumberOfRows());
	DataTable->AddColumn(column);
	for(int row=0; row<(int)DataTable->GetNumberOfRows(); ++row)
	{
		DataTable->SetValueByName(row, full_string.c_str(), 0);
	}
	for(int i=0; i < (int)kNeighborIDs.size(); ++i)
	{
		int Id = kNeighborIDs.at(i).at(0).first;
		double avg_dist = average(kNeighborIDs.at(i));
		for(int row=0; row<(int)DataTable->GetNumberOfRows(); ++row)
		{
			if(DataTable->GetValue(row,0).ToInt() == Id)
			{
				DataTable->SetValueByName(row, full_string.c_str(), vtkVariant(avg_dist));
				break;
			}
		}
	}*/
	//////////////////////////////////////////////////
	vtkSmartPointer<vtkTable> graphTable = KNObj->vectorsToGraphTable(kNeighborIDs);
	this->graphVisualize->setModels(this->getDataTable(), this->GetObjectSelection());
	this->graphVisualize->SetGraphTable(graphTable, "Source", "Target", "Distance", "Soma X Pos", "Soma Y Pos", "Soma Z Pos");
	this->graphVisualize->ShowGraphWindow();
  delete KNObj;
  KNObj = NULL;
}

double CellTraceModel::average(std::vector< std::pair<unsigned int, double> > ID)
{
	double dist = 0;
	for(int i=1; i<(int)ID.size(); ++i)
	{
		dist += ID.at(i).second;
	}
	double average = dist/(int)(ID.size()-1);
	return average;
}
int CellTraceModel::AddNewFeatureHeader(std::string NewHeader)
{
	for (int k = 0; k < this->AdditionalHeaders.size(); k++)
	{
		if (NewHeader == this->AdditionalHeaders[k])
		{
			return k;
		}
	}
	this->AdditionalHeaders.push_back(NewHeader);
	this->headers.push_back(QString(NewHeader.c_str()));
	return -1;
}

void CellTraceModel::CloseClusterManager()
{
	//this->CellClusterManager->close();
}

SelectiveClustering *CellTraceModel::GetCellSelectiveClustering()
{
	return CellClusterSelection;
}
