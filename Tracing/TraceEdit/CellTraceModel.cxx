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
	this->Cells.clear();
	this->graphVisualize = new GraphWindow();
	this->AdditionalHeaders.clear();
}
CellTraceModel::CellTraceModel(std::vector<CellTrace*> Cells)
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
}
void CellTraceModel::setCells(std::vector<CellTrace*> Cells)
{
	this->Cells.clear();
	this->Cells = Cells;
	this->SyncModel();
}
void CellTraceModel::SetupHeaders()
{
	this->headers.clear();
	this->headers.push_back("Root Trace");

	this->headers.push_back("Width X");
	this->headers.push_back("Height Y");
	this->headers.push_back("Depth Z");

	this->headers.push_back("Soma X Pos");
	this->headers.push_back("Soma Y Pos");
	this->headers.push_back("Soma Z Pos");
	this->headers.push_back("Soma Radii");
	this->headers.push_back("Soma Volume");
	this->headers.push_back("Soma Surface Area");

	this->headers.push_back("Skewness X");
	this->headers.push_back("Skewness Y");
	this->headers.push_back("Skewness Z");
	this->headers.push_back("Euclidean Skewness");

	this->headers.push_back("Segments");
	this->headers.push_back("Stems");
	this->headers.push_back("Branch Pt");
	this->headers.push_back("Leaf Nodes");

	this->headers.push_back("Min Diameter");
	this->headers.push_back("Average Diameter");
	this->headers.push_back("Max Diameter");

	this->headers.push_back("Min Diameter Power");
	this->headers.push_back("Average Diameter Power");
	this->headers.push_back("Max Diameter Power");

	this->headers.push_back("Total Volume");
	this->headers.push_back("Min Segment Volume");
	this->headers.push_back("Max Segment Volume");
	this->headers.push_back("Average Segment Volume");
	this->headers.push_back("Total Surface Area");
	this->headers.push_back("Min Segment Surface Area");
	this->headers.push_back("Average Segment Surface Area");
	this->headers.push_back("Max Segment Surface Area");
	this->headers.push_back("Total Segment Section Area");
	this->headers.push_back("Min Segment Section Area");
	this->headers.push_back("Max Segment Section Area");

	this->headers.push_back("Min Burk Taper");
	this->headers.push_back("Average Burk Taper");
	this->headers.push_back("Max Burk Taper");

	this->headers.push_back("Min Hillman Taper");
	this->headers.push_back("Average Hillman Taper");
	this->headers.push_back("Max Hillman Taper");

	this->headers.push_back("Total Euclidian Length");
	this->headers.push_back("Average Segment Euclidian Length");
	this->headers.push_back("Total Path Length");
	this->headers.push_back("Average Segment Path Length");

	this->headers.push_back("Min Stem Distance");
	this->headers.push_back("Average Stem Distance");
	this->headers.push_back("Max Stem Distance");

	this->headers.push_back("Min Contraction");
	this->headers.push_back("Average Contraction");
	this->headers.push_back("Max Contraction");

	this->headers.push_back("Total Fragmentation");
	this->headers.push_back("Min Fragmentation");
	this->headers.push_back("Average Fragmentation");
	this->headers.push_back("Max Fragmentation");

	this->headers.push_back("Min Daughter Ratio");
	this->headers.push_back("Average Daughter Ratio");
	this->headers.push_back("Max Daughter Ratio");

	this->headers.push_back("Min Parent Daughter Ratio");
	this->headers.push_back("Average Parent Daughter Ratio");
	this->headers.push_back("Max Parent Daughter Ratio");

	this->headers.push_back("Min Partition Asymmetry");
	this->headers.push_back("Average Partition Asymmetry");
	this->headers.push_back("Max Partition Asymmetry");

	this->headers.push_back("Min Rall Power");
	this->headers.push_back("Average Rall Power");
	this->headers.push_back("Max Rall Power");

	this->headers.push_back("Min Pk");
	this->headers.push_back("Average Pk");
	this->headers.push_back("Max Pk");

	this->headers.push_back("Min Pk Classic");
	this->headers.push_back("Average Pk Classic");
	this->headers.push_back("Max Pk Classic");

	this->headers.push_back("Min Pk 2");
	this->headers.push_back("Average Pk 2");
	this->headers.push_back("Max Pk 2");

	this->headers.push_back("Average Azimuth");
	this->headers.push_back("Min Azimuth");
	this->headers.push_back("Max Azimuth");
	this->headers.push_back("Average Elevation");
	this->headers.push_back("Min Elevation");
	this->headers.push_back("Max Elevation");

	this->headers.push_back("Average Local Bifurcation Amp");
	this->headers.push_back("Min Local Bifurcation Amp");
	this->headers.push_back("Max Local Bifurcation Amp");
	this->headers.push_back("Average Local Bifurcation Tilt");
	this->headers.push_back("Min Local Bifurcation Tilt");
	this->headers.push_back("Max Local Bifurcation Tilt");	
	
	this->headers.push_back("Average Remote Bifurcation Amp");
	this->headers.push_back("Min Remote Bifurcation Amp");
	this->headers.push_back("Max Remote Bifurcation Amp");
	this->headers.push_back("Average Remote Bifurcation Tilt");
	this->headers.push_back("Min Remote Bifurcation Tilt");
	this->headers.push_back("Max Remote Bifurcation Tilt");

	this->headers.push_back("Min Leaf Level");
	this->headers.push_back("Min Leaf Path Length");
	this->headers.push_back("Max Leaf Level");
	this->headers.push_back("Max Leaf Path Length");
	this->headers.push_back("Ave Leaf Level");
	this->headers.push_back("Ave Leaf Path Length");

	this->headers.push_back( "Total Terminal Compartments");
	this->headers.push_back( "Min Terminal Compartments");
	this->headers.push_back( "Average Terminal Compartments");
	this->headers.push_back( "Max Terminal Compartments");

	this->headers.push_back("Average Diameter Threshold");
	this->headers.push_back("Diameter Threshold Max");
	this->headers.push_back("Diameter Threshold Min");
	this->headers.push_back("Average Terminal Parent Diameter");
	this->headers.push_back("Terminal Parent Diameter Max");
	this->headers.push_back("Terminal Parent Diameter Min");

	this->headers.push_back("Min Terminal Hillman Thresh"); //extend for all branches
	this->headers.push_back("Average Terminal Hillman Thresh");
	this->headers.push_back("Max Terminal Hillman Thresh");

	this->headers.push_back("Trace File");

	/*this->headers.push_back("Prediction");
	this->headers.push_back("Confidence");*/
	this->headers.push_back("Distance to Device");

	int size = this->AdditionalHeaders.size();
	for (int k = 0; k < size; k++)
	{	
		this->headers.push_back(this->AdditionalHeaders[k]);
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
	this->SetupHeaders();
	if (size == 0)
	{
		for (int i = 0; i < (int) this->Cells.size(); i ++)
		{
			this->DataTable->InsertNextRow(this->Cells[i]->DataRow());
		}
	}
	else
	{
		//std::cout << "Size of extended features " << size << std::endl;
		for (int i = 0; i < (int) this->Cells.size(); i ++)
		{
			this->DataTable->InsertNextRow(this->Cells[i]->GetExtendedDataRow(size));
		}
	}
}

vtkSmartPointer<vtkTable> CellTraceModel::ConvertDefaultValueToNull(vtkSmartPointer<vtkTable> table)
{	
	vtkSmartPointer<vtkTable> convertedTable = vtkSmartPointer<vtkTable>::New();
	convertedTable->Initialize();
	vtkSmartPointer<vtkVariantArray> column = NULL;
	for(int i=0; i < table->GetNumberOfColumns(); ++i)
    {		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( table->GetColumnName(i));
		convertedTable->AddColumn(column);
    }

	for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
	{
		vtkSmartPointer<vtkVariantArray> row = table->GetRow( i);
		row = CellTrace::ConvertDefaultValueToNull(row);
		convertedTable->InsertNextRow(row);
	}
	return convertedTable;
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
	for (int j = 0; j < (int) this->Cells.size(); j ++)
	{
		CellBoundsTable->InsertNextRow(this->Cells[j]->BoundsRow());
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
		cells.push_back(Cells.at(*IDs_iterator));
	
	return cells;
}

std::set<long int> CellTraceModel::GetSelectedIDs()
{
	std::set<long int> allSelectedIDs, nextIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		int id = (int) *it;
		bool found = false;
		unsigned int j = 0;
		while(!found&&(j<this->Cells.size()))
		{
			if (id == this->Cells.at(j)->rootID())
			{
				nextIDs = this->Cells.at(j)->TraceIDsInCell();
				std::set<long>::iterator k;
				for (k = nextIDs.begin(); k != nextIDs.end(); k++)
				{
					allSelectedIDs.insert(*k);
				}
				found = true;
			}else
			{
				j++;
			}
		}//end while !found
	}//end for selected
	return allSelectedIDs;
}
std::vector<CellTrace*> CellTraceModel::GetSelectedCells()
{
	std::vector<CellTrace*> selectedCell;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{
		int id = (int) *it;
		bool found = false;
		unsigned int j = 0;
		while(!found&&(j<this->Cells.size()))
		{
			if (id == this->Cells.at(j)->rootID())
			{
				selectedCell.push_back(this->Cells.at(j));
				found = true;
			}else
			{
				j++;
			}
		}//end while !found
	}//end for selected
	return selectedCell;
}
unsigned int CellTraceModel::getCellCount()
{
	return (unsigned int) this->Cells.size();
}
CellTrace * CellTraceModel::GetCellAt( int i)
{
	CellTrace* currentCell;
	if (i < this->Cells.size())
	{
		currentCell = this->Cells.at(i);
	}
	else
	{
		currentCell = this->Cells.back();
	}
	this->Selection->select(currentCell->rootID());
	return currentCell;
}
CellTrace * CellTraceModel::GetCellAtNoSelection( int i)
{
	CellTrace* currentCell;
	if (i < this->Cells.size())
	{
		currentCell = this->Cells.at(i);
	}
	else
	{
		currentCell = this->Cells.back();
	}
	//this->Selection->select(currentCell->rootID());
	return currentCell;
}
void CellTraceModel::WriteCellCoordsToFile(const char *fileName)
{
	fstream outputTxt;
	outputTxt.open(fileName, fstream::out);
	for (int i = 0; i < this->Cells.size(); i++)
	{
		outputTxt << this->Cells.at(i)->BasicFeatureString() << "\n";
	}
	outputTxt.close();
	std::cout << "file written\n";
}
void CellTraceModel::createCellToCellGraph()
{
	std::map< unsigned int, std::vector<float> > centroidMap;
	for(unsigned int i = 0; i < this->Cells.size(); i++)
	{
		unsigned int id = this->Cells.at(i)->rootID();
		std::vector<float> cellCoord;
		cellCoord.push_back(this->Cells.at(i)->somaX);
		cellCoord.push_back(this->Cells.at(i)->somaY);
		cellCoord.push_back(this->Cells.at(i)->somaZ);
		centroidMap[id] = cellCoord;
	}//end for cell size
	kNearestObjects* KNObj = new kNearestObjects(centroidMap);
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
	QString QnewHeader = QString(NewHeader.c_str());
	for (int k = 0; k < this->AdditionalHeaders.size(); k++)
	{
		if (QnewHeader == this->AdditionalHeaders[k])
		{
			return k;
		}
	}
	this->AdditionalHeaders.push_back(QnewHeader);
	return -1;
}
