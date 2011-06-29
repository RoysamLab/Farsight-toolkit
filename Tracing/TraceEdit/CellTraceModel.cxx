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
	this->Cells.clear();
}
CellTraceModel::CellTraceModel(std::vector<CellTrace*> Cells)
{	
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->setCells(Cells);
}
CellTraceModel::~CellTraceModel()
{	
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
	this->headers.push_back("Segments");
	this->headers.push_back("Stems");
	this->headers.push_back("Branch Pt");
	this->headers.push_back("Leaf Nodes");

	this->headers.push_back("Min Leaf Level");
	this->headers.push_back("Min Leaf Path Length");
	this->headers.push_back("Max Leaf Level");
	this->headers.push_back("Max Leaf Path Length");
	this->headers.push_back("Ave Leaf Level");
	this->headers.push_back("Ave Leaf Path Length");
	this->headers.push_back("Total Euclidian Length");
	this->headers.push_back("Average Segment Euclidian Length");
	this->headers.push_back("Total Path Length");
	this->headers.push_back("Average Segment Path Length");

	this->headers.push_back("Total Fragmentation");
	this->headers.push_back("Min Fragmentation");
	this->headers.push_back("Average Fragmentation");
	this->headers.push_back("Max Fragmentation");

	this->headers.push_back("Min Burk Taper");
	this->headers.push_back("Average Burk Taper");
	this->headers.push_back("Max Burk Taper");

	this->headers.push_back("Min Hillman Taper");
	this->headers.push_back("Average Hillman Taper");
	this->headers.push_back("Max Hillman Taper");

	this->headers.push_back("Min Contraction");
	this->headers.push_back("Average Contraction");
	this->headers.push_back("Max Contraction");

	this->headers.push_back("Min Diameter");
	this->headers.push_back("Average Diameter");
	this->headers.push_back("Max Diameter");

	this->headers.push_back("Min Diameter Power");
	this->headers.push_back("Average Diameter Power");
	this->headers.push_back("Max Diameter Power");

	this->headers.push_back("Average Diameter Threshold");
	this->headers.push_back("Diameter Threshold Max");
	this->headers.push_back("Diameter Threshold Min");
	this->headers.push_back("Average Parent Diameter");
	this->headers.push_back("Parent Diameter Max");
	this->headers.push_back("Parent Diameter Min");

	this->headers.push_back( "Total Terminal Compartments");
	this->headers.push_back( "Min Terminal Compartments");
	this->headers.push_back( "Average Terminal Compartments");
	this->headers.push_back( "Max Terminal Compartments");

	this->headers.push_back("Total Volume");
	this->headers.push_back("Min Segment Volume");
	this->headers.push_back("Max Segment Volume");
	this->headers.push_back("Average Segment Volume");
	this->headers.push_back("Total Surface Area");
	this->headers.push_back("Max Segment Surface Area");
	this->headers.push_back("Average Segment Surface Area");
	this->headers.push_back("Min Segment Surface Area");

	this->headers.push_back("Average Local Bifurcation Amp");
	this->headers.push_back("Min Local Bifurcation Amp");
	this->headers.push_back("Max Local Bifurcation Amp");
	/*this->headers.push_back("Average Local Bifurcation Tilt");
	this->headers.push_back("Min Local Bifurcation Tilt");
	this->headers.push_back("Max Local Bifurcation Tilt");*/	
	
	this->headers.push_back("Average Remote Bifurcation Amp");
	this->headers.push_back("Min Remote Bifurcation Amp");
	this->headers.push_back("Max Remote Bifurcation Amp");
	/*this->headers.push_back("Average Remote Bifurcation Tilt");
	this->headers.push_back("Min Remote Bifurcation Tilt");
	this->headers.push_back("Max Remote Bifurcation Tilt");*/

	this->headers.push_back("Min Daughter Ratio");
	this->headers.push_back("Average Daughter Ratio");
	this->headers.push_back("Max Daughter Ratio");

	this->headers.push_back("Min Parent Daughter Ratio");
	this->headers.push_back("Average Parent Daughter Ratio");
	this->headers.push_back("Max Parent Daughter Ratio");

	this->headers.push_back("Min Partition Asymmetry");
	this->headers.push_back("Average Partition Asymmetry");
	this->headers.push_back("Max Partition Asymmetry");

	this->headers.push_back("Min Hillman Thresh");
	this->headers.push_back("Average Hillman Thresh");
	this->headers.push_back("Max Hillman Thresh");

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

	this->headers.push_back("Width X");
	this->headers.push_back("Height Y");
	this->headers.push_back("Depth Z");

	this->headers.push_back("Soma X");
	this->headers.push_back("Soma Y");
	this->headers.push_back("Soma Z");
	this->headers.push_back("Soma Radii");
	this->headers.push_back("Soma Volume");
	this->headers.push_back("Soma Surface Area");
	this->headers.push_back("Trace File");
	int numHeaders = (int)this->headers.size();
	std::cout<<numHeaders << "\t features computed\n";
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
	this->DataTable->Initialize();	
	this->Selection->clear();
	this->SetupHeaders();
	for (int i = 0; i < (int) this->Cells.size(); i ++)
	{
		this->DataTable->InsertNextRow(this->Cells[i]->DataRow());
	}
}

vtkSmartPointer<vtkTable> CellTraceModel::getDataTable()
{
	return this->DataTable;
}

vtkSmartPointer<vtkTable> CellTraceModel::getCellBoundsTable()
{
	vtkSmartPointer<vtkTable> CellBoundsTable = vtkSmartPointer<vtkTable>::New();
	CellBoundsTable->Initialize();
	std::vector<QString> BoundsHeaders;
	BoundsHeaders.clear();
	BoundsHeaders.push_back("Cell ID");
	BoundsHeaders.push_back("Soma X");
	BoundsHeaders.push_back("Soma Y");
	BoundsHeaders.push_back("Soma Z");
	BoundsHeaders.push_back("Min X");
	BoundsHeaders.push_back("Max X");
	BoundsHeaders.push_back("Min Y");
	BoundsHeaders.push_back("Max Y");
	BoundsHeaders.push_back("Min Z");
	BoundsHeaders.push_back("Max Z");
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
void CellTraceModel::SelectByRootTrace(std::vector<TraceLine*> roots)
{
	this->Selection->clear();
	std::set<long int> ID;
	unsigned int i = 0;
	for (i = 0; i < roots.size(); i++)
	{
		ID.insert((long)roots.at(i)->GetId());
		//std::cout<< "root" << roots.at(i)->GetId()<< " at" <<i<< std::endl;
	}
	this->Selection->select(ID);
}
std::set<long int> CellTraceModel::GetSelecectedIDs()
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
std::vector<CellTrace*> CellTraceModel::GetSelecectedCells()
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
int CellTraceModel::getCellCount()
{
	return (int) this->Cells.size();
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