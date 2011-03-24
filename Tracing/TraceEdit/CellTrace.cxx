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

#include "TraceBit.h"
#include "TraceLine.h"
#include "CellTrace.h"
CellTrace::CellTrace()
{
	this->clearAll();
}
CellTrace::CellTrace(std::vector<TraceLine*> Segments)
{
	this->clearAll();
	this->setTraces(Segments);
}
void CellTrace::setTraces(std::vector<TraceLine*> Segments)
{
	this->segments = Segments;
	unsigned int i = 0;
	this->NumSegments = (int) this->segments.size();
	this->stems = (int) this->segments[0]->GetBranchPointer()->size();
	TraceBit rootBit = this->segments[0]->GetTraceBitsPointer()->front();
	//set the soma point as the intital bounds
	this->somaX = rootBit.x;
	this->somaY = rootBit.y;
	this->somaZ = rootBit.z;
	this->maxX = rootBit.x;
	this->maxY = rootBit.y;
	this->maxZ = rootBit.z;
	this->minX = rootBit.x;
	this->minY = rootBit.y;
	this->minZ = rootBit.z;
	//std::cout << this->segments[0]->GetFileName()<< std::endl;
	//this->FileName = this->segments[0]->GetFileName();
	for(i = 0; i < this->segments.size(); i++)
	{
		this->IDs.insert(this->segments[i]->GetId());
		this->TotalPathLength += this->segments[i]->GetLength();
		this->TotalVolume += this->segments[i]->GetVolume();
		this->TotalEuclidianPath += this->segments[i]->GetEuclidianLength();
		int tempLevel = this->segments[i]->GetLevel();
		if (this->segments[i]->isLeaf())
		{
			TraceBit leafBit = this->segments[i]->GetTraceBitsPointer()->back();
			float lx = leafBit.x;
			float ly = leafBit.y;
			float lz = leafBit.z;

			if (lx > this->maxX)
			{
				this->maxX = lx;
			}else if ( lx < this->minX)
			{
				this->minX = lx;
			}

			if (ly > this->maxY)
			{
				this->maxY = ly;
			}else if ( ly < this->minY)
			{
				this->minY = ly;
			}

			if (lz > this->maxZ)
			{
				this->maxZ = lz;
			}else if ( lz < this->minZ)
			{
				this->minZ = lz;
			}

			this->terminalTips++;
			this->SumTerminalLevel += tempLevel;
			this->TerminalPathLength += this->segments[i]->GetPathLength();
			if(tempLevel > this->MaxTerminalLevel)
			{
				this->MaxTerminalLevel = tempLevel;
				this->maxTerminalPathLength = this->segments[i]->GetPathLength();
			}
			if(tempLevel < this->MinTerminalLevel)
			{
				this->MinTerminalLevel = tempLevel;
				this->minTerminalPathLength = this->segments[i]->GetPathLength();
			}
		}//end if leaf
		else if(!this->segments[i]->isRoot())
		{
			this->branchPoints++;
		}
	}//end for segment size
}
void CellTrace::setFileName(std::string newFileName)
{
	this->FileName = newFileName;
}
std::string CellTrace::GetFileName()
{
	std::string outputFileName;
	if (this->FileName != "file")
	{
		outputFileName = this->FileName;
	}
	else
	{
		std::stringstream newCellName; 
		newCellName<<"cell" << this->segments[0]->GetId();
		outputFileName = newCellName.str();
	}
	return outputFileName;
}
void CellTrace::getSomaCoord(double xyz[])
{
	xyz[0] = this->somaX;
	xyz[1] = this->somaY;
	xyz[2] = this->somaZ;
}
void CellTrace::getCellBounds(double bounds[])
{//min then max of x, y, z
	bounds[0] = this->minX;
	bounds[1] = this->maxX;
	bounds[2] = this->minY;
	bounds[3] = this->maxY;
	bounds[4] = this->minZ;
	bounds[5] = this->maxZ;
}
void CellTrace::clearAll()
{
	this->segments.clear();
	this->NumSegments = 0;
	this->stems = 0;
	this->branchPoints = 0;
	this->terminalTips = 0;
	this->MinTerminalLevel = 100; //something large for initial value
	this->MaxTerminalLevel = 0;
	this->SumTerminalLevel = 0;
	this->TotalEuclidianPath = 0;
	this->TotalPathLength = 0;
	this->TotalVolume = 0;
	this->TerminalPathLength = 0;
	this->FileName = "file";
	this->somaX = 0;
	this->somaY = 0;
	this->somaZ = 0;
	this->maxX = 0;
	this->maxY = 0;
	this->maxZ = 0;
	this->minX = 0;
	this->minY = 0;
	this->minZ = 0;
}
vtkSmartPointer<vtkVariantArray> CellTrace::DataRow()
{
	vtkSmartPointer<vtkVariantArray> CellData = vtkSmartPointer<vtkVariantArray>::New();
	CellData->InsertNextValue(this->segments[0]->GetId());
	CellData->InsertNextValue(this->NumSegments);
	CellData->InsertNextValue(this->stems);
	CellData->InsertNextValue(this->branchPoints);
	CellData->InsertNextValue(this->terminalTips);
	CellData->InsertNextValue(this->MinTerminalLevel);
	CellData->InsertNextValue(this->minTerminalPathLength);
	CellData->InsertNextValue(this->MaxTerminalLevel);
	CellData->InsertNextValue(this->maxTerminalPathLength);
	CellData->InsertNextValue(this->SumTerminalLevel /this->terminalTips);//average terminal level
	CellData->InsertNextValue(this->TerminalPathLength/this->terminalTips);//now average path to end
	CellData->InsertNextValue(this->TotalEuclidianPath);
	CellData->InsertNextValue(this->TotalPathLength);
	CellData->InsertNextValue(this->TotalPathLength/this->NumSegments);//average segment length
	CellData->InsertNextValue(this->TotalVolume);
	CellData->InsertNextValue(this->somaX);
	CellData->InsertNextValue(this->somaY);
	CellData->InsertNextValue(this->somaZ);
	CellData->InsertNextValue(this->GetFileName().c_str());
	//std::cout << this->FileName << std::endl;
	return CellData;
}
std::set<long int> CellTrace::TraceIDsInCell()
{
	return this->IDs;
}
int CellTrace::rootID()
{
	return this->segments[0]->GetId();
}
