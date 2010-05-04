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
	//this->rootBit = this->segments[0]->GetTraceBitsPointer()->front();
	for(i = 0; i < this->segments.size(); i++)
	{
		this->TotalPathLength += this->segments[i]->GetLength();
		this->TotalVolume += this->segments[i]->GetVolume();
		this->TotalEuclidian += this->segments[i]->GetEuclidianLength();
		int tempLevel = this->segments[i]->GetLevel();
		if (this->segments[i]->isLeaf())
		{
			this->terminalTips++;
			this->TerminalPathLength += this->segments[i]->GetPathLength();
			if(tempLevel > this->MaxTerminalLevel)
			{
				this->MaxTerminalLevel = tempLevel;
			}
			if(tempLevel < this->MinTerminalLevel)
			{
				this->MinTerminalLevel = tempLevel;
			}
		}//end if leaf
	}//end for segment size
}
void CellTrace::clearAll()
{
	this->segments.clear();
	this->NumSegments = 0;
	this->stems = 0;
	this->terminalTips = 0;
	this->MinTerminalLevel = 0;
	this->MaxTerminalLevel = 0;
	this->TotalEuclidian = 0;
	this->TotalPathLength = 0;
	this->TotalVolume = 0;
	this->TerminalPathLength = 0;
}
vtkSmartPointer<vtkVariantArray> CellTrace::DataRow()
{
	vtkSmartPointer<vtkVariantArray> CellData = vtkSmartPointer<vtkVariantArray>::New();
	CellData->InsertNextValue(this->segments[0]->GetId());
	CellData->InsertNextValue(this->NumSegments);
	CellData->InsertNextValue(this->stems);
	CellData->InsertNextValue(this->terminalTips);
	CellData->InsertNextValue(this->MinTerminalLevel);
	CellData->InsertNextValue(this->MaxTerminalLevel);
	CellData->InsertNextValue(this->TotalEuclidian);
	CellData->InsertNextValue(this->TotalPathLength);
	CellData->InsertNextValue(this->TotalVolume);
	CellData->InsertNextValue(this->TerminalPathLength/this->NumSegments);
	return CellData;
}