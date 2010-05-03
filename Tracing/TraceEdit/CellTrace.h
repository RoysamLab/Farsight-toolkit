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
#ifndef __CELLTRACE_H
#define __CELLTRACE_H
#include <vector>
#include "vtkSmartPointer.h"
class TraceBit;
class TraceLine;
class CellTrace
{
public:
	CellTrace();
	CellTrace(std::vector<TraceLine*> Segments);
	void setTraces(std::vector<TraceLine*> Segments);
private:
	void clearAll();
	std::vector<TraceLine*>  segments;
	int NumSegments, stems, terminalTips, MinTerminalLevel, MaxTerminalLevel;
	double TotalPathLength, TotalVolume, TotalEuclidian, TerminalPathLength; 
	//TraceBit rootBit;
};
#endif