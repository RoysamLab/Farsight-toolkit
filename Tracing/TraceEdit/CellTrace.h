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
#include <set>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <sstream>
#include <math.h>
class TraceBit;
class TraceLine;
class CellTrace
{
public:
	CellTrace();
	CellTrace(std::vector<TraceLine*> Segments);
	void setTraces(std::vector<TraceLine*> Segments);
	vtkSmartPointer<vtkVariantArray> DataRow();
	vtkSmartPointer<vtkVariantArray> BoundsRow();
	std::string BasicFeatureString();
	std::set<long int> TraceIDsInCell();
	int rootID();
	TraceLine * getRootTrace();
	void setFileName(std::string newFileName);
	std::string GetFileName();
	void getSomaCoord(double xyz[]);
	void getCellBounds(double bounds[]);
private:
	void clearAll();
	std::vector<TraceLine*>  segments;
	int NumSegments, stems, branchPoints,terminalTips, MinTerminalLevel, MaxTerminalLevel, SumTerminalLevel;
	double TotalPathLength, TotalVolume, TotalEuclidianPath, TerminalPathLength;
	double maxTerminalPathLength, minTerminalPathLength;
	float somaX, somaY, somaZ, maxX, maxY, maxZ, minX, minY, minZ; 
	float sectionAreaTotal, SectionAreaMin, SectionAreaMax, surfaceAreaTotal, SurfaceAreaMax, SurfaceAreaMin;
	double somaVolume, somaSurface, SomaRadii, DiamThresholdTotal, DiamThresholdMin, DiamThresholdMax;
	double TotalLastParentDiam, LastParentDiamMin, LastParentDiamMax;
	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double BifAmplLocal, BifAmpRemote, BifTiltLocal, BifTiltRemote;
	std::string FileName;
	std::set<long int> IDs;
	//TraceBit rootBit;
};
#endif