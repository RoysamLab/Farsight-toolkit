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
	vtkSmartPointer<vtkVariantArray> GetExtendedDataRow(int CheckAddFeatures);
	vtkSmartPointer<vtkVariantArray> BoundsRow();
	std::string BasicFeatureString();
	std::set<long int> TraceIDsInCell();
	unsigned int rootID();
	TraceLine * getRootTrace();
	void setFileName(std::string newFileName);
	std::string GetFileName();
	void getSomaCoord(double xyz[]);
	void getCellBounds(double bounds[]);
	void setDistanceToROI(double newDistance, double Coord_X , double Coord_Y, double Coord_Z);
	void SetClassifcation(int predicCol, double prediction, int confCol,double confidence);
	void addNewFeature(vtkVariant nextFeature);
	std::vector<TraceLine *> getSegments();
	static vtkSmartPointer<vtkVariantArray> ConvertDefaultValueToNull(vtkSmartPointer<vtkVariantArray> row);

private:
	void clearAll();
	void MaxMin(double NewValue, double &total, double &Min, double &Max);
	void MaxMin(float NewValue, float &total, float &Min, float &Max);
	void MaxMin(int NewValue, int &total, int &Min, int &Max);
	std::vector<TraceLine*>  segments;
	int terminalBifCount;
	
public:
	int NumSegments, stems, branchPoints,terminalTips;
	int MinTerminalLevel, MaxTerminalLevel, SumTerminalLevel;

	int TotalFragmentation, MaxFragmentation, MinFragmentation;
	double TotalBurkTaper, MaxBurkTaper, MinBurkTaper;
	double TotalHillmanTaper, MaxHillmanTaper, MinHillmanTaper;
	double TotalHillmanThresh, MaxHillmanThresh, MinHillmanThresh;
	double TotalContraction, MaxContraction, MinContraction;

	double TotalDiameter, MaxDiameter, MinDiameter;
	double TotalDiameterPower, MaxDiameterPower, MinDiameterPower;

	double TotalPathLength, minPathLength, MaxPathLength;
	double TotalVolume, minSegmentVolume, maxSegmentVolume;
	double TotalEuclidianPath, MinEuclidianPath, MaxEuclidianPath;
	double TerminalPathLength, maxTerminalPathLength, minTerminalPathLength;
	int TotalTerminalSegment, MaxTerminalSegment, MinTerminalSegment;
	float somaX, somaY, somaZ, maxX, maxY, maxZ, minX, minY, minZ, skewnessX, skewnessY, skewnessZ, euclideanSkewness; 

	float sectionAreaTotal, SectionAreaMin, SectionAreaMax, surfaceAreaTotal, SurfaceAreaMax, SurfaceAreaMin;
	double somaVolume, somaSurface, SomaRadii, DiamThresholdTotal, DiamThresholdMin, DiamThresholdMax;
	double MinStemDistance, MaxStemDistance, TotalStemDistance, EstimatedSomaRadius;

	double TotalLastParentDiam, LastParentDiamMin, LastParentDiamMax;

	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double BifAmplLocal, BifAmpRemote, BifTiltLocal, BifTiltRemote;

	double daughterRatioMin, parentDaughterRatioMin, partitionAsymmetryMin, rallPowerMin, PkMin, Pk_2Min, Pk_classicMin;
	double BifAmplLocalMin, BifAmpRemoteMin, BifTiltLocalMin, BifTiltRemoteMin;

	double daughterRatioMax, parentDaughterRatioMax, partitionAsymmetryMax, rallPowerMax, PkMax, Pk_2Max, Pk_classicMax;
	double BifAmplLocalMax, BifAmpRemoteMax, BifTiltLocalMax, BifTiltRemoteMax;

	double Azimuth, AzimuthMin, AzimuthMax;
	double Elevation, ElevationMin, ElevationMax;

	double DeviceDistance;
	double prediction, confidence;
	
private:
	//double daughterRatioAverage, parentDaughterRatioAverage, partitionAsymmetryAverage, rallPowerAverage, PkAverage, Pk_2Average, Pk_classicAverage;
	//double BifAmplLocalAverage, BifAmpRemoteAverage, BifTiltLocalAverage, BifTiltRemoteAverage;
	vtkSmartPointer<vtkVariantArray> CellData;
	std::string FileName;
	std::set<long int> IDs;
	//TraceBit rootBit;
};
#endif
