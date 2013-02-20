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
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include "vtkActor.h"

#include <vector>
#include <set>
#include <sstream>
#include <math.h>
#include "ConvexHull3D.h"

typedef itk::Image< unsigned char, 3 >   ImageType;
typedef itk::Image< float, 3> FloatImageType;

class TraceBit;
class TraceLine;
class CellTrace
{
public:
	CellTrace();
	CellTrace(std::vector<TraceLine*> Segments);
	void setTraces(std::vector<TraceLine*> Segments);
	vtkSmartPointer<vtkVariantArray> DataRow();
	vtkSmartPointer<vtkVariantArray> GetExtendedDataRow(std::vector<std::string> FeatureNames);
	vtkSmartPointer<vtkVariantArray> BoundsRow();
	std::string BasicFeatureString();
	std::set<long int> TraceIDsInCell();
	unsigned int rootID();
	TraceLine * getRootTrace();
	void setFileName(std::string newFileName);
	std::string GetFileName();
	void getSomaCoord(double xyz[]);
	void getCellBounds(double bounds[]);
	void setDistanceToROI(double newDistance, double Coord_X, double Coord_Y, double Coord_Z);
	void SetClassification(int predicCol, double prediction, int confCol,double confidence);
	void addNewFeature(std::string featureName, vtkVariant nextFeature);
	vtkVariant getFeature(std::string featureName);
	std::vector<TraceLine *> getSegments();
	std::vector<std::string> calculateConvexHull();
	vtkSmartPointer<vtkActor> GetDelaunayActor();
	vtkSmartPointer<vtkActor> GetEllipsoidActor();
	std::vector<std::string> calculateAnglesToDevice();
	std::string calculateDistanceToVessel(FloatImageType::Pointer distanceMap, ImageType::RegionType region);

	bool modified; //check if data needs to update

private:
	void clearAll();
	void MaxMin(double NewValue, double &total, double &Min, double &Max);
	void MaxMin(double NewValue, double &total, double &Min, double &Max, int &Count);
	void MaxMin(float NewValue, float &total, float &Min, float &Max);
	void MaxMin(int NewValue, int &total, int &Min, int &Max);
	void MaxMin(std::vector<double> NewValue, double &total, double &Min, double &Max, int &Count);
	std::vector<TraceLine*>  segments;
	int terminalBifCount;
	//int TerminalTriCount, notTerminalTriCount;
	//int TriCount;
	
public:
	int NumSegments, stems, branchPoints, terminalTips, actualBifurcations, branchingStem;
	int MinTerminalLevel, MaxTerminalLevel, SumTerminalLevel;

	int    FragmentationTotal, FragmentationMin, FragmentationMax;
	double BurkTaperTotal, BurkTaperMin, BurkTaperMax;
	double HillmanTaperTotal, HillmanTaperMin, HillmanTaperMax;
	double HillmanThreshTotal, HillmanThreshMin, HillmanThreshMax;
	double ContractionTotal, ContractionMin, ContractionMax;

	double DiameterTotal, DiameterMin, DiameterMax;
	double DiameterPowerTotal, DiameterPowerMin, DiameterPowerMax;

	double PathLengthTotal, PathLengthMin, PathLengthMax;
	double TotalVolume, SegmentVolumeMin, SegmentVolumeMax;
	double TotalEuclideanPath, MinEuclideanPath, MaxEuclideanPath;
	double TerminalPathLength, TerminalPathLengthMax, TerminalPathLengthMin;
	int	   TerminalSegmentTotal, TerminalSegmentMax, TerminalSegmentMin;
	float  somaX, somaY, somaZ, maxX, maxY, maxZ, minX, minY, minZ, skewnessX, skewnessY, skewnessZ, euclideanSkewness; 

	float  sectionAreaTotal, SectionAreaMin, SectionAreaMax, surfaceAreaTotal, SurfaceAreaMax, SurfaceAreaMin;
	double somaVolume, somaSurface, somaRadii, DiamThresholdTotal, DiamThresholdMin, DiamThresholdMax;
	double MinStemDistance, MaxStemDistance, TotalStemDistance, EstimatedSomaRadius;

	double TotalLastParentDiam, LastParentDiamMin, LastParentDiamMax;

	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double daughterLengthRatio, daughterLengthRatioMin, daughterLengthRatioMax;
	double BifAmplLocal, BifAmplRemote, BifTiltLocal, BifTiltRemote, BifTorqueLocal, BifTorqueRemote;
	double BifAmplLocalMin, BifAmplRemoteMin, BifTiltLocalMin, BifTiltRemoteMin, BifTorqueLocalMin, BifTorqueRemoteMin;
	double BifAmplLocalMax, BifAmplRemoteMax, BifTiltLocalMax, BifTiltRemoteMax, BifTorqueLocalMax, BifTorqueRemoteMax;
	int    BifTiltLocalCount, BifTiltRemoteCount, BifTorqueLocalCount, BifTorqueRemoteCount;

	double daughterRatioMin, parentDaughterRatioMin, partitionAsymmetryMin, rallPowerMin, PkMin, Pk_2Min, Pk_classicMin;
	double daughterRatioMax, parentDaughterRatioMax, partitionAsymmetryMax, rallPowerMax, PkMax, Pk_2Max, Pk_classicMax;

	double Azimuth, AzimuthMin, AzimuthMax;
	double Elevation, ElevationMin, ElevationMax;

	double TipToSomaEucDisTotal, TipToSomaEucDisMin, TipToSomaEucDisMax;
	double BranchPtToSomaEucDisTotal, BranchPtToSomaEucDisMin, BranchPtToSomaEucDisMax;
	double totalTipX, totalTipY, totalTipZ;
	double tipMagnitude, tipAzimuth, tipElevation;

	double DeviceDistance, DeviceAzimuth, DeviceElevation, cellDirectiontoDevice;
	double prediction, confidence;

	double convexHullMagnitude, convexHullAzimuth, convexHullElevation, convexHullArea, convexHullVol;
	
private:
	vtkSmartPointer<vtkVariantArray> CellData;
	std::string FileName;
	std::set<long int> IDs;

	bool delaunayCreated;
	vtkSmartPointer<vtkActor> delaunayActor;
	vtkSmartPointer<vtkActor> ellipsoidActor;

	std::vector<TraceBit> tips;
};
#endif
