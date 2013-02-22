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

#ifndef __TRACELINE_H
#define __TRACELINE_H

#define PI 3.14159265

#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include "vtkSmartPointer.h"
#include "vtkVariant.h"
#include "vtkImageData.h"

#include "TraceBit.h"
#include "StructuredObject.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIterator.h"

typedef itk::Image< ImageActorPixelType, Dimension >   ImageType;
typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > VolumeOfInterestFilterType;

/**
 * A TraceLine is a sequence of TraceBits that has pointers to two other
 * TraceLines
 **/
class TraceLine
{
public:
  typedef std::list<TraceBit> TraceBitsType;
	TraceLine();
	TraceLine(const TraceLine &t);
  ~TraceLine();
	bool modified;

	TraceLine *GetParent(int i);
	std::vector<TraceLine*> GetParents();
	int ParentSize();
	unsigned int GetParentID(int i);
	bool isParentLess();
	bool isMarked();
	void MarkLine();
	void UnmarkLine();
	void SetParent(TraceLine* p);
	void RemoveParents();
	void RemoveParent(int i);
	int GetParentNumber(TraceLine* p);

	int GetRootID();
	int GetLevel();
	int GetTerminalDegree() {return terminalDegree;}
	void setTerminalDegree(int degree);
	double GetPathLength()
		{return PathLength;}
	void calculateVol();
	void calculateBifFeatures();
	double GetLength() {return length;}
	double GetEuclideanLength();
	double GetBitDensity();
	double GetDistToParent();
	double GetFragmentationSmoothness();
	double GetRadii(){return radii;}
	double GetSomaVolume() {return somaVolume;}
	double GetSomaSurfaceArea() {return somaSurfaceArea;}
	double GetVolume() {return volume;}
	double GetBurkTaper() {return BurkTaper;}
	double GetHillmanTaper() {return HillmanTaper;}
	double GetHillmanThreshold() {return HillmanThreshold;}
	double GetSectionArea() {return sectionArea;}
	double GetSurfaceArea() {return surfaceArea;}
	double GetdaughterRatio() {return daughterRatio;}
	double GetparentDaughterRatio() {return parentDaughterRatio;}
	void SetParentDaughterRatio(double Ratio){parentDaughterRatio = Ratio;}
	double GetdaughterLengthRatio() {return daughterLengthRatio;}
	double GetpartitionAsymmetry() {return partitionAsymmetry;}
	double GetRallPower() {return rallPower;}
	double GetPk() {return Pk;}
	double GetPk_2() {return Pk_2;}
	double GetPk_classic() {return Pk_classic;}
	double GetBifAmplLocal() {return BifAmplLocal;}
	double GetBifAmplRemote() {return BifAmplRemote;}
	void setBifTiltLocal(std::vector<double> NewBifTiltLocal) {BifTiltLocal = NewBifTiltLocal;}
	std::vector<double> GetBifTiltLocal() {return BifTiltLocal;}
	void setBifTiltRemote(std::vector<double> NewBifTiltRemote) {BifTiltRemote = NewBifTiltRemote;}
	std::vector<double> GetBifTiltRemote() {return BifTiltRemote;}
	void setBifTiltLocalAvg(double NewBifTiltLocal) {BifTiltLocalAvg = NewBifTiltLocal;}
	double GetBifTiltLocalAvg() {return BifTiltLocalAvg;}
	void setBifTiltRemoteAvg(double NewBifTiltRemote) {BifTiltRemoteAvg = NewBifTiltRemote;}
	double GetBifTiltRemoteAvg() {return BifTiltRemoteAvg;}

	void setBifTorqueLocal(std::vector<double> NewBifTorqueLocal) {planeAngleLocal = NewBifTorqueLocal;}
	std::vector<double> GetBifTorqueLocal() {return planeAngleLocal;}
	void setBifTorqueRemote(std::vector<double> NewBifTorqueRemote) {planeAngleRemote = NewBifTorqueRemote;}
	std::vector<double> GetBifTorqueRemote() {return planeAngleRemote;}
	
	void setBifTorqueLocalAvg(double NewBifTorqueLocal) {BifTorqueLocalAvg = NewBifTorqueLocal;}
	double GetBifTorqueLocalAvg() {return BifTorqueLocalAvg;}
	void setBifTorqueRemoteAvg(double NewBifTorqueRemote) {BifTorqueRemoteAvg = NewBifTorqueRemote;}
	double GetBifTorqueRemoteAvg() {return BifTorqueRemoteAvg;}

	double GetDistanceToROI() {return DistanceToROI;}
	double GetAzimuthToROI() {return ROIAzimuth;}
	double GetElevationToROI() {return ROIElevation;}
	double GetTipToROI() {return somaROIAngle;}
	void SetDistanceToROI(double distance) {DistanceToROI = distance;}
	void SetDistanceToROICoord_X(double coord) {ROICoord_X = coord;}
	void SetDistanceToROICoord_Y(double coord) {ROICoord_Y = coord;}
	void SetDistanceToROICoord_Z(double coord) {ROICoord_Z = coord;}
	void CalculateDirectionToROI(TraceBit tipsPt);

	/*void SetClassification(double newPrediction, double newConfidence)
	{
		prediction = newPrediction;
		confidence = newConfidence;
	}
	double getPrediction()
	{
		return prediction;
	}
	double getConfidence()
	{
		return confidence;
	}*/

	void setRoot(int RootID, int traceLevel, double parentPath);
	void setRoot(int RootID);
	void AddBranch(TraceLine* b);
	TraceLine *GetBranch1();
	void SetBranch1(TraceLine* b0);
	TraceLine *GetBranch2();
	void SetBranch2(TraceLine* b1);
	bool isLeaf();
	bool isRoot();
	bool isFree();
	bool isBranch();
	bool isActualBifurcation();
	void setActualBifurcation(bool bifurcate);
	unsigned char GetType();
	void SetType(unsigned char t) ;
	void AddTraceBit(TraceBit tbit);
	void ExtendTrace(TraceBit tbit);
	bool removeLeadingBit();
	TraceBit removeLastBit();
	TraceBit removeFirstBit();
	TraceBit GetBitXFromEnd(int x);
	TraceBit GetBitXFromBegin(int x);
	TraceBitsType::iterator GetTraceBitIteratorBegin();
	TraceBitsType::iterator GetTraceBitIteratorEnd();
	TraceBitsType * GetTraceBitsPointer();
	void SetId(unsigned int lid);
	unsigned int GetId();
	int GetSize();
	void setTraceBitIntensities(vtkSmartPointer<vtkImageData> imageData, std::string ImageName);
	void setTraceBitWeightedIntensities(ImageType::Pointer input_image, std::string ImageName);

	void Print(std::ostream &c,int indent);

	std::vector<unsigned int> * GetMarkers();
	std::vector<TraceLine*> * GetBranchPointer();

	std::vector<double> Features; 
	//this loads in with rpi.xml files
	vtkVariant GetTraceFeature(std::string FeatureName);
	void SetTraceFeature(std::string FeatureName,vtkVariant FeatureValue);
	vtkVariant GetCellFeature(std::string FeatureName);
	void SetCellFeature(std::string FeatureName,vtkVariant FeatureValue);

	void setTraceColor(double newColor);
	double getTraceColor();
	void Getstats();
	bool EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
                 double &maxdist, double &angle);
	void EndPtDistVessel(TraceLine *Trace2, TraceBit &dir1, TraceBit &dir2, double &dist,
                 double &maxdist, double &angle);
	double CalculatePk(double Dp, double Da, double Db, double n);
	bool Orient(TraceLine * Trunk);
	bool Orient(TraceBit bit);
	void SetFileName(char * newFileName);
	const char * GetFileName();
	void getEndPtBounds(double bounds[]);
	std::string stats();	
	std::string RootCoord();	
	std::string statHeaders();

	double GetAzimuth();
	double GetElevation();
	double GetAngle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b);

	//Compartment Level Features
	double GetCompartmentCurvature();
	double Euclidean(TraceBit bit1, TraceBit bit2);
private:
	std::set<long int> bitIDs;

	double Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b);
	double Angle(TraceBit bit1, TraceBit vertex, TraceBit bit2);
	double AzimuthAngle(TraceBit vertex, TraceBit bit1);
	double ElevationAngle(TraceBit vertex, TraceBit bit1);
	void   Plane(TraceBit bit1, TraceBit vertex, TraceBit bit2, double vector[]);
	double PlaneAngle(double* plane1, double* plane2);
	double RallPower(double diamParent, double diamD1, double diamD2);
	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double BifAmplLocal, BifAmplRemote, BifTorqueLocal, BifTorqueRemote;
	std::vector<double> BifTiltLocal, BifTiltRemote;
	std::vector<double> planeAngleLocal, planeAngleRemote;
	double BifTiltLocalAvg, BifTiltRemoteAvg, BifTorqueLocalAvg, BifTorqueRemoteAvg;
	double traceColor, radii, sectionArea, length, volume, surfaceArea, PathLength, EuclideanD, DistToParent;
	double somaVolume, somaSurfaceArea;
	double BitDensity, BurkTaper, HillmanTaper, HillmanThreshold;
	double BifToSomaEucDistance;
	double daughterLengthRatio;
	bool actualBifurcation;
	//cell level features 
	double DistanceToROI, ROICoord_X, ROICoord_Y, ROICoord_Z, ROIAzimuth, ROIElevation, somaROIAngle;
	/*double prediction, confidence;*/


	std::string FileName;
	unsigned int m_id, root;
	int level, terminalDegree;
	std::vector<unsigned int> m_markers;
	unsigned char m_type;
	
	std::vector<TraceLine* >m_parent; 
	bool marked;

	std::vector<TraceLine* >m_branches;
	TraceBitsType m_trace_bits;

	std::map<std::string, vtkVariant> CellFeatures;
	std::map<std::string, vtkVariant> TraceFeatures;
};

#endif
