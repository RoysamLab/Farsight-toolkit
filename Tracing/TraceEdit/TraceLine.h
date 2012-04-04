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
#include <set>
#include <sstream>
#include "vtkSmartPointer.h"
#include "vtkVariant.h"
#include "vtkImageData.h"

#include "TraceBit.h"

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
	TraceLine *GetParent();
	unsigned int GetParentID();
	void SetParent(TraceLine* p);
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
	double GetVolume() {return volume;}
	double GetBurkTaper() {return BurkTaper;}
	double GetHillmanTaper() {return HillmanTaper;}
	double GetHillmanThreshold() {return HillmanThreshold;}
	double GetSectionArea() {return sectionArea;}
	double GetSurfaceArea() {return surfaceArea;}
	double GetdaughterRatio() {return daughterRatio;}
	double GetparentDaughterRatio() {return parentDaughterRatio;}
	void SetParentDaughterRatio(double Ratio){parentDaughterRatio = Ratio;}
	double GetpartitionAsymmetry() {return partitionAsymmetry;}
	double GetRallPower() {return rallPower;}
	double GetPk() {return Pk;}
	double GetPk_2() {return Pk_2;}
	double GetPk_classic() {return Pk_classic;}
	double GetBifAmplLocal() {return BifAmplLocal;}
	double GetBifAmpRemote() {return BifAmpRemote;}
	void setBifTiltLocal(double NewBifTiltLocal) {BifTiltLocal = NewBifTiltLocal;}
	double GetBifTiltLocal() {return BifTiltLocal;}
	void setBifTiltRemote(double NewBifTiltRemote) {BifTiltRemote = NewBifTiltRemote;}
	double GetBifTiltRemote() {return BifTiltRemote;}
	void setBifTorqueLocal(double NewBifTorqueLocal) {BifTorqueLocal = NewBifTorqueLocal;}
	double GetBifTorqueLocal() {return BifTorqueLocal;}
	void setBifTorqueRemote(double NewBifTorqueRemote) {BifTorqueRemote = NewBifTorqueRemote;}
	double GetBifTorqueRemote() {return BifTorqueRemote;}
	
	double GetBifTiltLocalBig() {return BifTiltLocalBig;}
	void setBifTiltLocalBig(double NewBifTiltLocal)			{BifTiltLocalBig = NewBifTiltLocal;}
	double GetBifTiltLocalSmall() {return BifTiltLocalSmall;}
	void setBifTiltLocalSmall(double NewBifTiltLocal)		{BifTiltLocalSmall = NewBifTiltLocal;}
	double GetBifTiltRemoteBig() {return BifTiltRemoteBig;}
	void setBifTiltRemoteBig(double NewBifTiltRemote)		{BifTiltRemoteBig = NewBifTiltRemote;}
	double GetBifTiltRemoteSmall() {return BifTiltRemoteSmall;}
	void setBifTiltRemoteSmall(double NewBifTiltRemote)		{BifTiltRemoteSmall = NewBifTiltRemote;}
	double GetBifTorqueLocalBig() {return BifTorqueLocalBig;}
	void setBifTorqueLocalBig(double NewBifTorqueLocal)		{BifTorqueLocalBig = NewBifTorqueLocal;}
	double GetBifTorqueLocalSmall() {return BifTorqueLocalSmall;}
	void setBifTorqueLocalSmall(double NewBifTorqueLocal)	{BifTorqueLocalSmall = NewBifTorqueLocal;}
	double GetBifTorqueRemoteBig() {return BifTorqueRemoteBig;}
	void setBifTorqueRemoteBig(double NewBifTorqueRemote)	{BifTorqueRemoteBig = NewBifTorqueRemote;}
	double GetBifTorqueRemoteSmall() {return BifTorqueRemoteSmall;}
	void setBifTorqueRemoteSmall(double NewBifTorqueRemote) {BifTorqueRemoteSmall = NewBifTorqueRemote;}

	void setBifTiltLocalTwoDaughter(double NewBifTiltLocalTwo)	 {BifTiltLocalTwoDaughter = NewBifTiltLocalTwo;}
	double GetBifTiltLocalTwoDaughter() {return BifTiltLocalTwoDaughter;}
	void setBifTiltRemoteTwoDaughter(double NewBifTiltRemoteTwo) {BifTiltRemoteTwoDaughter = NewBifTiltRemoteTwo;}
	double GetBifTiltRemoteTwoDaughter() {return BifTiltRemoteTwoDaughter;}
	void setBifTorqueLocalTwoDaughter(double NewBifTorqueLocalTwo)	{BifTorqueLocalTwoDaughter = NewBifTorqueLocalTwo;}
	double GetBifTorqueLocalTwoDaughter() {return BifTorqueLocalTwoDaughter;}
	void setBifTorqueRemoteTwoDaughter(double NewBifTorqueRemoteTwo){BifTorqueRemoteTwoDaughter = NewBifTorqueRemoteTwo;}
	double GetBifTorqueRemoteTwoDaughter() {return BifTorqueRemoteTwoDaughter;}

	double GetDistanceToROI() {return DistanceToROI;}
	void SetDistanceToROI(double distance) {DistanceToROI = distance;}
	void SetDistanceToROICoord_X(double coord) {ROICoord_X = coord;}
	void SetDistanceToROICoord_Y(double coord) {ROICoord_Y = coord;}
	void SetDistanceToROICoord_Z(double coord) {ROICoord_Z = coord;}

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
	void setTraceBitIntensities(vtkSmartPointer<vtkImageData> imageData);
	void Print(std::ostream &c,int indent);

	std::vector<unsigned int> * GetMarkers();
	std::vector<TraceLine*> * GetBranchPointer();

	std::vector<double> Features; 
	//this loads in with rpi.xml files

	std::vector<vtkVariant> GetCellFeatures()
	{
		return CellFeatures;
	}
	void addCellFeature(vtkVariant Feature)
	{
		CellFeatures.push_back(Feature);
	}
	void editCellFeature(vtkVariant NewValue, int pos)
	{
		if (pos == -1)
		{
			CellFeatures.push_back(NewValue);
		}
		else if (pos < CellFeatures.size())
		{
			CellFeatures[pos] = NewValue;
		}
	}

	void setTraceColor(double newColor);
	double getTraceColor();
	void Getstats();
	bool EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
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

	//Compartment Level Features
	double GetCompartmentCurvature();
	double Euclidean(TraceBit bit1, TraceBit bit2);
private:

	double Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b);
	double Angle(TraceBit bit1, TraceBit vertex, TraceBit bit2);
	double AzimuthAngle(TraceBit vertex, TraceBit bit1);
	double ElevationAngle(TraceBit vertex, TraceBit bit1);
	double* Plane(TraceBit bit1, TraceBit vertex, TraceBit bit2);
	double PlaneAngle(double* plane1, double* plane2);
	double RallPower(double diamParent, double diamD1, double diamD2);
	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double BifAmplLocal, BifAmpRemote, BifTiltLocal, BifTiltRemote, BifTorqueLocal, BifTorqueRemote;
	double BifTiltLocalTwoDaughter, BifTiltRemoteTwoDaughter, BifTorqueLocalTwoDaughter, BifTorqueRemoteTwoDaughter;
	double BifTiltLocalBig, BifTiltLocalSmall,  BifTiltRemoteBig, BifTiltRemoteSmall;
	double BifTorqueLocalBig, BifTorqueLocalSmall,  BifTorqueRemoteBig, BifTorqueRemoteSmall;
	double traceColor, radii, sectionArea, length, volume, surfaceArea, PathLength, EuclideanD, DistToParent;
	double BitDensity, BurkTaper, HillmanTaper, HillmanThreshold;
	double BifToSomaEucDistance;
	bool actualBifurcation;
	//cell level features 
	double DistanceToROI, ROICoord_X, ROICoord_Y, ROICoord_Z;
	/*double prediction, confidence;*/


	std::string FileName;
	unsigned int m_id, root;
	int level, terminalDegree;
	std::vector<unsigned int> m_markers;
	unsigned char m_type;
	TraceLine *m_parent;
	std::vector<TraceLine* >m_branches;
	TraceBitsType m_trace_bits;

	std::vector<vtkVariant> CellFeatures;
};

#endif
