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

#include <list>
#include <vector>
#include <set>
#include <sstream>
#include "vtkSmartPointer.h"
#include "vtkImageData.h"

class TraceBit;

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
	int GetParentID();
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
	double GetEuclidianLength();
	double GetBitDensity();
	double GetDistToParent();
	double GetFragmentationSmoothness();
	double GetRadii(){return radii;}
	double GetVolume() {return volume;}
	double GetBurkTaper() {return BurkTaper;}
	double GetHillmanTaper() {return HillmanTaper;}
	double GetSectionArea() {return sectionArea;}
	double GetSurfaceArea() {return surfaceArea;}
	double GetdaughterRatio() {return daughterRatio;}
	double GetparentDaughterRatio() {return parentDaughterRatio;}
	double GetpartitionAsymmetry() {return partitionAsymmetry;}
	double GetrallPower() {return rallPower;}
	double GetPk() {return Pk;}
	double GetPk_2() {return Pk_2;}
	double GetPk_classic() {return Pk_classic;}
	double GetBifAmplLocal() {return BifAmplLocal;}
	double GetBifAmpRemote() {return BifAmpRemote;}
	double GetBifTiltLocal() {return BifTiltLocal;}
	double GetBifTiltRemote() {return BifTiltRemote;}
	void setRoot(int RootID, int traceLevel, double parentPath);
	void AddBranch(TraceLine* b);
	TraceLine *GetBranch1();
	void SetBranch1(TraceLine* b0);
	TraceLine *GetBranch2();
	void SetBranch2(TraceLine* b1);
	bool isLeaf();
	bool isRoot();
	bool isFree();
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
	void SetId(int lid);
	int GetId();
	int GetSize();
	void setTraceBitIntensities(vtkSmartPointer<vtkImageData> imageData);
	void Print(std::ostream &c,int indent);

	std::vector<unsigned int> * GetMarkers();
	std::vector<TraceLine*> * GetBranchPointer();
	std::vector<double> Features;
	void setTraceColor(double newColor);
	double getTraceColor();
	void Getstats();
	bool EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
                 double &maxdist, double &angle);
	double CalculatePk(double Dp, double Da, double Db, double n);
	bool Orient(TraceLine * Trunk);
	bool Orient(TraceBit bit);
	void SetFileName(char * newFileName);
	char * GetFileName();
	void getEndPtBounds(double bounds[]);
	std::string stats();	
	std::string statHeaders();

private:

	double Euclidian(TraceBit bit1, TraceBit bit2);
	double Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b);
	double Angle(TraceBit bit1, TraceBit vertex, TraceBit bit2);
	double daughterRatio, parentDaughterRatio, partitionAsymmetry, rallPower, Pk, Pk_2, Pk_classic;
	double BifAmplLocal, BifAmpRemote, BifTiltLocal, BifTiltRemote;
	double traceColor, radii, sectionArea, length, volume, surfaceArea, PathLength, EuclidianD, DistToParent;
	double BitDensity, BurkTaper, HillmanTaper;
	char * FileName; 
	int m_id, root, level, terminalDegree;
	std::vector<unsigned int> m_markers;
	unsigned char m_type;
	TraceLine *m_parent;
	std::vector<TraceLine* >m_branches;
	TraceBitsType m_trace_bits;
};

#endif
