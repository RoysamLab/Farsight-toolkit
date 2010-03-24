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
	double GetPathLength()
		{return PathLength;}
	void calculateVol();
	double GetLength() {return length;}
	double GetEuclidianLength();
	double GetFragmentationSmoothness();
	double GetRadii(){return radii;}
	double GetVolume() {return volume;}
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
	TraceBit removeLastBit();
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
	bool Orient(TraceLine * Trunk);
	bool Orient(TraceBit bit);
	std::string stats();	

private:

	double Euclidian(TraceBit bit1, TraceBit bit2);
	double Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b);

	double traceColor, radii, length, volume, PathLength, EuclidianD;
	int m_id, root, level;
	std::vector<unsigned int> m_markers;
	unsigned char m_type;
	TraceLine *m_parent;
	std::vector<TraceLine* >m_branches;
	TraceBitsType m_trace_bits;
};

#endif
