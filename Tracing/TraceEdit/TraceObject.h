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

#ifndef __TRACEOBJECT_H
#define __TRACEOBJECT_H
#include <vector>
#include <vtksys/hash_map.hxx> /* Platform independent hashmap */
#include "vtkSmartPointer.h"

class TraceBit;
class TraceLine;
class TraceGap;
class branchPT;

class vtkPoints;
class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;

/* A Trace object is a list of root TraceLines*/
class TraceObject
{
public:
	TraceObject();
	TraceObject(const TraceObject &T);
	~TraceObject();
	double getSmallLineColor()
	{
		return this->smallLineColor;
	};
	double getMergeLineColor()
	{
		return this->mergeLineColor;
	};
	void setSmallLineColor(double set)
	{
		this->smallLineColor=set;
	};
	void setMergeLineColor(double set)
	{
		this->mergeLineColor=set;
	};
	void setLUT(int num);
	double getTraceLUT(unsigned char type);
//	I/O functions
	bool ReadFromSWCFile(char * filename);
	bool ReadFromRPIXMLFile(char * filename);
	bool ReadFromFeatureTracksFile(char *filename, int type_offset);
	bool ReadFromFeatureTracksFileForKymograph(char *filename,int type_offset);
	bool WriteToSWCFile(const char * filename);
	void WriteToVTKFile(const char * filename);
//	operators
	int getNewLineId();
	void splitTrace(int selectedCellId);
	void ReverseSegment(TraceLine*);
	void RemoveTraceLine(TraceLine*);
	void FixPointMarkers(TraceLine* tline);
	void mergeTraces(unsigned long long int eMarker, unsigned long long int sMarker);
	void CreatePolyDataRecursive(TraceLine* , vtkSmartPointer<vtkFloatArray> , vtkSmartPointer<vtkPoints> ,vtkSmartPointer<vtkCellArray>);
	void FindMinLines(int smallSize);
	int createGapLists(std::vector<TraceLine*> traceList);
//	public data
	vtkSmartPointer<vtkPolyData> GetVTKPolyData();
	vtkSmartPointer<vtkPolyData> generateBranchIllustrator();
	void Print(std::ostream &c);

	std::vector<TraceLine*>* GetTraceLinesPointer();
	std::vector<TraceLine*> GetTraceLines();
	std::vector<TraceBit> CollectTraceBits();
	std::vector<int> SmallLines;
	std::vector<TraceGap*> Gaps;
	std::vector<std::string> FeatureHeaders;

	double gapTol;
	int gapMax;

	vtksys::hash_map< unsigned int, unsigned long long int > hashp;
	vtksys::hash_map< unsigned int, unsigned long long int > hashc;
	
private:
	std::vector<TraceLine*> trace_lines;
	vtkSmartPointer<vtkPolyData> PolyTraces;
	double smallLineColor, mergeLineColor;	
  void CollectTraceBitsRecursive(std::vector<TraceBit> &vec,TraceLine *l);
  void CollectIdsRecursive(std::vector<int> &ids, TraceLine* tline);
  void LinearTraceLinesRecursive(std::vector<TraceLine*> &allLine, TraceLine* tline);
  void CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells,TraceLine *tline);
  void CollectSegmentMidPointsRecursive(vtkSmartPointer<vtkPoints>p, vtkSmartPointer<vtkCellArray> cells, vtkSmartPointer<vtkFloatArray> da,TraceLine* tline);
};

#endif

