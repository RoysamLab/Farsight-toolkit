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
#include <set>
#include <vtksys/hash_map.hxx> /* Platform independent hashmap */
#include "vtkSmartPointer.h"
#include "vtkImageData.h"

class TraceBit;
class TraceLine;
class TraceGap;
class branchPT;

class vtkPoints;
class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;

typedef struct {
	TraceBit* theBit;
	TraceLine* theBitLine;
	TraceLine* theAnitBitLine;
	TraceLine* theAnitBitLineAdjusted;
	double slope[3];
	//The plane equation is ((slope[0]*x + slope[1]*y)/slope[3])*t def
	bool Front;
} EndPointInfo;
typedef struct {
	TraceBit* arrayOBits[2];
	double slope[3];
	double intersect[3];
	//The plane equation is ((slope[0]*x + slope[1]*y)/slope[3])*t def
	bool Front;
} adjPointLines;
const int timeSize = 20; //this should be user en
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
	void leastSquaresThreeDependentVar(EndPointInfo* endPoint, std::vector<TraceBit*>*  ExtrapolationLines);
	void Calculatebranches();
	std::vector<EndPointInfo>* findingBranchingPoints();
	std::vector<EndPointInfo>* findingBranchingPointintercepts(std::vector<EndPointInfo>* PossibleBranches);
	void setLUT(int num);
	double getTraceLUT(unsigned char type);
//	I/O functions
	bool ReadFromSWCFile(char * filename);
	bool ReadFromRPIXMLFile(char * filename);
	bool ReadFromSuperellipseXML(char * filename);
	void ReadFromVTKFile(char * filename);
	bool ReadFromFeatureTracksFile(char *filename, int type_offset);
	bool ReadFromFeatureTracksFileForKymograph(char *filename,int type_offset);
	bool WriteToSWCFile(const char * filename);
	void WriteToVTKFile(const char * filename);
	void SetBranchPoints(std::vector<branchPT*> Branches);
	void SetTraceOffset(double ntx, double nty, double ntz);
	void ImageIntensity(vtkSmartPointer<vtkImageData> imageData);
//	operators
	int getNewLineId();
  int GetMaximumBitId();
	void splitTrace(int selectedCellId);
	void ReverseSegment(TraceLine*);
	void RemoveTraceLine(TraceLine*);
	void FixPointMarkers(TraceLine* tline);
	void mergeTraces(unsigned long long int eMarker, unsigned long long int sMarker);
	void CreatePolyDataRecursive(TraceLine* , vtkSmartPointer<vtkFloatArray> , vtkSmartPointer<vtkPoints> ,vtkSmartPointer<vtkCellArray>);
	void FindMinLines(int smallSize);
	void cleanTree();
	void explode(TraceLine* parent);
	int createGapLists(std::vector<TraceLine*> traceList);
	int solveParents(std::vector<int> ids);
  void SetCombineShortVTKLines(bool b) { this->CombineShortVTKLines = b; }
//	public data
	vtkSmartPointer<vtkPolyData> GetVTKPolyData();
	vtkSmartPointer<vtkPolyData> generateBranchIllustrator();
	void Print(std::ostream &c);

	std::vector<TraceLine*>* GetTraceLinesPointer();
	std::vector<TraceLine*> GetTraceLines();
	std::vector<TraceBit> CollectTraceBits();
	std::set<long int> SmallLines;
	std::vector<TraceGap*> Gaps;
	std::vector<branchPT*> BranchPoints;
	std::vector<std::string> FeatureHeaders;
  vtksys::hash_map<unsigned int, unsigned long long int> hash_load;

	double gapTol;
	int gapMax;

	vtksys::hash_map< unsigned int, unsigned long long int > hashp;
	vtksys::hash_map< unsigned int, unsigned long long int > hashc;
	
private:
	std::vector<TraceLine*> trace_lines;
	vtkSmartPointer<vtkPolyData> PolyTraces;
	double smallLineColor, mergeLineColor;	
	double tx,ty,tz;
	int unsolvedBranches;
	bool isParent(int id);
  void CollectTraceBitsRecursive(std::vector<TraceBit> &vec,TraceLine *l);
  void CollectIdsRecursive(std::vector<int> &ids, TraceLine* tline);
  void LinearTraceLinesRecursive(std::vector<TraceLine*> &allLine, TraceLine* tline);
  void CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells,TraceLine *tline);
  void CollectSegmentMidPointsRecursive(vtkSmartPointer<vtkPoints>p, vtkSmartPointer<vtkCellArray> cells, vtkSmartPointer<vtkFloatArray> da,TraceLine* tline);

  //helper functions & objects for converting .vtk polydata into a TraceObject
  void ConvertVTKDataToTraceLines();
  void FindVTKTraceEnds();
  int VTKLineIsTraceEnd(int rootID);
  int VTKLineContainsPoint(int cellID, double point[3]);
  void ConvertVTKLineToTrace(int cellID, int parentTraceLineID, double *endPoint);
  vtkSmartPointer<vtkPolyData> VTKData;
	std::vector< std::pair<int, double *> > VTKTraceEnds;
  vtkSmartPointer<vtkIdTypeArray> DisconnectedVTKLines;
  //WARNING: this is only guaranteed to be accurate in the context of the above function.
  int NextTraceBitID;  
  bool CombineShortVTKLines;
	std::vector<branchPT*> branchPTsInProgress;

	
};

#endif

