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
#include <map>
#include <vtksys/hash_map.hxx> /* Platform independent hashmap */
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkIdTypeArray.h"

class TraceBit;
class TraceLine;
class TraceGap;
class branchPT;
class CellTrace;

class vtkPoints;
class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;


//needed because gcc doesn't have a built-in method to hash unsigned long long ints
struct hashulli
{
	size_t operator()(const unsigned long long int __x) const
	{
		return __x;
	}
	size_t operator()(const unsigned long long int __x, const unsigned long long int __y)
	{
		return __x == __y;
	}
	const static size_t bucket_size = 4;
	const static size_t min_buckets = 8;
};

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
	void setFalseLineColor(double set)
	{
		this->falseLineColor=set;
	};
	void setLUT(int num);
	double GetTraceLUT(TraceLine *line);
	double GetNodeTypeLUT(unsigned char type);
	double GetTreeLUT(int whichTree);
	std::vector<CellTrace*> CalculateCellFeatures();
	//  I/O functions
	bool ReadFromSWCFile(char * filename);
	bool ReadFromRPIXMLFile(char * filename);
	//  bool ReadFromSuperellipseXML(char * filename);//use convert to swc instead
	void ReadFromVTKFile(char * filename);
	bool ReadFromFeatureTracksFile(char *filename, int type_offset);
	bool ReadFromFeatureTracksFileForKymograph(char *filename,int type_offset);
	bool WriteToSWCFile(const char * filename);
	bool WriteToSWCFile(std::vector<TraceLine *> selectedLines, const char * filename);
	void WriteToVTKFile(const char * filename);
	void SetBranchPoints(std::vector<branchPT*> Branches);
	void SetTraceOffset(double ntx, double nty, double ntz);
	void ImageIntensity(vtkSmartPointer<vtkImageData> imageData);
	//  operators
	int getNewLineId();
	int GetMaximumBitId();
	void splitTrace(int selectedCellId);
	void ReverseSegment(TraceLine*);
	void RemoveTraceLine(TraceLine*);
	void FixPointMarkers(TraceLine* tline);
	void mergeTraces(unsigned long long int eMarker, unsigned long long int sMarker);
	void CreatePolyDataRecursive(TraceLine* , vtkSmartPointer<vtkFloatArray> , vtkSmartPointer<vtkPoints> ,vtkSmartPointer<vtkCellArray>);
	void FindMinLines(int smallSize);
	void FindFalseSpines(int maxBit, int maxLength);
	void FindFalseBridges(int maxBit);
	void FindHalfBridges(int maxBit, int DtoParent);
	void cleanTree();
	void Shave(TraceLine* starting, int smallerThan);
	bool BreakOffBranch(TraceLine* branch, bool keep);
	void explode(TraceLine* parent);
	int createGapLists(std::vector<TraceLine*> traceList);
	int createBranchPtFromList(std::vector<TraceLine*> traceList);
	int solveParents(std::vector<int> ids);
	bool isParent(int id);
	void SetCombineShortVTKLines(bool b) { this->CombineShortVTKLines = b; }
	//manual tracing tools
	void createSomaFromPT(double pt[],std::vector<TraceLine*> stems);
	TraceBit CreateBitAtCoord(double pt[]);
	TraceLine* CreateTraceFromBit(TraceBit firstBit);
	void ExtendTraceTo(TraceLine* tline, double pt[]);
	//  public data
	vtkSmartPointer<vtkPolyData> GetVTKPolyData();
	vtkSmartPointer<vtkPolyData> generateBranchIllustrator();
	void Print(std::ostream &c);

	TraceLine* findTraceByID(int id);
	std::vector<TraceLine*>* GetTraceLinesPointer();
	std::vector<TraceLine*> GetTraceLines();
	std::vector<TraceBit> CollectTraceBits();
	std::vector<int> GetTreeIDs( TraceLine* root);
	std::vector<int> GetTreeIDs( std::vector<TraceLine*> roots);
	std::set<long int> SmallLines;
	std::set<long int> FalseSpines;
	std::set<long int> FalseBridges;
	std::set<long int> HalfBridges;
	std::vector<TraceGap*> Gaps;
	std::vector<branchPT*> BranchPoints;
	std::vector<std::string> FeatureHeaders;
	vtksys::hash_map<unsigned int, unsigned long long int> hash_load;

	double gapTol;
	int gapMax;
	bool AutoSolveBranchOrder;

	vtksys::hash_map< unsigned int, unsigned long long int > hashp;
	vtksys::hash_map< unsigned int, unsigned long long int > hashc;
	std::vector<TraceBit> debug_points;// ADDED BY ARUN ON 22nd August 2010 for debugging
	//the following functions are used by "color by trees" mode
	void SetColorByTrees(bool b) { this->ColorByTrees = b; }
	bool GetColorByTrees() { return this->ColorByTrees; }
	std::vector<TraceLine*> GetAllRoots();
	void UpdateRootToTree();
	void RecolorTraces();
	void RecolorTrace(TraceLine *line);

private:
	std::vector<TraceLine*> trace_lines;
	std::vector<CellTrace*> Cells;
	vtkSmartPointer<vtkPolyData> PolyTraces;
	double smallLineColor, mergeLineColor, falseLineColor;  
	double tx,ty,tz;
	int unsolvedBranches;
	void CollectTraceBitsRecursive(std::vector<TraceBit> &vec,TraceLine *l);
	void CollectIdsRecursive(std::vector<int> &ids, TraceLine* tline);
	int LinearTraceLinesRecursive(std::vector<TraceLine*> &allLine, TraceLine* tline);
	void CollectBranchPointsRecursive(vtkSmartPointer<vtkPoints> p, vtkSmartPointer<vtkCellArray> cells,TraceLine *tline);
	void CollectSegmentMidPointsRecursive(vtkSmartPointer<vtkPoints>p, vtkSmartPointer<vtkCellArray> cells, vtkSmartPointer<vtkFloatArray> da,TraceLine* tline);

	void ParseFileName(char* fullName);
	std::vector<std::string> ParsedName;
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

	//used for "color by trees" mode
	bool ColorByTrees;
	std::map<int, int> RootToTree;
};

#endif

