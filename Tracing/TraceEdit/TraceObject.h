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
#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkIdTypeArray.h"

#include "vtkCellArray.h"
#include "vtkProperty.h"
#include "vtkDataSetMapper.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkMath.h"
#include "vtkDelaunay3D.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"

#include "itkImage.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#ifdef USE_BALL_TRACER
#include "Tracing/ftkVesselTracer/ftkVesselTracer.h"
#endif

#ifdef USE_GT_CLUSTERING
#include "PatternAnalysis/ftkGTClustering/ftkGTClustering.h"
#endif

class TraceBit;
class TraceLine;
class TraceGap;
class branchPT;
class CellTrace;

class vtkPoints;
class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;

#define TRACE_TYPE_TREE 1
#define TRACE_TYPE_GRAPH 2
#define TRACE_TYPE_BOTH 3 // For future use

typedef itk::Image< unsigned char, 3 >   ImageType;
typedef itk::Image< float, 3> FloatImageType;

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
	std::map< int ,CellTrace*> CalculateCellFeatures();
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
	
	void AddTraceLine(int selectedCellId1, int selectedCellId2);
	TraceLine* AddTraceLine(double p1[],double p2[], double color);
	void AddEndTraceLine(int selectedCellId1, int selectedCellId2);
	void AddExtensionToTraceLine(int selectedCellID, double point[]);

	void addTrace(TraceLine* traceToAdd);
	bool removeTrace(TraceLine* traceToRemove);
	void markRootAsModified(int RootID);

	void ImageIntensity(vtkSmartPointer<vtkImageData> imageData);
	void ImageWeightedIntensity(ImageType::Pointer intensityImage);
	//  operators
	int getNewLineId();
	int GetMaximumBitId();
	void splitTrace(int selectedCellId);
	void ReverseSegment(TraceLine* tline);
	void RemoveTraceLine(TraceLine* tline);
	void FixPointMarkers(TraceLine* tline);
	void mergeTraces(unsigned long long int eMarker, unsigned long long int sMarker);
	void UnmarkLines(TraceLine* gtline);
	void CreatePolyDataRecursive(TraceLine* tline, vtkSmartPointer<vtkFloatArray> point_scalars, vtkSmartPointer<vtkPoints> line_points,vtkSmartPointer<vtkCellArray> line_cells);
	void CreatePolyDataRecursive(TraceLine* tline, vtkSmartPointer<vtkPoints> line_points,vtkSmartPointer<vtkCellArray> line_cells);
	void FindMinLines(int smallSize);
	void FindFalseSpines(int maxBit, int maxLength);
	void FindFalseBridges(int maxBit);
	void FindHalfBridges(int maxBit, int DtoParent);
	void DetectCandidateGapsVBT();
	void RunClusteringGTRepDyn();
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
	vtkSmartPointer<vtkPolyData> GetVTKPolyData(bool bSetScalar = true);
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
	std::set<long int> RetracingLines;
	std::map<unsigned int, TraceBit> RetracingPoints;
	std::vector<TraceGap*> Gaps;
	std::vector<unsigned int>RootIDs;
	std::vector<branchPT*> BranchPoints;
	std::vector<std::string> FeatureHeaders;
	vtksys::hash_map<unsigned int, unsigned long long int> hash_load;

	double gapTol;
	double gapAngleTol;
	int gapMax;
	bool AutoSolveBranchOrder;

	double maxGapLength;
	double minGapAngle;
	double minGapMutualAngle;
	double alphaForClustering;

	FloatImageType::Pointer gx, gy, gz, vesselnessImg, inputImageVBT;

	vtksys::hash_map< unsigned int, unsigned long long int > hashp;
	vtksys::hash_map< unsigned int, unsigned long long int > hashc;
	//std::vector<TraceBit> debug_points;// ADDED BY ARUN ON 22nd August 2010 for debugging
	//the following functions are used by "color by trees" mode
	void SetColorByTrees(bool b) { this->ColorByTrees = b; }
	bool GetColorByTrees() { return this->ColorByTrees; }
	std::vector<TraceLine*> GetAllRoots();
	void UpdateRootToTree();
	void RecolorTraces();
	void RecolorTrace(TraceLine *line);

	void SetTraceTypeGeneric(int type);
	void SetTraceTypeGeneric(std::string type);
	int GetTraceTypeGeneric();
	void SetGVFImagesVBT(); 
	void SetVesselnessImageVBT();
	void SetInputImageVBT();
	void RemoveCandidateAndClusteredTraceLines();

private:
	std::vector<TraceLine*> trace_lines;
	std::map< int ,CellTrace*> Cells;
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
	
	int TraceTypeGeneric;
	std::vector<TraceGap*> CandidateGaps;
	std::vector<TraceGap*> ClusteredGaps;
	std::vector<TraceLine*> CandidateTLines;
	std::vector<TraceLine*> ClusteredTLines;
	
#ifdef USE_BALL_TRACER
	ftkVesselTracer *VBT;
#endif

#ifdef USE_GT_CLUSTERING
	ftkGTClustering *GTCluster;
#endif
};

#endif

