

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <functional>
#include <stdlib.h>
#include <queue>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

#include <omp.h>

#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkImageActor.h"
#include "vtkImageData.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkIndent.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkBoostPrimMinimumSpanningTree.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPolyLine.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkLine.h"
#include "vtkUnsignedCharArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkTree.h"
#include "vtkVersion.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkVectorText.h"
#include "vtkBoostConnectedComponents.h"
#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkGraph.h"
#include "vtkIntArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkDirectedGraph.h"
#include "vtkOutEdgeIterator.h"
#include "vtkInEdgeIterator.h"

#include "itkTimeProbe.h"
#include "itkStatisticsImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkDivideImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImageFileWriter.h"

#include "Common.h"

/**
 * This class implements a vessel tracer using a spherical structuring element for tracing the vasuclature 
 * from 3D microscopy images. This is a C++ implementation of Amit's vessel tracing work. 
 *
 * Major steps of the algorithm are:
 * 
 * Pre-processing:
 *  1. Calculate gradient vector field
 *  2. Median filter (and flat-field correction?) of the input data
 *  3. Save images as *.mhd files
 * 
 * Part 1:
 *	1. Generate primary nodes (initialization using a grid) - using energy minimization for the spherical model
 *	2. Generate secondary nodes - using ODFs and additional constraints
 *	Result: Final set of nodes
 *
 * Part 2:
 *	1. Generate affinity graph
 *	2. Create MST from graph
 *	3. Prune the tree - based on heuristics
 * 
 * Post processing:
 * 1. Complete loops and generate a vascular mask (using a cylinder model)
 * Result: Final tracing of vasculature 
 */

typedef float PixelType;
typedef itk::Image<PixelType, 3> ImageType3D;
typedef unsigned char RenderPixelType;
typedef itk::Image<RenderPixelType, 3> RenderImageType3D;
typedef itk::ImageToVTKImageFilter<RenderImageType3D> ITKToVTKConnectorType; 
typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
typedef itk::ImageDuplicator<ImageType3D> DuplicatorType;
typedef itk::MaximumProjectionImageFilter<RenderImageType3D, RenderImageType3D> MaxProjectionFilterType;
typedef itk::MinimumProjectionImageFilter<RenderImageType3D, RenderImageType3D> MinProjectionFilterType;
typedef itk::ShiftScaleImageFilter<ImageType3D,ImageType3D> ShiftScaleFilterType;
typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VolumeOfInterestFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
typedef itk::DivideImageFilter<ImageType3D, ImageType3D, ImageType3D> DivideImageFilterType;
typedef itk::InvertIntensityImageFilter<ImageType3D> InvertImageFilterType;
typedef itk::ImageFileWriter<RenderImageType3D> ImageWriter;

typedef std::vector<int> VectorType1D;
typedef std::vector<VectorType1D> VectorType2D;
typedef std::vector<VectorType2D> VectorType3D;
typedef std::pair<int, double> queue_element;

struct SphericalBinInfo{

	//ArrayType3D BinIndex;
	VectorType3D binIndexVec;
	VectorType2D nbr;
	std::vector<std::vector<double> > binCenters;
	int indexLength;
	int angleCount;
	int angleInrement;
	int nLastIndicesOfInterest;
	double histSmoothingFactor;
	double minSphHistCount;

	void initByDefaultValues(void);
};

struct PreprocessingParameters{

	int medianFilterRadius;
	int anisDiffusionNIter;
	int anisDiffusionConductance;
	double smoothingSigma;
	int iterNGVF;
	void initByDefaultValues(void);
};

struct NodeDetectionParameters{

	int gridSpacing; // Spacing for uniform square grid 

	int iterNPrimaryNode; // Number of iterations required for the the energy minimization

	double increaseLikelihoodThreshold; // The likelhood of nodes less than this value is increased 
	double discardNodeLikelihoodThreshold; // Nodes with likelihood less than this value are discarded
	// If the number of iterations (while fitting a sphere to the node) are below this value, only region based term is used to update the energy
	int iterNForOnlyRegionBasedTerm;
	int iterNMinimum; // Minimum number of iterations
	double regionBasedTermWeight;
	double edgeBasedTermWeight;
	double minimumAccumulatedParamChange; // If the change in parameters falls below this value, exit model fitting
	int iterNMonitorParamChange; // The number of iterations for which the chnages in parameters are monitored and accumulated 
	
	double dtX, dtY, dtZ, dtScale; // The rate of change of node parameters
	int dirX, dirY, dirZ, dirScale; // The direction of change of the node parameters (positive, negetive or zero)
	std::vector<double> chX, chY, chZ, chScale; // The magnitude of change of the node parameters
	int currentMonitoredIter; // The current iteration for which the changes in parameters is recorded. Ranges from 1 to itetNMonitorParamChange 

	int maxVesselWidth, minVesselWidth; // Bounds on the node scale (in pixels)

	double likelihoodThresholdPrimary; // Nodes with likelihood less than this value will be discarded
	// A node will be discarded if its distance from any other node is lesser than this number times the distance between their centers
	double distanceThresholdPrimary;
	double distanceThresholdSecondary; // Same as distanceThresholdPrimary, but applied in the scheduling for secondary nodes

	double traceQualityThreshold;
	int maxQueueIter;
	double maxBranchAngle; // Maximum angle between branches of any two vessels
	double branchingThreshold;
	// Rate of change of secondary parameters
	double secondaryParamRateReduction;
	double dtXSecondary;
	double dtYSecondary;
	double dtZSecondary;
	double dtScaleSecondary;
	int iterNMinimumSecondary;
	double regionBasedTermWeightSecondary;
	double edgeBasedTermWeightSecondary;
	double secondarySearchConstraint;
	// Rates at which the direction of model fitting is reversed if energy and parameters point in opposite directions
	double primaryReversePositionRate;
	double primaryReverseScaleRate;
	double secondaryReversePositionRate;
	double secondaryReverseScaleRate;
	double maxTraceCost;
	double traceLengthCost;
	double primaryNodeSearchRadFactor; // Defines the factor of primary node radius which contributes to the initial position of secondary nodes
	double infTraceQuality;
	int maxQueueSize;

	void initByDefaultValues(void);
};

struct GraphAndMSTPartameters{

	double affinityRadThresh; // Distance threshold for selecting neighbouring nodes
	int NBinsAffinity; // number of bins required for binning neighborhood nodes based on angle
	double maxEdgeWeight; 
	double minBranchAngle; // minimum branching angle in radian
	int maxNBranches; 
	int maxTreeNodes;

	void initByDefaultValues(void);
};

struct AllParameters{
	
	PreprocessingParameters preProcessingParams;
	SphericalBinInfo oriBin;

	NodeDetectionParameters nodeDetectionParams;

	GraphAndMSTPartameters graphAndMSTParams;
	
	void initByDefaultValues(void);
};

struct GlobalStatisticsInput{
	
	double volumeMean;
	double volumeStd;
	double volumeMax;
	double volumeMin;
};

struct AffinityEdge{
	
	int from; 
	int to;
	double weight;

	AffinityEdge();
	AffinityEdge(int, int, double);
};

struct Tree{

	int ID;
	int start;
	int NNodes;

	Tree();
	Tree(int, int, int);
};

class SWCNodeVessel{

public:
	int ID;
	itk::Index<3> position;
	double scale;
	double likelihood;
	int label;
	bool isLeaf;
	bool isBifurgation;
	bool isTrifurgation;
	bool isMultifurgation;
	bool isOrphan;
	bool isRoot;

	std::vector<int> parents;
	std::vector<int> children;

	SWCNodeVessel();
};

struct Node{

	double x;
	double y;
	double z;
	PixelType intensity;
	double scale;
	double likelihood;
	bool isValid;
	double nHoodScale;

	static const int DEFALUT_SCALE = 4.0;
	static const int MIN_LIKELIHOOD = 0;
	
	// The surrounding (and inclusive) region of a node which contributes to the energy functional. 
	// This region is like a band around the sphere's circumference.
	double bandNhood;
	int bandArea;
	int foregroundArea;
	int backgroundArea;
	std::vector<double> xNormalizedInBand;
	std::vector<double> yNormalizedInBand;
	std::vector<double> zNormalizedInBand;
	std::vector<double> intensityInBand;
	std::vector<double> gxInBand;
	std::vector<double> gyInBand;
	std::vector<double> gzInBand;
	PixelType meanForegroundIntensity;
	PixelType meanBackgroundIntensity;
	int exitIter; // The number of iterations completed before exiting model fitting

	double traceQuality; // Quality of the trace to which the node belongs (eq. 5.6 Amit thesis)
	std::vector<double> parentID;
	int parentIDLength; // Last parentIDLength nodes are used for computing the trace quality
	// Indicate whether the node is CURRENTLY a primary or a secondary node
	bool isPrimary; 
	bool isSecondary;
	std::vector<double> sphHistRegionBased; // Spherical histogram binned by directions for determining the search space for secondary nodes
	// Vectors for storing the optimal directions of the above histogram 
	std::vector<double> dirX;
	std::vector<double> dirY; 
	std::vector<double> dirZ;
	double nHoodSecondaryMultiplier;
	double nHoodScaleSecondary; 

	double secondaryNodeSearchRad; // Defines the radius for the search space of secondary nodes
	// The coordinates at which the secondary node for this primary node will be initialized
	double xInitSecondary;
	double yInitSecondary;
	double zInitSecondary;


	// Stores the affinity connections for a node based on angular binning for neighbors (see Amit thesis pg. 126)
	std::vector<double> connectedNodesBinned;
	
	//Stores the filtered affinity connections for each node (affinity graph)
	std::vector<std::pair<int, double> > connectedNodesAffinity;
	
	int ID; 
	std::vector<int> branchIDs;
	int NBranches;

	int forestLabel;
	bool isRoot;

	Node();
	Node(double, double, double, PixelType);

	static double ComputeNorm(Node);
	static void ComputeDistanceBetweenNodes(Node, Node, Node&);
	void NormalizeNode(double);
	void InvertNodeDir(void);
	static double DotProduct(Node, Node);
	void InitDefaultParamsBeforeOptimization(void);
	void InitDefaultParamsBeforeODFRecursion(void);
};

class compareNodes{

public:
	compareNodes(const int& type=2){
		this->comparisonType = type;
	}
	bool operator()(Node& n1, Node& n2){
		if(this->comparisonType == 1)
			return (n1.likelihood > n2.likelihood);
		if(this->comparisonType == 2)
			return (n1.traceQuality > n2.traceQuality);
	}

private:
	/* Type 1: Comparison based on likelihood 
	 * Type 2: Comparison based on quality (of the trace to which the node belongs)
	 */
	int comparisonType; 
};
typedef std::priority_queue<Node, std::vector<Node>, compareNodes> PriorityQueueType;

class arrayElement{

public:
	arrayElement(double element, int index){
		this->element = element;
		this->index = index;
	}
	static bool CompareArrayElementsDescendPreserveIndex(arrayElement a1, arrayElement a2){
		return(a1.element > a2.element);
	}
	double getElement(void){
		return this->element;
	}
	int getIndex(void){
		return this->index;
	}
private:
	double element;
	int index;
};

class ftkVesselTracer{

public:

	ftkVesselTracer();
	ftkVesselTracer(std::string, bool, bool);
	~ftkVesselTracer();

	/** Preprocessing of data: 1. Reading the data from TIFF files 2. Median filtering 3. Edge enhancing
	 * 4. GVF computations 5. Saving all data and GVF as MHD files
	 * (data path with extension, empty data pointer)
	 */
	int PreprocessData(std::string, ImageType3D::Pointer&);
	
	/** Load gx, gy, gz, data and oriBin for further processing.
	 * (data path)
	 */
	void LoadPreprocessedData(std::string);

	/** Returns maximum intensity projection image
	 */
	void ComputeIntensityProjectionImages(void);

	/** Render MIP image
	 */
	void RenderMaximumProjectionImage(void);

	/** Compute the spherical bins and related info in preprocessing
	 * (indexLength)
	 */
	void SphericalBinPreprocess(void);

	/** Compute all primary nodes
	 */
	void ComputeAllPrimaryNodes(void);

	/** Compute seeds to initialize the tracing. This is done using a grid of a given spacing.
	 */
	void ComputeSeeds(void);

	/** Fit a sphere at all nodes (seeda) using the energy minimization technique (Amit thesis pp. 112)
	 * and sort nodes using the likelihood.
	 */
	void FitSphereAndSortNodes(void);
	
	/** Visualize nodes on the data 3D
	 * (nodes vector, visualize nodes as a point?)
	 */
	void VisualizeNodesWithData3D(std::vector<Node>, bool);

	/** Return the sign of the value - THIS SHOULD MOVE TO COMMON
	 * (value)
	 */
	int inline GetSign(double);

	/* Compute all secondary nodes
	 */
	void ComputeAllSecondaryNodes(void);

	/* Compute all secondary nodes
	 */
	void ComputeAllSecondaryNodes2(void);

	/* Write a vector of nodes to a text file for further processing
	 * (Vector of nodes to write, file path with extension)
	 */
	static void writeNodesToFile(std::vector<Node>&, std::string);

	/* Read nodes from text file 
	 */
	void ReadNodesFromTextFile(const std::string&);

	/* Create MST 
	 */
	void CreateMinimumSpanningForest(void);

	/* Create affinity graph
	 */
	void CreateAffinityGraph(void);

	/* Visualise affinity graph
	 */
	void VisualizeAffinityGraph(bool render_with_data);

	/* Compute MST using Boost/VTK
	 * NOT IMPLEMENTED YET
	 */
	void ComputeMinimumSpanningTreeBoost(void);

	/* Compute minimum spanning forest with loop detection.
	 * Detects the loops while generating a forest
	 */
	void ComputeMinimumSpanningForestWithLoopDetection(void);
	
	/* Fetch the tree at the given ID 
	 * (ID, connected IDs vector)
	 */
	void GetTree(int, std::vector<int>&);

	/* Visualize Minimum spanning forest
	 */
	void VisualizeMinimumSpanningForest(bool render_with_data);

	/* Find a element in the list upto the given index
	 */
	static bool IsInList(std::vector<int>, int, int);

	/* Print the current forest
	 */
	void PrintForest(void);

	/* Fit sphere at a single node
	 * (Node object)	
	 */
	void FitSphereAtNode(Node&);

	/* Fit sphere at a single node
	 * (Node object, image_ptr, gx ptr, gy, ptr, gz ptr)	
	 */
	void FitSphereAtNode(Node&, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer);

	void InitNodeDetectionParamsDefault(void);

	void PopulateSWCNodeContainer(void);

	void WriteSWCFileVessel(const std::string);  

	void PrintToySWCFile(void);

private:

	AllParameters allParams;
	GlobalStatisticsInput globalStatsInput;

	std::string data_folder_path;
			
	ImageType3D::Pointer originalData;
	ImageType3D::Pointer inputData;
	ImageType3D::Pointer normalizedInputData;
	ImageType3D::Pointer gx;
	ImageType3D::Pointer gy;
	ImageType3D::Pointer gz;

	RenderImageType3D::Pointer originalDataForRendering;
	RenderImageType3D::Pointer inputDataForRendering;
	RenderImageType3D::Pointer maximumProjectionImage;
	RenderImageType3D::Pointer minimumProjectionImage;
	RenderImageType3D::Pointer primaryNodesImage, secondaryNodesImage;
	
	std::vector<Node> initialSeeds;
	std::vector<Node> primaryNodes;
	std::vector<Node> primaryNodesAfterHitTest;
	std::vector<Node> allNodes;
	std::vector<Node> allForestNodes;

	std::vector<AffinityEdge> edges;
	std::vector<AffinityEdge> loops;
	std::vector<Tree> forest;
	std::vector<SWCNodeVessel> SWCNodeVessel_vec;

	/* Update the node appearance 
	 * (Node object)
	 */
	void UpdateAppearanceVectorized(Node&);

	/* Update the node appearance 
	 * (Node object, image ptr, gx ptr, gy ptr, gz ptr)
	 */
	void UpdateAppearanceVectorized(Node&, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer);

	/* Update the node model
	 * (Node object, iteration number)
	 */
	void UpdateModel(Node&, int);

	/* Define the criteria for exiting the model fitting at seeds
	 * (Node object, iteration number)
	 */
	bool ExitModelFitting(Node&, int);

	/* Sort the nodes based on likelihood and filter them based on 
	 * Euclidian distance (nodes which are too close are filtered out)
	 */
	void SortAndFilterPrimaryNodes(void);

	/* While scheduling the tracer, determines wheather or not the image region is 
	 * occupied by a current node
	 * (Node object)
	 */
	int TraceHitTest(Node); 

	/* Compute the directions for secondary nodes for a primary node
	 * (Node ref, direction hist ref)
	 */
	void ComputeSecondaryNodeDirections(Node&, std::vector<double>&);
	
	/* Fit the model using additional constraints for secondary nodes for the given direction vector
	 * (Primary/anchor node ref, Secondary node ref, direction vec)
	 */
	void FitSphereAtNodeSecondary(Node&, Node&, std::vector<double>);

	/* Update the node model for secondary nodes
	 * (Node ref, anchor node ref, iteration number)
	 */
	void UpdateModelSecondary(Node&, Node&, int);

	/* Compute the cost of the trace for the given node
	 * (Node ref)
	 */
	double computeTraceQuality(Node&);

	/* Compute the bin for affinity graph
	 * (direction node)	
	 */
	int ComputeAffinityBin(Node);

	/* Relabel nodes when there is a collision in node IDs
	 * (old ID, new ID)
	 */
	void RelabelForestNodes(int, int);

	/* Check if already labeled nodes are neighbors in the tree
	 * (node 1 ID, node 2 ID)
	 */
	bool CheckNeighbors(int, int);

	/* Sort the queue based on trace quality and return the best trace index and quality
	 */
	void GetBestTrace(std::vector<queue_element>&, queue_element&);
};
