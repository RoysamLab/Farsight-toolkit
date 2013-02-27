#ifndef _ftkVesselTracer_h_
#define _ftkVesselTracer_h_
 
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
#include <cmath>

#ifdef _OPENMP
	#include <omp.h>
#endif

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
#include "vtkTable.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkImageCanvasSource2D.h"
#include "vtkTIFFWriter.h"

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
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkPoint.h"
#include "itkLineIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_numeric_traits.h>
#include <vnl/vnl_math.h>
#include <vnl/algo/vnl_convolve.h>
#include <vnl/algo/vnl_gaussian_kernel_1d.h>

#include <vcl_complex.h>
#include <mbl/mbl_stats_nd.h>

#include "boost/multi_array.hpp"
#include "boost/lexical_cast.hpp"

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
typedef unsigned char RenderPixelType;
typedef unsigned short LabelPixelType;
typedef itk::Image<PixelType, 3> ImageType3D;

typedef itk::SymmetricSecondRankTensor<double, 3> HessianPixelType;
typedef itk::Image< HessianPixelType, 3 > HessianImageType3D;

typedef itk::Image<RenderPixelType, 3> RenderImageType3D;
typedef itk::Image<LabelPixelType, 3> LabelImageType3D; 
typedef itk::ImageToVTKImageFilter<RenderImageType3D> ITKToVTKConnectorType; 
typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
typedef itk::ImageDuplicator<ImageType3D> DuplicatorType;
typedef itk::MaximumProjectionImageFilter<RenderImageType3D, RenderImageType3D> MaxProjectionFilterType;
typedef itk::MinimumProjectionImageFilter<RenderImageType3D, RenderImageType3D> MinProjectionFilterType;
typedef itk::ShiftScaleImageFilter<ImageType3D,ImageType3D> ShiftScaleFilterType;
typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VBTVolumeOfInterestFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
typedef itk::DivideImageFilter<ImageType3D, ImageType3D, ImageType3D> DivideImageFilterType;
typedef itk::InvertIntensityImageFilter<ImageType3D> InvertImageFilterType;
typedef itk::ImageFileWriter<RenderImageType3D> ImageWriter;
//typedef itk::HessianToObjectnessMeasureImageFilter<PixelType, 3> ObjectnessFilterType;
typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType3D, HessianImageType3D > MultiScaleHessianFilterType;
typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType3D, ImageType3D > ObjectnessFilterType;


typedef itk::AddImageFilter<ImageType3D> VBTAddImageFilterType;
typedef itk::MultiplyImageFilter<ImageType3D> MultiplyImageFilterType;
typedef itk::NormalizeImageFilter<ImageType3D, ImageType3D> NormalizeImageFilterType;
typedef itk::BresenhamLine<3> LineType3D;
typedef itk::LineIterator<RenderImageType3D> LineIteratorType3D;
typedef itk::SignedMaurerDistanceMapImageFilter<RenderImageType3D, ImageType3D> SignedMaurerDistanceMapImageFilterType;
typedef itk::BinaryThresholdImageFilter<LabelImageType3D, RenderImageType3D> ThresholdFilterType;
typedef itk::LabelGeometryImageFilter<LabelImageType3D, RenderImageType3D> LabelGeometryFilterType;
typedef itk::LabelOverlapMeasuresImageFilter<RenderImageType3D> LabelOverlapFilterType;
typedef itk::ConnectedComponentImageFilter<RenderImageType3D, RenderImageType3D> ConnectedComponentFilterType;

typedef std::vector<int> VectorType1D;
typedef std::vector<VectorType1D> VectorType2D;
typedef std::vector<VectorType2D> VectorType3D;
typedef std::pair<int, double> queue_element;

typedef boost::multi_array<int, 3> ArrayType3D;

struct SphericalBinInfo{

	//USE BOOST 3D ARRAY TO IMPLEMENT THIS CLASS

	//ArrayType3D BinIndex;
	VectorType3D binIndexVec;
	VectorType2D nbr;
	VectorType3D nbr2D;
	std::vector<std::vector<double> > binCenters;
	vnl_matrix<double> gaussian_prior;
	int indexLength;
	int angleCount;
	int angleInrement;
	int nLastIndicesOfInterest;
	double histSmoothingFactor;
	double minSphHistCount;
	double gaussian_sigma;
	int gaussian_win;

	void initByDefaultValues(void);
};

class VesselnessMeasures{
	
public:
	float alpha;
	float beta;
	float gamma;

	float sigma_min;
	float sigma_max;
	int sigma_intervals;
	int vesselness_type; //0: Blobness, 1: Vesselness, 2: Plateness
	
	float noiseness;
	float ballness;
	float plateness;
	float vesselness;

	VesselnessMeasures();
	VesselnessMeasures(float alpha, float beta, float gamma);
	VesselnessMeasures(float sigma_min, float sigma_max, float sigma_intervals, int obj_type);
};

struct PreprocessingParameters{

	int medianFilterRadius;
	int anisDiffusionNIter;
	int anisDiffusionConductance;
	double smoothingSigma;
	int iterNGVF;
	VesselnessMeasures vesselness_measures;

	void initByDefaultValues(void);
};

struct VBTNodeDetectionParameters{

	int gridSpacing; // Spacing for uniform square grid 

	int iterNPrimaryVBTNode; // Number of iterations required for the the energy minimization

	double increaseLikelihoodThreshold; // The likelhood of nodes less than this value is increased 
	double discardVBTNodeLikelihoodThreshold; // VBTNodes with likelihood less than this value are discarded
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

	double likelihoodThresholdPrimary; // VBTNodes with likelihood less than this value will be discarded
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
	double primaryVBTNodeSearchRadFactor; // Defines the factor of primary node radius which contributes to the initial position of secondary nodes
	double infTraceQuality;
	int maxQueueSize;

	double vesselnessThershold;
	double vesselnessWeight;

	double gaussianPriorSigma; // In angle units

	double traceLengthCostRetracing;
	double traceQualityThresholdRetracing;

	void initByDefaultValues(void);
};

struct GraphAndMSTPartameters{

	double affinityRadThresh; // Distance threshold for selecting neighbouring nodes
	int NBinsAffinity; // number of bins required for binning neighborhood nodes based on angle
	double maxEdgeWeight; 
	double minBranchAngle; // minimum branching angle in radian
	int maxNBranches; 
	int maxTreeVBTNodes;

	void initByDefaultValues(void);
};

struct AllParameters{
	
	PreprocessingParameters preProcessingParams;
	SphericalBinInfo oriBin;

	VBTNodeDetectionParameters nodeDetectionParams;

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

struct VBTTree{

	int ID;
	int start;
	int NVBTNodes;

	VBTTree();
	VBTTree(int, int, int);
};

class SWCVBTNodeVessel{

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

	SWCVBTNodeVessel();
};

// THIS COUPLE OF NEXT CLASSES TO BE MODIFIED AND TRANSFERRED TO TRACE EDITOR
class VesselNetworkFeatures{

public:
	double totalVesselLength;
	double totalVesselVolume;
	double totalVesselness;
	int nBifurgations;
	int nTrifurgations;
	int nMultifurfations;
	double meanRadius;
	double meanVesselSize;
	double meanVesselness;
	double vesselVolumePercentage;

	// Not worked on computing these
	double totalTortuosity;
	double meanTortuosity;
	int nLoops;
	double collaterality;
	double meanLoopRadius;
	
	// Add associative features


	VesselNetworkFeatures();
};

class VesselSegmentFeatures{

public:
	double segmentLength;
	double segmentVolume;
	double segmentVesselness;
	double meanRadius;
	double meanVesselSize;
	double meanVesselness;
	double meanCurvature;
	double totalTortuosity;
	double meanTortuosity;

	VesselSegmentFeatures();
};

class VesselVBTNodeFeatures{

public:
	int ID;
	itk::Index<3> position;
	double scale;
	double likelihood;
	int forestLabel;
	bool isLeaf;
	bool isRoot;
	bool isBifurgation;
	bool isTrifurgation;
	double traceQuality;
	int nODFModes;
	std::vector<double> ODFModesX;
	std::vector<double> ODFModesY;
	std::vector<double> ODFModesZ;
	
	VesselVBTNodeFeatures();
};

class ODFFeatures{

public:
	int nModes;
	std::vector<double> ODFModeVals;
	double std;
	double mean;
	double energy;
};

struct VBTNode{

	double x;
	double y;
	double z;
	PixelType intensity;
	double scale;
	double likelihood;
	double vesselness_likelihood;
	bool isValid;
	double nHoodScale;
	double vesselness;

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
	std::vector<double> vesselnessInBand;
	std::vector<double> gxInBand;
	std::vector<double> gyInBand;
	std::vector<double> gzInBand;
	PixelType meanForegroundIntensity;
	PixelType meanBackgroundIntensity;
	PixelType meanForegroundVesselness;
	PixelType meanBackgroundVesselness;
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

	double secondaryVBTNodeSearchRad; // Defines the radius for the search space of secondary nodes
	// The coordinates at which the secondary node for this primary node will be initialized
	double xInitSecondary;
	double yInitSecondary;
	double zInitSecondary;


	// Stores the affinity connections for a node based on angular binning for neighbors (see Amit thesis pg. 126)
	std::vector<double> connectedVBTNodesBinned;
	
	//Stores the filtered affinity connections for each node (affinity graph)
	std::vector<std::pair<int, double> > connectedVBTNodesAffinity;
	
	int ID; 
	std::vector<int> branchIDs;
	int NBranches;

	int forestLabel;
	bool isRoot;
	bool isLeaf;
	bool isBifurgation;
	bool isTrifurgation;
	bool isMultifurgation;
	bool isOrphan;
	std::vector<int> children;
	std::vector<int> parents;
	int ODF_modes;

	VesselVBTNodeFeatures nodeFeatures;
	ODFFeatures odfFeatures;

	itk::Index<3> gridNdx;

	VBTNodeDetectionParameters nodeDetectionParams;

	VBTNode();
	VBTNode(double, double, double, PixelType);

	static double ComputeNorm(VBTNode);
	static void ComputeDistanceBetweenVBTNodes(VBTNode, VBTNode, VBTNode&);
	void NormalizeVBTNode(double);
	void InvertVBTNodeDir(void);
	static double DotProduct(VBTNode, VBTNode);
	void InitDefaultParamsBeforeOptimization(void);
	void InitDefaultParamsBeforeODFRecursion(void);
	void SetLocationFromArray(double p[]);
	itk::Index<3> GetLocationAsITKIndex(void);
	void SetLocationFromITKIndex(itk::Index<3> idx);
};

class compareVBTNodes{

public:
	compareVBTNodes(const int& type=2){
		this->comparisonType = type;
	}
	bool operator()(VBTNode& n1, VBTNode& n2){
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

typedef std::priority_queue<VBTNode, std::vector<VBTNode>, compareVBTNodes> PriorityQueueType;

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


class IntrinsicFeatureVector_VT{

public:

	itk::Index<3> centroid;
	//HeapVBTNode_astro weightedCentroid;
	unsigned short int ID;
	double volume;
	double boundingBoxVolume;
	int integratedIntensity;
	double meanIntensity;
	double varianceIntensity;
	double eccentricity;
	double elongation;
	double meanSurfaceGradient;
	double radiusVariation;
	double shapeMeasure; // Ratio of surface voxels to total voxels
	double energy;
	double entropy;
	double inverseDiffMoment;
	double inertia;
	double clusterShade;
	double clusterProminence;
	
	IntrinsicFeatureVector_VT();
};

class AssociativeFeatureVector_VT{

public:

	double astro_total;
	double astro_avg;
	double astro_surr;
	double micro_total;
	double micro_avg;
	double micro_surr;
	double neuro_total;
	double neuro_avg;
	double neuro_surr;
	double vessel_total;
	double vessel_avg;
	double vessel_surr;

	double minRootDist;
	double maxRootDist;
	double meanRootDist;
	double varRootDist;
	double nRoots;

	double distToVessel;
	double alignmentWithVessel;
	double overlapWithVessel;
	//double minVesselDist;
	//double maxVesselDist;
	double meanVesselDist;
	double varVesselDist;
	int nVesselPoints;

	AssociativeFeatureVector_VT();

};

class NucleiObject_VT{

public:
	IntrinsicFeatureVector_VT intrinsicFeatures;
	AssociativeFeatureVector_VT associativeFeatures;
	
	// 1: Astrpcytes	2: Microglia	3:Neurons	4: Endotheliels
	int classValue; 
	double confidenceMeasure;

	NucleiObject_VT();
};
class ftkVesselTracer{

public:

	ftkVesselTracer();
	ftkVesselTracer(std::string, bool, bool, int);
	ftkVesselTracer(std::string, ImageType3D::Pointer, bool, bool, int);
	~ftkVesselTracer();

	/** Preprocessing of data: 1. Reading the data from TIFF files 2. Median filtering 3. Edge enhancing
	 * 4. GVF computations 5. Saving all data and GVF as MHD files
	 * (data path with extension, empty data pointer)
	 */
	int PreprocessData(std::string, ImageType3D::Pointer&, bool, bool);

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
	void ComputeAllPrimaryVBTNodes(void);

	/** Compute seeds to initialize the tracing. This is done using a grid of a given spacing.
	 */
	void ComputeSeeds(void);

	/** Fit a sphere at all nodes (seeds) using the energy minimization technique (Amit thesis pp. 112)
	 * and sort nodes using the likelihood.
	 */
	void FitSphereAndSortVBTNodes(void);

	/** Visualize nodes on the data 3D
	 * (nodes vector, visualize nodes as a point?)
	 */
	void VisualizeVBTNodesWithData3D(std::vector<VBTNode>, bool);

	/** Return the sign of the value - THIS SHOULD MOVE TO COMMON
	 * (value)
	 */
	int inline GetSign(double);

	/* Compute all secondary nodes
	 */
	void ComputeAllSecondaryVBTNodes(void);

	/* Compute all secondary nodes
	 */
	void ComputeAllSecondaryVBTNodes2(void);

	/* Compute all secondary nodes for retracing
	 */
	void ComputeAllSecondaryVBTNodesRetracing(void);

	/* Write a vector of nodes to a text file for further processing
	 * (Vector of nodes to write, file path with extension)
	 */
	static void writeVBTNodesToFile(std::vector<VBTNode>&, std::string);

	/* Read nodes from text file 
	 */
	void ReadVBTNodesFromTextFile(const std::string&);

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
	void GetVBTTree(int, std::vector<int>&);

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
	 * (VBTNode object)	
	 */
	void FitSphereAtVBTNode(VBTNode&);

	/* Fit sphere at a single node
	 * (VBTNode object, image_ptr, gx ptr, gy, ptr, gz ptr)	
	 */
	void FitSphereAtVBTNode(VBTNode&, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer);
	
	/* Internal function	
	 */
	void InitVBTNodeDetectionParamsDefault(void);

	/* Generate a labelled, directed graph from the unlabelled graph. Also generates preliminary node-based features
	 */
	void PopulateSWCVBTNodeContainerAndComputeVBTNodeFeatures(void);

	/* Write the SWC file to disk
	 */
	void WriteSWCFileVessel(void);  

	/* Print a toy SWC file to cout
	 */
	void PrintToySWCFile(void);

	/* Compute vesselness image (part of preprocessing)
	 */
	void ComputeVesselnessImage(VesselnessMeasures&, std::string&, ImageType3D::Pointer&);

	/* Pipeline function for smart tracing
	 */
	void SmartRetrace(void);

	/* Estimate points to begin retracing from
	 */
	void ComputeRetracingStartPoints(void);

	/* Write vessel skeleton image to disk
	 */
	void WriteSkeletonImage(void);

	/* Write vessel skeleton image to disk
	 */
	void WriteSkeletonImageFromVTK(void);

	/* Write the segmentation mask to disk as a label image
	 */
	void WriteSegmentationMask(void);
	
	/* Write node properties to file
	 */
	void WriteVBTNodeFeaturesFile(void);
	
	/* Compute vessel network features
	 */
	void ComputeVesselNetworkFeatures(void);

	/* Fit spheres along given line and return the cost of tracing
	 */
	std::vector<VBTNode> FitSpheresOnTraceLine(double p1[], double p2[]);

	/* Return the costs for tracing given nodes
	 */
	bool ComputeTracingCosts(double p1[], double p2[], double& tracing_cost, double& vesselness_cost, double& scale_var, double& vesselness_var);

	/* Call this function before calling any other processing function in this class
	*/
	void NormalizeAndRescaleData();

	/* Compute ODF on an isolated node
	*/
	void ComputeODFNoPrior(VBTNode& node);

	/* Compute some features from ODFs
	*/
	void ComputeODFFeatures(VBTNode& node);

	void SetInputImage(ImageType3D::Pointer input_img);
	void SetGVFImages(ImageType3D::Pointer gx, ImageType3D::Pointer gy, ImageType3D::Pointer gz);
	void SetVesselnessImage(ImageType3D::Pointer vesselness_img);
	void Set_useVesselness(int value);

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
	ImageType3D::Pointer VesselnessImage;
	ImageType3D::Pointer InputImageRetracing;

	RenderImageType3D::Pointer originalDataForRendering;
	RenderImageType3D::Pointer inputDataForRendering;
	RenderImageType3D::Pointer maximumProjectionImage;
	RenderImageType3D::Pointer minimumProjectionImage;
	RenderImageType3D::Pointer primaryVBTNodesImage, secondaryVBTNodesImage;
	RenderImageType3D::Pointer retracingStartPointsImage;
	RenderImageType3D::Pointer skeletonImage, segmentationMaskImage;

	std::vector<VBTNode> initialSeeds;
	std::vector<VBTNode> primaryVBTNodes;
	std::vector<VBTNode> primaryVBTNodesAfterHitTest;
	std::vector<VBTNode> allVBTNodes;
	std::vector<VBTNode> allForestVBTNodes;
	std::vector<VBTNode> retracingStartPoints;

	//std::map<itk::Index<3>, VBTNode, compareIndex> nodeGridMap;
	//ArrayType3D nodeGridArray;

	std::vector<AffinityEdge> edges;
	std::vector<AffinityEdge> loops;
	std::vector<VBTTree> forest;
	std::vector<SWCVBTNodeVessel> SWCVBTNodeVessel_vec;

	// 0: Use vesselness in selecting primary seeds only
	// 1: Use vesselness in case0 and tracing cost function
	// 2: Use vesselness in case1 and sphere cost function
	int useVesselness; 
	
	VesselNetworkFeatures networkFeatures;

	/* Update the node appearance 
	 * (VBTNode object)
	 */
	void UpdateAppearanceVectorized(VBTNode&);

	/* Update the node appearance 
	 * (VBTNode object, image ptr, gx ptr, gy ptr, gz ptr)
	 */
	void UpdateAppearanceVectorized(VBTNode&, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer, ImageType3D::Pointer);

	/* Update the node model
	 * (VBTNode object, iteration number)
	 */
	void UpdateModel(VBTNode&, int);

	/* Define the criteria for exiting the model fitting at seeds
	 * (VBTNode object, iteration number)
	 */
	bool ExitModelFitting(VBTNode&, int);

	/* Sort the nodes based on likelihood and filter them based on 
	 * Euclidian distance (nodes which are too close are filtered out)
	 */
	void SortAndFilterPrimaryVBTNodes(void);

	/* While scheduling the tracer, determines wheather or not the image region is 
	 * occupied by a current node
	 * (VBTNode object)
	 */
	int TraceHitTest(VBTNode); 

	/* Compute the directions for secondary nodes for a primary node
	 * (VBTNode ref, direction hist ref)
	 */
	void ComputeSecondaryVBTNodeDirections(VBTNode&, std::vector<double>&);
	
	/* Fit the model using additional constraints for secondary nodes for the given direction vector
	 * (Primary/anchor node ref, Secondary node ref, direction vec)
	 */
	void FitSphereAtVBTNodeSecondary(VBTNode&, VBTNode&, std::vector<double>);

	/* Update the node model for secondary nodes
	 * (VBTNode ref, anchor node ref, iteration number)
	 */
	void UpdateModelSecondary(VBTNode&, VBTNode&, int);

	/* Compute the cost of the trace for the given node
	 * (VBTNode ref)
	 */
	double computeTraceQuality(VBTNode&);

	/* Compute the bin for affinity graph
	 * (direction node)	
	 */
	int ComputeAffinityBin(VBTNode);

	/* Relabel nodes when there is a collision in node IDs
	 * (old ID, new ID)
	 */
	void RelabelForestVBTNodes(int, int);

	/* Check if already labeled nodes are neighbors in the tree
	 * (node 1 ID, node 2 ID)
	 */
	bool CheckNeighbors(int, int);

	/* Sort the queue based on trace quality and return the best trace index and quality
	 */
	void GetBestTrace(std::vector<queue_element>&, queue_element&);
};

class VesselBasedNucleiFeatures{

public:

	VesselBasedNucleiFeatures();
	~VesselBasedNucleiFeatures();

	/* Initialize file locations
	 * (feature table, skeleton image, nuclei label image, vessel segmentation mask)
	 */
	VesselBasedNucleiFeatures(const std::string, const std::string, const std::string, const std::string, const std::string);
	
	/* Read nuclei feature table
	 * (file path)
	 */
	void ReadNucleiFeatureTable(const std::string);

	/* Write the updated nuclei feature table to disk
	 */
	void WriteNucleiFeatureTable(void);
	
	/* Read skeleton image and compute distance map
	 * (file path)
	 */
	void ReadSkeletonImage(const std::string);

	/* Read the nuclei label image
	 * (file path)
	 */
	void ReadNucleiLabelImage(const std::string);

	/* Read the segmentation mask image
	 * (file path)
	 */
	void ReadSegmentationMask(const std::string);

	/* Compute vessel-based features for nuclei
	 */
	void ComputeVesselFeaturesForNuclei(void);

	/* Set the inside region
	 */
	void SetInsideRegionToGlobal(void);

	/* Compute distance map image from the vessel skeleton
	 */
	void ComputeSkeletonDistanceMap(void);	

	/* Compute distance map image from the nuclei label image
	 */
	void ComputeNucleiDistanceMap(void);	

	/* Read node properties from file
	 * (file path)
	 */
	void ReadVBTNodePropertiesFile(const std::string);

private:

	std::string data_folder_path;

    std::vector<vtkSmartPointer<vtkTable> > nucFeaturesTable;

	RenderImageType3D::Pointer skeletonImage;
	RenderImageType3D::Pointer vesselMaskImage;
	LabelImageType3D::Pointer nucLabelImage;
	ImageType3D::Pointer skeletonDistanceMap;
	ImageType3D::Pointer nucDistanceMap;
	RenderImageType3D::Pointer nucBinaryImage;

	RenderImageType3D::RegionType insideRegion;

	std::vector<NucleiObject_VT> nucleiObjects;
	std::vector<VBTNode> skeletonVBTNodes;
};

#endif //_ftkVesselTracer_h
