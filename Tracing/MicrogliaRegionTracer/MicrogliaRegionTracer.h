#ifndef MICROGLIAREGIONTRACER_H
#define MICROGLIAREGIONTRACER_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodIterator.h"
#include "itkPathIterator.h"
#include "itkPathConstIterator.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
//#include "itkMinimumMaximumImageCalculator.h"
#include "itkPolyLineParametricPath.h"
#include "itkIntTypes.h"
#include "itkImageDuplicator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "InsightJournalFilters/MinimalPath/itkSpeedFunctionToPathFilter.h"
#include "InsightJournalFilters/MinimalPath/itkImageToPathFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"


#include <cstring>
#include <vector>
#include <cstddef>
#include <cmath>

#include "Cell.h"				//Simple class to hold seed coordinates
#include "ROIGrabber.h"
#include "LoG.h"

#include "Tree.h"
#include "time.h"

#include "itkVector.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkLogImageFilter.h"

#ifdef _OPENMP
	#include "omp.h"
#endif

class MicrogliaRegionTracer
{
private:
	typedef Cell::ImageType						ImageType;
	typedef Cell::LoGImageType					LoGImageType;
	//typedef Cell::HessianImageType			HessianImageType;
	//typedef Cell::HessianTensorType			HessianTensorType;
	typedef Cell::VesselnessImageType			VesselnessImageType;
	typedef Cell::DistanceImageType				DistanceImageType;
	typedef itk::Image< unsigned char, 3 >		MaskedImageType;

	typedef itk::PolyLineParametricPath< 3 >	PathType;

private:
	std::vector<Cell*> cells;
	ROIGrabber* roi_grabber;
	std::string soma_filename;
	double aspect_ratio;

public:
	MicrogliaRegionTracer(const std::string & joint_transforms_filename, const std::string & img_path, const std::string & anchor_filename, const std::string & soma_filename);
	~MicrogliaRegionTracer() {};

	void LoadCellPoints(const std::string & image_filename);

	void	Trace();


	void	CalculateCandidatePixels(Cell* cell);
	void	CreateIsometricImage(Cell* cell);
	void	RidgeDetection(Cell* cell);
	void	VesselnessDetection(Cell* cell);

	void		BuildTree(Cell* cell);
	double**	BuildAdjacencyGraph(Cell* cell);
	double		CalculateDistance(itk::uint64_t k, itk::uint64_t l, Cell* cell);	//THIS IS NOT THE EUCLIDEAN DISTANCE
	Tree*		BuildMST1(Cell* cell, double** AdjGraph);

	void		SmoothTree(Cell* cell, Tree* tree);
	void		SmoothSegments(Cell* cell, Tree* tree, Node* start_node);
	void		SmoothSegments2(Cell* cell, Tree* tree);
	void		SmoothPath(Cell* cell, Tree* tree, Node* start_node, Node* end_node, PathType::Pointer path );

	void		CreateSpeedImage(Cell* cell);

	void		WriteTreeToSWCFile(Tree* tree, Cell* cell, std::string filename, std::string filename_local);	
	void		WriteLinkToParent(Node* node, itk::uint64_t tree_depth, Cell* cell, std::ofstream &traceFile, std::ofstream &traceFile_local);
	
	double		CalculateEuclideanDistance(ImageType::IndexType node1, ImageType::IndexType node2);

};

#endif