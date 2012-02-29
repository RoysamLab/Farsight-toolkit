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
#include "itkMinimumMaximumImageCalculator.h"
#include "itkPolyLineParametricPath.h"
#include "itkIntTypes.h"

#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkShiftScaleImageFilter.h"

#include <cstring>
#include <vector>
#include <cstddef>
#include <cmath>

#include "Cell.h"				//Simple class to hold seed coordinates
#include "ROIGrabber.h"
#include "LoG.h"
//#include "Vesselness.h"

#include "Tree.h"
#include "time.h"

#ifdef _OPENMP
	#include "omp.h"
#endif

class MicrogliaRegionTracer
{
private:
	typedef Cell::ImageType						ImageType;
	typedef Cell::LoGImageType					LoGImageType;
	//typedef Cell::HessianImageType				HessianImageType;
	//typedef Cell::HessianTensorType				HessianTensorType;
	typedef Cell::VesselnessImageType			VesselnessImageType;
	typedef Cell::DistanceImageType				DistanceImageType;
	typedef itk::Image< float, 3 >				VoronoiImageType;
	typedef itk::Image< unsigned char, 3 >		MaskedImageType;

private:
	std::vector<Cell*> cells;
	ROIGrabber* roi_grabber;

public:
	MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename);
	~MicrogliaRegionTracer();
	
	void LoadCellPoints(std::string image_filename, std::string soma_filename);
	
	void Trace();

	void	CalculateCandidatePixels(Cell* cell);
	void	RidgeDetection(Cell* cell);
	void	RidgeDetection2(Cell* cell);
	void	VesselnessDetection(Cell* cell);
	double	RunHessian( LoGImageType::Pointer log_image, itk::NeighborhoodIterator<LoGImageType> neighbor_iter, double max_intensity_multiscale);
	double	ComputeVesselness( double ev1, double ev2, double ev3, double maximum_intensity );

	void		BuildTree(Cell* cell);
	double**	BuildAdjacencyGraph(Cell* cell);
	double		CalculateDistance(itk::uint64_t k, itk::uint64_t l, Cell* cell);
	Tree*		BuildMST1(Cell* cell, double** AdjGraph);
	
	void					TraceSkeletonImage(Cell* cell);
	ImageType::IndexType	FindNearestCriticalPointToCentroid(Cell* cell);
	void					Trace2();
	void					FollowSkeleton(Cell* cell, ImageType::IndexType index, ImageType::Pointer visited_image, itk::int64_t parent_id, itk::int64_t &swc_line_number, std::ofstream &traceFile, std::ofstream &traceFile_local);

};

#endif