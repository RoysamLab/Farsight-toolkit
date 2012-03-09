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

#include <cstring>
#include <vector>
#include <cstddef>
#include <cmath>

#include "Cell.h"				//Simple class to hold seed coordinates
#include "ROIGrabber.h"
#include "LoG.h"

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
	std::string soma_filename;

public:
	MicrogliaRegionTracer::MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename, std::string soma_filename);
	~MicrogliaRegionTracer();
	
	void LoadCellPoints(std::string image_filename);
	
	void Trace();

	void	CalculateCandidatePixels(Cell* cell);
	void	RidgeDetection(Cell* cell);
	void	VesselnessDetection(Cell* cell);
	
	void		BuildTree(Cell* cell);
	double**	BuildAdjacencyGraph(Cell* cell);
	double		CalculateDistance(itk::uint64_t k, itk::uint64_t l, Cell* cell);
	Tree*		BuildMST1(Cell* cell, double** AdjGraph);
	
	void					TraceSkeletonImage(Cell* cell);
	ImageType::IndexType	FindNearestCriticalPointToCentroid(Cell* cell);
	void					WriteTreeToSWCFile(Tree* tree, Cell* cell, std::string filename, std::string filename_local);	
	void					WriteLinkToParent(Node* node, itk::uint64_t tree_depth, Cell* cell, std::ofstream &traceFile, std::ofstream &traceFile_local);

};

#endif