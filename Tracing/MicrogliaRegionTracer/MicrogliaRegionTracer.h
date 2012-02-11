#ifndef MICROGLIAREGIONTRACER_H
#define MICROGLIAREGIONTRACER_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodIterator.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkIntTypes.h"
#include "itkMaskNegatedImageFilter.h"

#include <fstream>
#include <cstring>
#include <vector>
#include <cstddef>

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
	typedef Cell::ImageType ImageType;
	typedef Cell::LoGImageType LoGImageType;
	typedef Cell::VesselnessImageType VesselnessImageType;
	typedef itk::Image<unsigned char, 3> MaskedImageType;

private:
	std::vector<Cell*> cells;
	ROIGrabber* roi_grabber;

public:
	MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename);
	~MicrogliaRegionTracer();

	ImageType::Pointer GetMaskedImage(MaskedImageType::Pointer mask, ImageType::Pointer image);
	ImageType::Pointer GetMaskedImage(std::string filename, ImageType::Pointer image);
	
	void LoadCellPoints(std::string image_filename, std::string soma_filename);

	void WriteImage(std::string filename, ImageType::Pointer image);
	void WriteVesselnessImage(std::string filename, VesselnessImageType::Pointer image);
	
	void Trace();

	void CalculateCandidatePixels(Cell* cell);
	void RidgeDetection(Cell* cell);
	double RunHessian( LoGImageType::Pointer log_image, itk::NeighborhoodIterator<LoGImageType> neighbor_iter, double max_intensity_multiscale);
	double ComputeVesselness( double ev1, double ev2, double ev3, double maximum_intensity );
	
	void BuildTree(Cell* cell);
	double** BuildAdjacencyGraph(Cell* cell);
	double CalculateDistance(itk::uint64_t k, itk::uint64_t l, Cell* cell);
	Tree* BuildMST1(Cell* cell, double** AdjGraph);
	Tree* BuildMST2(Cell* cell, double** AdjGraph);
	
};

#endif