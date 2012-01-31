#ifndef MICROGLIAREGIONTRACER_H
#define MICROGLIAREGIONTRACER_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodIterator.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkIntTypes.h"

#include <fstream>
#include <cstring>
#include <vector>
#include <cstddef>


#include "Seed.h"				//Simple class to hold seed coordinates
#include "ROIGrabber.h"
#include "LoG.h"
#include "Tree.h"


#include "time.h"

class MicrogliaRegionTracer
{
public:
	typedef fregl_roi::ImageType ImageType;
	typedef LoG::LoGImageType LoGImageType;
	typedef itk::Image<float, 3> VesselnessImageType;

private:
	ImageType::Pointer image;
	std::vector<Seed*> seeds;
	ROIGrabber* roi_grabber;

public:
	MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename);
	~MicrogliaRegionTracer();

	void LoadImage(ImageType::Pointer image);
	void LoadImage(std::string filename);

	void LoadSeedPoints(std::string filename);

	void WriteImage(std::string filename, ImageType::Pointer image);
	void WriteVesselnessImage(std::string filename, VesselnessImageType::Pointer image);
	void WriteSeedImages();
	
	void Trace();

	void CalculateCandidatePixels(Seed* seed, std::vector<ImageType::IndexType> &critical_points_vector);
	void RidgeDetection(std::vector<LoGImageType::Pointer> log_seedimage_vector, ImageType::SizeType size, std::vector<ImageType::IndexType> &critical_points_vector);
	double RunHessian( LoGImageType::Pointer log_image, itk::NeighborhoodIterator<LoGImageType> neighbor_iter);
	double ComputeVesselness( double ev1, double ev2, double ev3 );
	
	void BuildTree(Seed* seed, std::vector<ImageType::IndexType> &critical_points_vector);
	double** BuildAdjacencyGraph(std::vector<ImageType::IndexType> critical_points_vector);
	double MicrogliaRegionTracer::CalculateDistance(itk::uint64_t k, itk::uint64_t l, std::vector<ImageType::IndexType> critical_points_vector);
	Tree* BuildMST(std::vector<ImageType::IndexType> critical_points_vector, double** AdjGraph);
private:

	
};

#endif