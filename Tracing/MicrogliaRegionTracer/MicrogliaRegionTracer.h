#ifndef MicrogliaRegionTracer_H
#define MicrogliaRegionTracer_H

#include "itkIntTypes.h"

#include "itkImage.h"
#include "itkPolyLineParametricPath.h"

#include <cstring>
#include <vector>

#include "Cell.h"
#include "ROIGrabber.h"
#include "Tree.h"

#ifdef _OPENMP
	#include "omp.h"
#endif

class MicrogliaRegionTracer
{
private:
	typedef Cell::ImageType													ImageType;
	typedef Cell::LoGImageType												LoGImageType;
	typedef Cell::VesselnessImageType										VesselnessImageType;
	typedef Cell::DistanceImageType											DistanceImageType;
	typedef itk::Image< unsigned char, 3 >									MaskedImageType;
	typedef Cell::GVFImageType												GVFImageType;

	typedef itk::PolyLineParametricPath< 3 >								PathType;

private:
	std::vector<Cell> cells;
	std::string joint_transforms_filename;
	std::string image_series_pathname;
	std::string anchor_image_filename;
	std::string soma_image_filename;
	double aspect_ratio;
    
    bool aspect_ratio_is_set;
    
    const int itk_default_num_threads; //initialized in contructor
private:
    void                ReadAGroupOfROIs(const int num_cells_in_group, const int group_num, const int num_openmp_threads, ROIGrabber & roi_grabber);
    void                TraceAGroupOfCells(const int num_cells_in_group, const int group_num, const int num_openmp_threads);
	
	void				CalculateSeedPoints(Cell & cell);
	void				RidgeDetectionByLoG(Cell & cell);
	void				SeedDetectionByVesselness(Cell & cell);
	void				SeedAdjustment(Cell & cell, const VesselnessImageType::Pointer & thresholded_vesselness_img);

	void				BuildTree(Cell & cell);
	float**				BuildAdjacencyGraph(Cell & cell);
	float				CalculateDistance(Cell & cell, itk::uint64_t k, itk::uint64_t l);
	Tree*				BuildMST1(Cell & cell, float** AdjGraph);

	void				SmoothTree(Cell & cell, Tree* smoothed_tree);
	void				SmoothSegments(Cell & cell, Tree* smoothed_tree, Node* start_node);
	void				ReplaceTreeSegmentWithPath(Cell & cell, Tree* smoothed_tree, PathType::Pointer speed_path, Node* start_node, Node* end_node);
	PathType::Pointer	SmoothPath(Cell & cell, Tree* smoothed_tree, Node* start_node, Node* end_node, PathType::Pointer path );

	void				PruneTree(Cell & cell, Tree* smoothed_tree);
	void				CheckSegment(Cell & cell, Node* start_node, Tree * smoothed_tree);
	bool				ShouldPruneSegment(Cell & cell, PathType::Pointer path);

	void				CreateSpeedImage(Cell & cell);
    
public:
                        explicit MicrogliaRegionTracer();
                        ~MicrogliaRegionTracer();
    
    void				Trace();
    
    void				SetJointTransformsFile(const std::string & joint_transforms_filename);
	void				SetImageSeriesPath(const std::string & image_series_pathname);
	void				SetAnchorImage(const std::string & anchor_image_filename);
	void				LoadSeedPoints(const std::string & seedpoints_filename);
	void				SetSomaImage(const std::string & soma_image_filename);
	void				SetAspectRatio(const float aspect_ratio);
};

#endif