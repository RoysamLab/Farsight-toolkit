#include "MicrogliaRegionTracer.h"

#include "itkImageDuplicator.h"
#include "InsightJournalFilters/MinimalPath/itkSpeedFunctionToPathFilter.h"
#include "InsightJournalFilters/MinimalPath/itkImageToPathFilter.h"

#include "itkNeighborhoodIterator.h"
#include "itkPathIterator.h"
#include "itkPathConstIterator.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorRescaleIntensityImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "BoundingBoxFromTraceCalculator.h"

#include <cmath>
#include <algorithm>

#define MASK 0
#define PI 3.1415926535897932384626433832795

//typedefs
typedef struct
{
	double x;
	double y;
	double z;
	double saliency;
} VesselnessPointType;

//Non-member function prototypes
float CalculateEuclideanDistance(Cell::ImageType::IndexType node1, Cell::ImageType::IndexType node2);
bool VesselnessPointTypeComparator(VesselnessPointType first_pt, VesselnessPointType second_pt);

MicrogliaRegionTracer::MicrogliaRegionTracer() :
    aspect_ratio_is_set(false),
    itk_default_num_threads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads())
{
#ifdef _OPENMP
	std::cerr << "OpenMP detected!" << std::endl;
#endif
}

MicrogliaRegionTracer::~MicrogliaRegionTracer()
{
}

void MicrogliaRegionTracer::SetJointTransformsFile(const std::string & joint_transforms_filename)
{
	this->joint_transforms_filename = joint_transforms_filename;
}

void MicrogliaRegionTracer::SetImageSeriesPath(const std::string & image_series_pathname)
{
	this->image_series_pathname = image_series_pathname;
}

void MicrogliaRegionTracer::SetAnchorImage(const std::string & anchor_image_filename)
{
	this->anchor_image_filename = anchor_image_filename;
}

void MicrogliaRegionTracer::LoadSeedPoints(const std::string & seedpoints_filename)
{
	if (!aspect_ratio_is_set)
        throw std::runtime_error("The aspect ratio must be set before attemping to load the seed points");
    
    std::cout << "Opening centroid file named: " << seedpoints_filename << std::endl;

	std::ifstream seed_point_file;
	seed_point_file.open(seedpoints_filename.c_str());

	float cellX, cellY, cellZ;
	itk::uint64_t num_seed_points = 0;
	if (seed_point_file.is_open())
	{
		while (seed_point_file >> cellX >> cellY >> cellZ)
		{	
			Cell cell(cellX, cellY, cellZ, this->aspect_ratio);
			this->cells.push_back(cell);
			++num_seed_points;
		}

		if (num_seed_points == 0)
			throw std::runtime_error("Error, you probably didn't supply a correct seed point file");

		std::cout << "Number of seed points read: " << num_seed_points << std::endl;
	}
	else
	{
		std::cerr << "Error opening centroid file: " << seedpoints_filename << std::endl;
	}	
}

void MicrogliaRegionTracer::SetSomaImage(const std::string & soma_image_filename)
{
	this->soma_image_filename = soma_image_filename;
}

void MicrogliaRegionTracer::SetAspectRatio(const float aspect_ratio)
{	
    this->aspect_ratio = aspect_ratio;
    this->aspect_ratio_is_set = true;
}

/* This is the main loop where tracing takes place */
void MicrogliaRegionTracer::Trace()
{
	//Create groups of n cells at a time (where n is the number of logical processors)
    //This is needed because OpenMP doesn't play well with ITK
#ifdef _OPENMP
    int num_openmp_threads = omp_get_num_procs();
#else
    int num_openmp_threads = 1;    //Default 1 thread in a group
#endif
    
    int num_groups = std::ceil(cells.size() / (float) num_openmp_threads);
	std::cerr << "Number of groups of cells to process: " << num_groups << std::endl;
    std::cerr << "Number of cells in each group (except the last group): " << num_openmp_threads << std::endl;
    std::cerr << "Number of cells to process in the last group: " << cells.size() % num_openmp_threads << std::endl;
    
    ROIGrabber roi_grabber(this->joint_transforms_filename, this->image_series_pathname, this->anchor_image_filename);

    //This for loop are to handle the case of parallelizing the processing of the
    //microglia because the fregl_roi class does not like to be accessed by different OpenMP
    //threads even though if it guarded by a omp critical section.
    //First we read num_thread number of cells sequentially, then we process the group in parallel.
    for (int group_num = 0; group_num < num_groups; group_num++)
    {
        
        int num_cells_in_group = (cells.size() / ((group_num+1) * num_openmp_threads)) ? num_openmp_threads : cells.size() % num_openmp_threads; //sets num_local_threads to the number of threads if we are not on the last group, otherwise set it to the number of remaining threads

        ReadAGroupOfROIs(num_cells_in_group, group_num, num_openmp_threads, roi_grabber);  //Read a bunch of cells WITHOUT OpenMP        
        TraceAGroupOfCells(num_cells_in_group, group_num, num_openmp_threads);             //Then we trace a group of cells using OpenMP
    }
}

void MicrogliaRegionTracer::ReadAGroupOfROIs(const int num_cells_in_group, const int group_num, const int num_openmp_threads, ROIGrabber & roi_grabber)
{
    
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads( itk_default_num_threads );	//Change default number of threads to itk_default_num_threads for ITK filters in fregl_roi
    
    //In this for loop we read num_thread ROIs without OpenMP (see comment before outer for loop)
    for (int local_thread_num = 0; local_thread_num < num_cells_in_group; local_thread_num++)
    {
        int global_thread_num = group_num * num_openmp_threads + local_thread_num;
        
        Cell & cell = cells[global_thread_num];
        std::cerr << "Reading cell #: " << global_thread_num << std::endl;
        
        //Setup the size of the initial ROI
        ImageType::SizeType roi_size;
        roi_size[0] = 200;
        roi_size[1] = 200;
        roi_size[2] = 100;
        
        cell.SetRequestedSize(roi_size);
        
        //Grab the initial cellimage
        //std::cout << "Grabbing ROI for cell" << std::endl;
        ImageType::IndexType shift_index;
        ImageType::Pointer temp_cell_image = roi_grabber.GetROI(cell, roi_size, shift_index);
        
        //Set the origin of the image to be {0, 0, 0}
        ImageType::PointType origin = temp_cell_image->GetOrigin();
        cell.SetOrigin(origin);
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        temp_cell_image->SetOrigin(origin);
        
        //Duplicate the image so that when roi_grabber gets destroyed the roi's dont get destroyed with it
        typedef itk::ImageDuplicator< ImageType > DuplicatorType;
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(temp_cell_image);
        ftk::TimeStampOverflowSafeUpdate( duplicator.GetPointer() );
        cell.image = duplicator->GetOutput();
        cell.image->DisconnectPipeline();
        
        //Shift_index is just to demonstrate how much the center point of the local point of the image if it got cropped by the left, top, closest point of the image
        cell.SetShiftIndex(shift_index);
        //std::cout << cell.GetShiftIndex()[0] << " " << cell.GetShiftIndex()[1] << " " << cell.GetShiftIndex()[2] << std::endl;
        
        //Make the file name of the raw cell image
        std::stringstream cell_filename_stream;
        cell_filename_stream << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << ".TIF";	//X_Y_Z.TIF
        
        //Write the cell image
        Cell::WriteImage(cell_filename_stream.str(), cell.image);
        
        roi_size = cell.image->GetLargestPossibleRegion().GetSize();	//The size of the returned image may not be the size of the image that entered because clipping at edges
        cell.SetSize(roi_size);
    }
}

void MicrogliaRegionTracer::TraceAGroupOfCells(const int num_cells_in_group, const int group_num, const int num_openmp_threads)
{
#ifdef _OPENMP
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads( std::max<float>(1, num_openmp_threads / cells.size()));	//We trace multiple cells with OpenMP, so reduce the number of ITK threads so that we don't cause thread contention#
#endif
    
	std::cout << "Tracing a group of cells using " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() << " itk threads and " << num_cells_in_group << " openmp threads" << std::endl;

    #pragma omp parallel for
    for (int local_thread_num = 0; local_thread_num < num_cells_in_group; local_thread_num++)
    {
        int global_thread_num = group_num * num_openmp_threads + local_thread_num;
        Cell & cell = cells[global_thread_num];
        
        //Get the mask
#if MASK	//development only
        cell.GetMask(this->soma_filename);
#endif
        std::cout << "Calculating candidate pixels for a new cell" << std::endl;
        CalculateSeedPoints(cell);
        
        std::cout << "Detected " << cell.critical_points_queue.size() << " critical points" << std::endl;
        
        std::cout << "Tree Building" << std::endl;
        BuildTree(cell);
    }
}

/* This function determines the candidate pixels (pixels which we connect to form the tree) */
void MicrogliaRegionTracer::CalculateSeedPoints(Cell & cell)
{
	cell.CreateIsotropicImage();	//Must be called before CreateVesselnessImage and CreateLoGImage
	cell.CreateVesselnessImage();
	
	cell.CreateLoGImage();
	//std::cout << "Ridge Detection By LoG" << std::endl;
	//RidgeDetectionByLoG(cell);

	SeedDetectionByVesselness(cell);

	std::ostringstream critical_points_filename_stream;
	critical_points_filename_stream << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << "_critical.nrrd";
	Cell::WriteImage(critical_points_filename_stream.str(), cell.critical_point_image);
}

/* This function does Ridge detection (Laplacian of Gaussian implementation) */
void MicrogliaRegionTracer::RidgeDetectionByLoG( Cell & cell )
{
	//Add the seed to the critical points
	itk::Index<3> seed_index = {{	cell.GetRequestedSize()[0]/2 + cell.GetShiftIndex()[0], 
									cell.GetRequestedSize()[1]/2 + cell.GetShiftIndex()[1], 
									cell.GetRequestedSize()[2]/2 + cell.GetShiftIndex()[2] }};

	cell.critical_points_queue.push_front(seed_index);
	
	//Make a new image to store the ridge pixels
	ImageType::SizeType size = cell.multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	cell.ridge_image = LoGImageType::New();
	ImageType::IndexType start;
	start.Fill(0);
	ImageType::RegionType region(start, size);
	cell.ridge_image->SetRegions(region);
	cell.ridge_image->Allocate();
	cell.ridge_image->FillBuffer(0);

    //Non-maximal suppression
	//Make a iterator for the image and a neighborhood around the current point we are visiting
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator< LoGImageType > ridge_img_iter(cell.ridge_image, cell.ridge_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > LoG_neighbor_iter(rad, cell.multiscale_LoG_image, cell.multiscale_LoG_image->GetLargestPossibleRegion());

	while(!LoG_neighbor_iter.IsAtEnd()) 
	{
        LoGImageType::PixelType a = LoG_neighbor_iter.GetPixel(4);
		
		if (a > 0.01) //Ignore pixels less than 1.0 % response
		{
			bool isInBounds;
			LoGImageType::PixelType a1 = LoG_neighbor_iter.GetPixel(0, isInBounds);
			if (!isInBounds)
				a1 = 0;
			LoGImageType::PixelType a2 = LoG_neighbor_iter.GetPixel(1, isInBounds);
			if (!isInBounds)
				a2 = 0;
			LoGImageType::PixelType a3 = LoG_neighbor_iter.GetPixel(2, isInBounds);
			if (!isInBounds)
				a3 = 0;
			LoGImageType::PixelType a4 = LoG_neighbor_iter.GetPixel(3, isInBounds);
			if (!isInBounds)
				a4 = 0;
			LoGImageType::PixelType a5 = LoG_neighbor_iter.GetPixel(5, isInBounds);	//note the change in indexing
			if (!isInBounds)
				a5 = 0;
			LoGImageType::PixelType a6 = LoG_neighbor_iter.GetPixel(6, isInBounds);
			if (!isInBounds)
				a6 = 0;
			LoGImageType::PixelType a7 = LoG_neighbor_iter.GetPixel(7, isInBounds);
			if (!isInBounds)
				a7 = 0;
			LoGImageType::PixelType a8 = LoG_neighbor_iter.GetPixel(8, isInBounds);
			if (!isInBounds)
				a8 = 0;
        
			bool ridge_pixel = false;
        
			if ((a1 + a4 + a6 < a2 + a + a7 && a2 + a + a7 > a3 + a5 + a8) ||   //pixel on vertical ridge
				(a1 + a2 + a3 < a4 + a + a5 && a4 + a + a5 > a6 + a7 + a8) ||   //pixel on horizontal ridge
				(a1 + a2 + a4 < a3 + a + a6 && a3 + a + a6 > a5 + a7 + a8) ||   //pixel on 45 degree ridge
				(a2 + a3 + a5 < a1 + a + a8 && a1 + a + a8 > a4 + a6 + a7)      //pixel on 135 degree ridge
				)
				ridge_pixel = true;
        
			if (ridge_pixel)
				ridge_img_iter.Set(a);
		}

		++ridge_img_iter;
		++LoG_neighbor_iter;
	}


	//Make a new image to store the critical points	
	ImageType::SizeType critical_point_image_size = cell.multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	cell.critical_point_image = ImageType::New();
	ImageType::IndexType critical_point_image_start;
	critical_point_image_start.Fill(0);
	ImageType::RegionType critical_image_region(critical_point_image_start, critical_point_image_size);
	cell.critical_point_image->SetRegions(critical_image_region);
	cell.critical_point_image->Allocate();
	cell.critical_point_image->FillBuffer(0);

	//Make a iterator for the image and a neighborhood around the current point we are visiting
	itk::Size<3> ridge_neighbor_iter_rad = {{3,3,3}};
	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell.critical_point_image, cell.critical_point_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > ridge_neighbor_iter(ridge_neighbor_iter_rad, cell.ridge_image, cell.ridge_image->GetLargestPossibleRegion());

	while(!ridge_neighbor_iter.IsAtEnd()) 
	{
		unsigned int neighborhood_size = (ridge_neighbor_iter_rad[0] * 2 + 1) * (ridge_neighbor_iter_rad[1] * 2 + 1) * (ridge_neighbor_iter_rad[2] * 2 + 1);
		unsigned int center_pixel_offset_index = neighborhood_size / 2;

		LoGImageType::PixelType center_pixel_intensity = ridge_neighbor_iter.GetPixel(center_pixel_offset_index);
		
		if (center_pixel_intensity >= 0.005)	//Must have greater than 0.5% response from LoG to even be considered for local maximum (so we don't choose points that are part of noise)
		{	
			bool local_maximum = true;
			
			for (unsigned int neighborhood_index = 0; neighborhood_index < neighborhood_size; ++neighborhood_index)
			{
				if (neighborhood_index != center_pixel_offset_index)
				{
					bool isInBounds; //true if the pixel is in bounds
					LoGImageType::PixelType neighbor_pixel_log_intensity = ridge_neighbor_iter.GetPixel(neighborhood_index, isInBounds);

					if (isInBounds && neighbor_pixel_log_intensity > center_pixel_intensity)
						local_maximum = false;
				}
			}

			if (local_maximum)
			{
				critical_point_img_iter.Set(255);
				cell.critical_points_queue.push_back(critical_point_img_iter.GetIndex());
			}
		}

		++critical_point_img_iter;
		++ridge_neighbor_iter;
	}
}

void MicrogliaRegionTracer::SeedDetectionByVesselness(Cell & cell)
{
	//Take the log() of the vesselness
	typedef itk::LogImageFilter< VesselnessImageType, VesselnessImageType > LogImageFilterType;
	LogImageFilterType::Pointer log_image_filter = LogImageFilterType::New();
	log_image_filter->SetInput(cell.vesselness_image);

	try
	{
		log_image_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Exception in log_image_filter in SeedDetectionByVesselness(): " << err << std::endl;
	}

	VesselnessImageType::Pointer log_vesselness_img = log_image_filter->GetOutput();
	
	Cell::WriteImage("log_vesselness.nrrd", log_vesselness_img);

	//Clamp the pixel values
	itk::ImageRegionIterator< VesselnessImageType > log_vesselness_iter(log_vesselness_img, log_vesselness_img->GetLargestPossibleRegion());

	while (!log_vesselness_iter.IsAtEnd())
	{
		VesselnessImageType::PixelType value = log_vesselness_iter.Get();

		if ( value <= -20.0 ) //pixel value is less than -20
			log_vesselness_iter.Set( -20.0 ); //HACK!: Set it to -20.0 so that we can actually do image processing without running into FP maximum range problems

		++log_vesselness_iter;
	}

	Cell::WriteImage("clamped_log_vesselness.nrrd", log_vesselness_img);
	
	//Threshold the image to get rid of BG pixels
	typedef itk::OtsuMultipleThresholdsImageFilter< VesselnessImageType, ImageType > ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetInput(log_vesselness_img);
	threshold_filter->SetNumberOfThresholds(2);
	threshold_filter->SetNumberOfHistogramBins(1024);
	try
	{
		threshold_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Exception in threshold_filter in SeedDetectionByVesselness(): " << err << std::endl;
	}

	cell.threshold_label_mask = threshold_filter->GetOutput();
	Cell::WriteImage("threshold_label_mask.nrrd", cell.threshold_label_mask);

	
	//Prune the label image to only take the pixels with label "2" (foreground)
	//Make an image to hold the foreground pixels
	ImageType::Pointer foreground_vesselness_mask = ImageType::New();
	foreground_vesselness_mask->SetRegions(cell.threshold_label_mask->GetLargestPossibleRegion());
	foreground_vesselness_mask->Allocate();
	foreground_vesselness_mask->FillBuffer(0);

	itk::ImageRegionConstIterator< ImageType > threshold_label_img_iter(cell.threshold_label_mask, cell.threshold_label_mask->GetLargestPossibleRegion());
	itk::ImageRegionIterator< ImageType > foreground_vesselness_mask_iter(foreground_vesselness_mask, foreground_vesselness_mask->GetLargestPossibleRegion());

	while (!threshold_label_img_iter.IsAtEnd())
	{
		ImageType::PixelType label = threshold_label_img_iter.Get();

		if (label == 2)
			foreground_vesselness_mask_iter.Set(255);
		
		++threshold_label_img_iter;
		++foreground_vesselness_mask_iter;
	}

	Cell::WriteImage("foreground_vesselness_mask.nrrd", foreground_vesselness_mask);	

	//Use the mask to get the thresholded vesselness image
	typedef itk::MaskImageFilter< VesselnessImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer mask_image_filter = MaskImageFilterType::New();
	mask_image_filter->SetInput(cell.vesselness_image);
	mask_image_filter->SetMaskImage(foreground_vesselness_mask);

	try
	{
		mask_image_filter->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Exception in mask_image_filter in SeedDetectionByVesselness(): " << err << std::endl;
	}

	VesselnessImageType::Pointer thresholded_vesselness_image = mask_image_filter->GetOutput();

	Cell::WriteImage("thresholded_vesselness_image.nrrd", thresholded_vesselness_image);

	//Adjust the seed points
	SeedAdjustment(cell, thresholded_vesselness_image);

	//Non-maximal suppression
	VesselnessImageType::Pointer non_maximal_vesselness_img = VesselnessImageType::New();
	non_maximal_vesselness_img->SetRegions(cell.vesselness_image->GetLargestPossibleRegion());
	non_maximal_vesselness_img->Allocate();
	non_maximal_vesselness_img->FillBuffer(0);

	cell.critical_point_image = ImageType::New();
	cell.critical_point_image->SetRegions(cell.vesselness_image->GetLargestPossibleRegion());
	cell.critical_point_image->Allocate();
	cell.critical_point_image->FillBuffer(0);
		
	itk::Size< 3 > rad = {{1,1,1}};
	itk::ConstNeighborhoodIterator< VesselnessImageType > vesselness_img_iter(rad, cell.adjusted_seed_img, cell.adjusted_seed_img->GetLargestPossibleRegion());
	itk::ImageRegionIterator< VesselnessImageType> non_maximal_vesselness_img_iter(non_maximal_vesselness_img, non_maximal_vesselness_img->GetLargestPossibleRegion());
	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell.critical_point_image, cell.critical_point_image->GetLargestPossibleRegion());
		
	//Make a iterator for the image and a neighborhood around the current point we are visiting
	while (!vesselness_img_iter.IsAtEnd())
	{
		unsigned int neighborhood_size = (rad[0] * 2 + 1) * (rad[1] * 2 + 1) * (rad[2] * 2 + 1);
		unsigned int center_pixel_offset_index = neighborhood_size / 2;

		VesselnessImageType::PixelType center_pixel_intensity = vesselness_img_iter.GetPixel(center_pixel_offset_index);

		bool local_maximum = true;
		for (unsigned int neighborhood_index = 0; neighborhood_index < neighborhood_size; ++neighborhood_index)
		{
			if (neighborhood_index != center_pixel_offset_index)
			{
				bool isInBounds; //true if the pixel is in bounds
				VesselnessImageType::PixelType neighbor_pixel_log_intensity = vesselness_img_iter.GetPixel(neighborhood_index, isInBounds);

				if (isInBounds && neighbor_pixel_log_intensity >= center_pixel_intensity)
					local_maximum = false;
			}
		}

		if (local_maximum)
		{
			non_maximal_vesselness_img_iter.Set(center_pixel_intensity);
			critical_point_img_iter.Set(255);
			cell.critical_points_queue.push_back(non_maximal_vesselness_img_iter.GetIndex());
		}
		
		++critical_point_img_iter;
		++vesselness_img_iter;
		++non_maximal_vesselness_img_iter;
	}
	Cell::WriteImage("non_maximal.nrrd", non_maximal_vesselness_img);
}

bool VesselnessPointTypeComparator(VesselnessPointType first_pt, VesselnessPointType second_pt)
{
	return (first_pt.saliency > second_pt.saliency);	
}

void MicrogliaRegionTracer::SeedAdjustment(Cell & cell, const VesselnessImageType::Pointer & thresholded_seeds_img)
{
	std::cout << "Entering seed adjustment" << std::endl;
	std::vector< VesselnessPointType > seed_points;
	
	//Extract the seeds from the thresholded_seeds_img and put them into a vector
	itk::ImageRegionConstIterator< VesselnessImageType > thresholded_seeds_img_iter(thresholded_seeds_img, thresholded_seeds_img->GetLargestPossibleRegion());
	while (!thresholded_seeds_img_iter.IsAtEnd())
	{
		VesselnessImageType::IndexType index = thresholded_seeds_img_iter.GetIndex();
		VesselnessImageType::PixelType saliency = thresholded_seeds_img_iter.Get();

		if (saliency != 0.0)
		{
			VesselnessPointType seed;
			seed.x = index[0];
			seed.y = index[1];
			seed.z = index[2];
			seed.saliency = saliency;

			seed_points.push_back(seed);
		}
		++thresholded_seeds_img_iter;
	}
	std::cout << "Extracted " << seed_points.size() << " seeds from the thresholded_seeds_img" << std::endl;

	////Sort the vector of seeds by their saliency, most salient at the front of the vector
	//std::sort(seed_points.begin(), seed_points.end(), VesselnessPointTypeComparator);	

	typedef itk::VectorRescaleIntensityImageFilter< GVFImageType, GVFImageType > VectorRescaleFilterType;
	VectorRescaleFilterType::Pointer vector_rescale_filter = VectorRescaleFilterType::New();
	vector_rescale_filter->SetInput(cell.gvf_image);
	vector_rescale_filter->SetOutputMaximumMagnitude(1.0);	

	try
	{
		vector_rescale_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "Error in vector_rescale_filter in SeedAdjustment(): " << err << std::endl;
	}
	
	//Create an interpolator for the GVF field
	typedef itk::VectorLinearInterpolateImageFunction< GVFImageType, double > VectorInterpolatorType;
	VectorInterpolatorType::Pointer vector_interpolator = VectorInterpolatorType::New();
	vector_interpolator->SetInputImage(vector_rescale_filter->GetOutput());

	//Create an interpolator for the vesselness_image to get the new saliencies
	typedef itk::LinearInterpolateImageFunction< VesselnessImageType, double > InterpolatorType;
	InterpolatorType::Pointer vesselness_interpolator = InterpolatorType::New();
	vesselness_interpolator->SetInputImage(cell.vesselness_image);
	
	//Actually move the seeds around
	std::cout << "Moving GVF seeds around" << std::endl;
	unsigned int num_iterations = 1000;
	unsigned int iteration_num = 0;
	double learning_rate = 0.1;

	VectorInterpolatorType::IndexType end_index = vector_interpolator->GetEndIndex();

	while (iteration_num < num_iterations)
	{
		for (std::vector< VesselnessPointType >::iterator seed_points_iter = seed_points.begin(); seed_points_iter != seed_points.end(); ++seed_points_iter)
		{
			VesselnessPointType & seed = *seed_points_iter;
			
			VectorInterpolatorType::ContinuousIndexType gvf_index;
			gvf_index[0] = seed.x;
			gvf_index[1] = seed.y;
			gvf_index[2] = seed.z;

			GVFImageType::PixelType gradient_vector = vector_interpolator->EvaluateAtContinuousIndex(gvf_index);
			
			seed.x += learning_rate * gradient_vector[0];
			seed.y += learning_rate * gradient_vector[1];
			seed.z += learning_rate * gradient_vector[2];

			//Make sure the seeds are in-bounds
			seed.x = std::min< double >(seed.x, end_index[0]);
			seed.y = std::min< double >(seed.y, end_index[1]);
			seed.z = std::min< double >(seed.z, end_index[2]);

			seed.x = std::max< double >(seed.x, 0.0);
			seed.y = std::max< double >(seed.y, 0.0);
			seed.z = std::max< double >(seed.z, 0.0);

			InterpolatorType::ContinuousIndexType vesselness_index;
			vesselness_index[0] = seed.x;
			vesselness_index[1] = seed.y;
			vesselness_index[2] = seed.z;

			
			seed.saliency = vesselness_interpolator->EvaluateAtContinuousIndex(vesselness_index);
		}
		++iteration_num;
	}
	std::cout << "Done moving the seeds around" << std::endl;

	//Write the moved seeds back out to an image
	cell.adjusted_seed_img = VesselnessImageType::New();
	cell.adjusted_seed_img->SetRegions(cell.vesselness_image->GetLargestPossibleRegion());
	cell.adjusted_seed_img->Allocate();
	cell.adjusted_seed_img->FillBuffer(0);

	std::vector< VesselnessPointType >::iterator adjusted_seeds_iter;

	for (std::vector< VesselnessPointType >::iterator seed_points_iter = seed_points.begin(); seed_points_iter != seed_points.end(); ++seed_points_iter)
	{
		const VesselnessPointType seed = *seed_points_iter;
		
		VectorInterpolatorType::ContinuousIndexType seed_index;
		seed_index[0] = seed.x;
		seed_index[1] = seed.y;
		seed_index[2] = seed.z;

		VesselnessImageType::IndexType converted_index;
		vector_interpolator->ConvertContinuousIndexToNearestIndex(seed_index, converted_index);

		double saliency = cell.adjusted_seed_img->GetPixel(converted_index); //Get the vesselness for the converted pixel

		if (seed.saliency > saliency)
			cell.adjusted_seed_img->SetPixel(converted_index, seed.saliency);

	}

	Cell::WriteImage("adjusted_seed_img.nrrd", cell.adjusted_seed_img);
}

/* After the candidates pixels are calculated, this function connects all the candidate pixels into a minimum spanning tree based on a cost function */
void MicrogliaRegionTracer::BuildTree(Cell & cell)
{
	std::cout << "Building Adjacency Graph" << std::endl;
	float** AdjGraph = BuildAdjacencyGraph(cell);

	std::cout << "Building Minimum Spanning Tree" << std::endl;
	Tree* tree = BuildMST1(cell, AdjGraph);

	std::ostringstream swc_filename_stream, swc_filename_stream_local;
	swc_filename_stream << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << "_tree.swc";
	swc_filename_stream_local << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << "_tree_local.swc";

	cell.WriteTreeToSWCFile(tree, swc_filename_stream.str(), swc_filename_stream_local.str());

	Tree * smoothed_tree = new Tree(*tree); //Copy the original tree into smoothed_tree
	
	cell.CreateSpeedImage();

	for (int k = 0; k < 1; ++k)
	{
		//Smooth the tree via Fast Marching
		SmoothTree(cell, smoothed_tree);
	
		//Prune the leaves of the tree if there are too many blank pixels (according to threshold_mask)
		PruneTree(cell, smoothed_tree);
	}

	//Smooth the tree one final time
	SmoothTree(cell, smoothed_tree);
	
	std::ostringstream swc_filename_stream_smoothed, swc_filename_stream_local_smoothed;
	swc_filename_stream_smoothed << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << "_tree_smoothed.swc";
	swc_filename_stream_local_smoothed << cell.getX() << "_" << cell.getY() << "_" << cell.getZ() << "_tree_local_smoothed.swc";

	cell.WriteTreeToSWCFile(smoothed_tree, swc_filename_stream_smoothed.str(), swc_filename_stream_local_smoothed.str());
    
    BoundingBox bounding_box = BoundingBoxFromTraceCalculator::GetBoundingBoxFromTree(*smoothed_tree);
    std::cout << bounding_box << std::endl;

	for (int k = 0; k < cell.critical_points_queue.size(); k++)
		delete[] AdjGraph[k];
	delete[] AdjGraph;

	delete tree;
	delete smoothed_tree;
}

/* This function generates a matrix of all possible pairs of candidate pixels */
float** MicrogliaRegionTracer::BuildAdjacencyGraph(Cell & cell)
{
	float** AdjGraph = new float*[cell.critical_points_queue.size()];
	for (int k = 0; k < cell.critical_points_queue.size(); k++)
		AdjGraph[k] = new float[cell.critical_points_queue.size()];

	for (itk::int64_t k = 0; k < (itk::int64_t)cell.critical_points_queue.size(); k++)
	{
		//std::cout << "Calculating distance for node " << k << std::endl;
		for (itk::uint64_t l = 0; l < cell.critical_points_queue.size(); l++)
		{
			AdjGraph[k][l] = CalculateDistance(cell, k, l);
		}
	}

	std::cout << "Done calculating distance graph" << std::endl;
	return AdjGraph;
}

/* This is the distance part of the cost function */
float MicrogliaRegionTracer::CalculateDistance(Cell & cell, itk::uint64_t node_from, itk::uint64_t node_to)
{
	ImageType::IndexType node1 = cell.critical_points_queue[node_from];
	ImageType::IndexType node2 = cell.critical_points_queue[node_to];

	ImageType::IndexType trace_vector;
	trace_vector[0] = node2[0] - node1[0];
	trace_vector[1] = node2[1] - node1[1];
	trace_vector[2] = node2[2] - node1[2];

	float mag_trace_vector = sqrt(pow(trace_vector[0], 2.0f) + pow(trace_vector[1], 2.0f) + pow(trace_vector[2], 2.0f));
	
	if ((node_from == node_to) || (mag_trace_vector < 0.0001f))	//Very likely that these two points are equal so we shouldn't try to connect them
		return std::numeric_limits< float >::max();

#if MASK
	//If we are connecting from the root node
	if (node_from == 0)
	{	
		if (cell.soma_label_image->GetPixel(node1) == cell.soma_label_image->GetPixel(node2)) //If it is inside the soma our centroid belongs to, we give it 0 weight
			return 0;
		else if (cell.soma_label_image->GetPixel(node2) != 0)											//We are trying to connect to a pixel in another soma, so we give it infinite weight
			return std::numeric_limits< float >::max();	
	}
#endif

	////Filter out the path between the pair of nodes if there are too many low-valued pixels in between
	//PathType::Pointer path = PathType::New();
	//path->Initialize();

	////std::cout << "Start Index: " << node1 << " " << "End Index: " << node2 << std::endl;
	//
	//path->AddVertex(node1);
	//path->AddVertex(node2);

	//typedef itk::PathConstIterator< ImageType, PathType > PathIteratorType;
	//PathIteratorType path_iter(cell.threshold_label_mask, path);

	//itk::uint64_t path_length = 0;
	//itk::uint64_t num_blank_pixels = 0;
	//
	//while (!path_iter.IsAtEnd())
	//{
	//	//std::cout << "Path iterator position: " << path_iter.GetPathPosition() << std::endl;
	//	//std::cout << "Path iterator index: " << path_iter.GetIndex() << std::endl;
	//	if (path_iter.Get() == 0)
	//		++num_blank_pixels;
	//	++path_iter;
	//	++path_length;
	//}
	//
	//if (num_blank_pixels >= 20)
	//{
	//	//std::cerr << "Rejected a path" << std::endl;
	//	return std::numeric_limits< float >::max();	//Thresholding based on jumping large gaps
	//}
	//else
		return mag_trace_vector;
}

/* This function generates the minimum spanning tree (Prim's algorithm implementation, the output is a Tree structure */
Tree* MicrogliaRegionTracer::BuildMST1(Cell & cell, float** AdjGraph)
{	
	Tree* tree = new Tree();
	ImageType::IndexType root_index = cell.critical_points_queue.front();
	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], 1));

	//Make a copy of the AdjGraph called TreeAdjGraph which we will use to build the MST (Important that we keep the original AdjGraph to measure tortuosity)
	float** TreeAdjGraph = new float*[cell.critical_points_queue.size()];
	for (int k = 0; k < cell.critical_points_queue.size(); k++)
	{
		TreeAdjGraph[k] = new float[cell.critical_points_queue.size()];
		for (int l = 0; l < cell.critical_points_queue.size(); l++)
			TreeAdjGraph[k][l] = AdjGraph[k][l];
	}

	//Root node should have infinite weights to connect to
	for (int m = 0; m < cell.critical_points_queue.size(); m++)
		TreeAdjGraph[m][0] = std::numeric_limits< float >::max();

	//Prim's algorithm
	//for each node but the last, find the minimum weight unconnected node and connect it to the tree
	for (itk::uint64_t l = 0; l < cell.critical_points_queue.size() - 1; l++)
	{
		itk::uint64_t minimum_connected_node_id = 0;
		itk::uint64_t minimum_node_index_to_id = 0;
		float minimum_node_cost = std::numeric_limits< float >::max();
		Node* minimum_connected_node = NULL;
		std::vector< Node* > member_nodes = tree->GetMemberNodes();
		std::vector< Node* >::iterator member_nodes_iter;
        
		for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)	//For each node that is already part of the Tree in Prim's algorithm
		{
			Node* connected_node = *member_nodes_iter;

			for (itk::uint64_t k = 0; k < cell.critical_points_queue.size(); k++) //Search through all the unconnected nodes to this connected node
			{
				itk::uint64_t connected_node_id = connected_node->getID() - 1; //Node IDs are 1 greater than than vector indices, so we subtract 1 here
				
                //Weights for tortuosity and distance
				float alpha = 0.33;			//this affects the cost weight for the angle
				float beta = 1.0 - alpha;	//this affects the cost weight for the distance
								
				float node_to_candidate_length = TreeAdjGraph[connected_node_id][k];
				
                float cost;
				if (connected_node->GetParent() != NULL)	//is not the root node (has a parent)
				{
					itk::uint64_t parent_id = connected_node->GetParent()->getID() - 1;

					float parent_to_candidate_length = AdjGraph[parent_id][k];
					float parent_to_node_length = AdjGraph[parent_id][connected_node_id];
					
					//Law of cosines solved for the missing angle. PI then is subtracted from that, since it is deviation from a straight line.
					//Theta ranges between 0 and PI with 0 being the least deviation from a straight-line with the parent trace
					float theta = PI - acos((pow(parent_to_node_length, 2.0) + pow(node_to_candidate_length, 2.0) - pow(parent_to_candidate_length, 2.0)) / (2 * parent_to_node_length * node_to_candidate_length)); 
					cost = node_to_candidate_length * (alpha * theta / PI + beta);
				}
				else	//is the root node (has no parent)
				{	
					cost = beta * node_to_candidate_length;
				}
				
				if (cost < minimum_node_cost) //from l (connected point) to k (unconnected point) if the current cost is less than the minimum cost
				{
					minimum_connected_node = connected_node;
					minimum_connected_node_id = connected_node_id;			
					minimum_node_index_to_id = k;
					minimum_node_cost = cost;
				}
			}	
		}

		if (minimum_node_cost >= 200)
			break;	//Minimum distance way too far
        
        assert(minimum_node_cost < std::numeric_limits< float >::max() && minimum_connected_node != NULL);
		
		std::cout << "Found new edge from " << minimum_connected_node_id << " to " << minimum_node_index_to_id << " Location: " << cell.critical_points_queue[minimum_connected_node_id][0] << " " << cell.critical_points_queue[minimum_connected_node_id][1] << " " << cell.critical_points_queue[minimum_connected_node_id][2] << " " << cell.critical_points_queue[minimum_node_index_to_id][0] << " " << cell.critical_points_queue[minimum_node_index_to_id][1] << " "  << cell.critical_points_queue[minimum_node_index_to_id][2] << " cost: " << minimum_node_cost << std::endl;

		ImageType::IndexType point_index;
		point_index[0] = cell.critical_points_queue[minimum_node_index_to_id][0] + cell.getX() - cell.GetRequestedSize()[0]/2 - cell.GetShiftIndex()[0];
		point_index[1] = cell.critical_points_queue[minimum_node_index_to_id][1] + cell.getY() - cell.GetRequestedSize()[1]/2 - cell.GetShiftIndex()[1];
		point_index[2] = cell.critical_points_queue[minimum_node_index_to_id][2] + cell.getZ() - cell.GetRequestedSize()[2]/2 - cell.GetShiftIndex()[2];
		
		ImageType::IndexType point_local_index;
		point_local_index[0] = cell.critical_points_queue[minimum_node_index_to_id][0];
		point_local_index[1] = cell.critical_points_queue[minimum_node_index_to_id][1];
		point_local_index[2] = cell.critical_points_queue[minimum_node_index_to_id][2];

		
		for (int m = 0; m < cell.critical_points_queue.size(); m++)
			TreeAdjGraph[m][minimum_node_index_to_id] = std::numeric_limits< float >::max();						//Node already has parents, we dont want it to have more parents
		TreeAdjGraph[minimum_node_index_to_id][minimum_connected_node_id] = std::numeric_limits< float >::max();	//So the parent doesn't become the child of its child

		//Connect the unconnected point
		Node* new_connected_node = new Node(point_local_index[0], point_local_index[1], point_local_index[2], minimum_node_index_to_id + 1); //vector indices are 1 less than SWC IDs, so we add one here
		new_connected_node->SetParent(minimum_connected_node);
		minimum_connected_node->AddChild(new_connected_node);
		tree->AddNode(new_connected_node);
	}

	//Update next available ID
	cell.next_available_ID = cell.critical_points_queue.size() + 1;	//We used up IDs [1, critical_pts_queue.size()], so here we start at critical_points_queue.size + 1
    
    for (int k = 0; k < cell.critical_points_queue.size(); k++)
		delete[] TreeAdjGraph[k];
	delete[] TreeAdjGraph;
    
	return tree;
}

/* This function smoothes a Tree structure */
void MicrogliaRegionTracer::SmoothTree(Cell & cell, Tree* smoothed_tree)
{
	std::cout << "Smoothing the tree" << std::endl;
	SmoothSegments(cell, smoothed_tree, smoothed_tree->GetRoot());
}

/* The Tree segments are traversed here and SmoothPath is called on each segment */
void MicrogliaRegionTracer::SmoothSegments(Cell & cell, Tree* smoothed_tree, Node* start_node) //Smooth AND Prune
{
	std::vector< Node* > start_node_children = start_node->GetChildren();
	if (start_node_children.size() == 0)
		return;	//start_node has no children so it is a leaf node so there is no segment

	ImageType::IndexType start_node_index;
	start_node_index[0] = start_node->x;
	start_node_index[1] = start_node->y;
	start_node_index[2] = start_node->z;


	//Smooth all the paths that is connected to this start_node
	std::vector< Node* >::iterator start_node_children_iter;
	for (start_node_children_iter = start_node_children.begin(); start_node_children_iter != start_node_children.end(); ++start_node_children_iter)
	{
		float segment_length = 0.0;
		
		//Make a new path for each possible segment from this start_node
		PathType::Pointer path = PathType::New();
		path->Initialize();
		path->AddVertex(start_node_index);

		Node* start_node_child = *start_node_children_iter;	//Keep position of start_node_child so we can determine along which branch to smooth later
		Node* child_node = start_node_child;
		std::vector< Node* > child_node_children = start_node_child->GetChildren();
		ImageType::IndexType last_node_index = start_node_index;

		//std::cout << "Visiting path: " << start_node->getID() << " [" << start_node->x << ", " << start_node->y << ", " << start_node->z << "] ";
		
        //Keep going down the segment until we hit a branch point or the leaf node
		while (child_node_children.size() < 2 && child_node_children.size() != 0)
		{
			//std::cout << child_node->getID() << " [" << child_node->x << ", " << child_node->y << ", " << child_node->z << "] ";
			ImageType::IndexType child_node_index;
			child_node_index[0] = child_node->x;
			child_node_index[1] = child_node->y;
			child_node_index[2] = child_node->z;
			path->AddVertex(child_node_index);
			segment_length += CalculateEuclideanDistance(last_node_index, child_node_index);

			//Remove the child_node from the tree
			smoothed_tree->RemoveNode(child_node);

			//Updating variables to point to the next node along the segment
			last_node_index = child_node_index;
			child_node = child_node_children.front();
			child_node_children = child_node->GetChildren();
		}
				
		//If we got here, it means that we have at least 2 children or no children, so we should add ourselves to the end of the path, smooth it, and call SmoothSegments at the end node
		Node* end_node = child_node;
		ImageType::IndexType end_node_index;
		end_node_index[0] = end_node->x;
		end_node_index[1] = end_node->y;
		end_node_index[2] = end_node->z;
		path->AddVertex(end_node_index);
		//std::cout << end_node->getID() << " [" << end_node->x << ", " << end_node->y << ", " << end_node->z << "] " << std::endl;

		//Break the links at each node between the start node, start node child and the end node
		start_node->RemoveChild(start_node_child);
		Node* current_node = start_node_child;
		while (current_node != end_node)
		{
			assert(current_node->GetChildren().size() > 0);			//Can't break the child off of a leaf node, this means we are not checking for the right end node
			assert(current_node->GetChildren().size() < 2);			//We are breaking a link incorrectly if it has 2+ children
			
			Node* next_node = current_node->GetChildren().front();
			current_node->RemoveChild(next_node);

			current_node = next_node;								//Update current node to point to the next node
		}
		
		//If we have no children and we are short, then prune
		/*if (child_node_children.size() == 0 && segment_length < 5.0)
			continue;*/
		

		PathType::Pointer speed_path = SmoothPath(cell, smoothed_tree, start_node, end_node, path);	//Get the smoothed path
		
		ReplaceTreeSegmentWithPath(cell, smoothed_tree, speed_path, start_node, end_node);			//Replace the segment in the tree with the smoothed path

		SmoothSegments(cell, smoothed_tree, end_node);		//Call SmoothSegments on the next segment
	}
}

/* This function takes takes in start_node and end_node and a path in between and returns the smoothed path */
MicrogliaRegionTracer::PathType::Pointer MicrogliaRegionTracer::SmoothPath(Cell & cell, Tree* smoothed_tree, Node* start_node, Node* end_node, PathType::Pointer path )
{
	//Extract path from speed image (the cubed distance map)
	// Create interpolator
	typedef itk::SpeedFunctionToPathFilter< DistanceImageType, PathType >	SpeedPathFilterType;
	typedef SpeedPathFilterType::CostFunctionType::CoordRepType CoordRepType;
	typedef itk::LinearInterpolateImageFunction<DistanceImageType, CoordRepType> InterpolatorType;
	InterpolatorType::Pointer interp = InterpolatorType::New();

	// Create cost function
	SpeedPathFilterType::CostFunctionType::Pointer cost = SpeedPathFilterType::CostFunctionType::New();
	cost->SetInterpolator( interp );

	// Create optimizer
	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetNumberOfIterations( 1000000 );
	optimizer->SetMaximumStepLength( 0.01 );
	optimizer->SetMinimumStepLength( 0.01 );
	optimizer->SetRelaxationFactor( 0.00 );

	// Create path filter
	SpeedPathFilterType::Pointer pathFilter = SpeedPathFilterType::New();
	pathFilter->SetInput( cell.speed_image );
	pathFilter->SetCostFunction( cost );
	pathFilter->SetOptimizer( optimizer );
	pathFilter->SetTerminationValue( 10.0 );

	// Add path information
	SpeedPathFilterType::PathInfo info;
	SpeedPathFilterType::PointType start;
	start[0] = start_node->x; start[1] = start_node->y; start[2] = start_node->z;
	//std::cout << "Adding Smoothing startpoint index to: " << start << std::endl;
	

	const PathType::VertexListType *vertex_list = path->GetVertexList();

	for (unsigned int i = 1; i < vertex_list->Size() - 1; ++i)
	{
		PathType::ContinuousIndexType vertex = vertex_list->GetElement(i);
		//std::cout << "Adding smoothing waypoint index to: " << vertex << std::endl;
		info.AddWayPoint(vertex);
	}

	SpeedPathFilterType::PointType end;
	end[0] = end_node->x; end[1] = end_node->y; end[2] = end_node->z;
	//std::cout << "Adding smoothing endpoint index to: " << end << std::endl;
		
	info.SetStartPoint( end ); //the reversal of start and end is on purpose. There seems to be a bug in SpeedFunctionToPathFilter
	info.SetEndPoint( start );
	
	pathFilter->AddPathInfo( info );
	try
	{
		pathFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "pathFilter exception: " << err << std::endl;
		return NULL;
	}
	
	PathType::Pointer speed_path = pathFilter->GetOutput();

	return speed_path;

	//HACK: Seeks like the pathFilter potentially returns paths that will round to an index outside the image (ie. 99.5406....). Thus we take original vertex list and copy it over while clamping the values
	ImageType::RegionType image_region = cell.image->GetLargestPossibleRegion();
	ImageType::SizeType image_size = image_region.GetSize();
	PathType::Pointer clamped_speed_path = PathType::New();

	const PathType::VertexListType * smoothed_vertices_list = speed_path->GetVertexList();
		
	for (unsigned int i = 0; i < smoothed_vertices_list->Size(); ++i)
	{
		PathType::ContinuousIndexType vertex = smoothed_vertices_list->GetElement(i);
		
		for (int k = 0; k < 3; ++k)
		{
			vertex[k] = std::min< double >(image_size[k], vertex[k]);
		}
		clamped_speed_path->AddVertex(vertex);
	}

	return clamped_speed_path;
}

/* This function modifies the tree passed to it and replaces the path between start_node and end_node with the speed_path */
void MicrogliaRegionTracer::ReplaceTreeSegmentWithPath(Cell & cell, Tree* smoothed_tree, PathType::Pointer speed_path, Node* start_node, Node* end_node)
{
	////Break the links between start_node and end_node and delete all the nodes in between
	//Node * delete_node = end_node->GetParent();
	//while (delete_node != start_node)
	//{
	//	smoothed_tree->RemoveNode(delete_node);
	//	Node * delete_node_parent = delete_node->GetParent();
	//	delete delete_node;
	//	delete_node = delete_node_parent;
	//	
	//}
	
	//Convert path output back to tree format
	if ( speed_path->GetVertexList()->Size() <= 2 )
	{
		std::cerr << "WARNING: Path contains not enough points!" << std::endl;
		start_node->AddChild(end_node);
		end_node->SetParent(start_node);
		return;
	}
	else
	{
		std::cerr << "PATH CONTAINS " << speed_path->GetVertexList()->Size() << " SMOOTHED POINTS!" << std::endl;
	}

	const PathType::VertexListType *smoothed_vertex_list = speed_path->GetVertexList();

	Node* current_node = start_node;
	for (unsigned int i = 0; i < smoothed_vertex_list->Size(); ++i)
	{
		PathType::ContinuousIndexType path_index = smoothed_vertex_list->GetElement(i);

		Node* new_node = new Node(path_index[0], path_index[1], path_index[2], cell.next_available_ID );

		//std::cout << "Adding smoothed path index: " << path_index << std::endl;

		//Update next available ID
		++(cell.next_available_ID);

		//Connect this new node to the current_node
		current_node->AddChild(new_node);
		new_node->SetParent(current_node);
		smoothed_tree->AddNode(new_node);

		current_node = new_node;
	}
	//std::cout << "Connecting last index: " << end_node->x << " " << end_node->y << " " << end_node->z << std::endl;
	current_node->AddChild(end_node);
	end_node->SetParent(current_node);
}

/* Calculate the Euclidean Distance */
float CalculateEuclideanDistance(Cell::ImageType::IndexType node1, Cell::ImageType::IndexType node2)
{
    Cell::ImageType::IndexType trace_vector;
	trace_vector[0] = node2[0] - node1[0];
	trace_vector[1] = node2[1] - node1[1];
	trace_vector[2] = node2[2] - node1[2];

	return sqrt(pow(trace_vector[0], 2.0f) + pow(trace_vector[1], 2.0f) + pow(trace_vector[2], 2.0f));
}

void MicrogliaRegionTracer::PruneTree(Cell & cell, Tree* smoothed_tree)
{
	std::cout << "Pruning the smoothed tree" << std::endl;
	CheckSegment(cell, smoothed_tree->GetRoot(), smoothed_tree);	
}

void MicrogliaRegionTracer::CheckSegment(Cell & cell, Node* start_node, Tree * smoothed_tree)
{
	typedef std::vector< Node *> NodeVectorType;
	NodeVectorType start_node_children = start_node->GetChildren();

	//For each of our children, make a new path and populate the path
	for (NodeVectorType::iterator start_node_children_iter = start_node_children.begin(); start_node_children_iter != start_node_children.end(); ++start_node_children_iter)
	{
		PathType::Pointer path = PathType::New();

		PathType::ContinuousIndexType start_node_index;
		start_node_index[0] = start_node->GetX();
		start_node_index[1] = start_node->GetY();
		start_node_index[2] = start_node->GetZ();

		//Add the start_node into the path
		path->AddVertex(start_node_index);
		
		Node * start_node_child = *start_node_children_iter;
		
		Node * node = start_node_child;
		NodeVectorType node_children = node->GetChildren();
		
		//std::cout << "Added path with start node ID: " << start_node->getID();
		//Add nodes until you reach a node that has 0 or more than 1 child
		while (node_children.size() == 1)
		{
			node = node_children.front();

			PathType::ContinuousIndexType node_index;
			node_index[0] = node->GetX();
			node_index[1] = node->GetY();
			node_index[2] = node->GetZ();

			path->AddVertex(node_index);

			node_children = node->GetChildren();
			//std::cout << " node ID: " << node->getID();
		}
		

		//Add the end_node to the path
		Node* end_node = node;
		PathType::ContinuousIndexType end_node_index;
		end_node_index[0] = end_node->GetX();
		end_node_index[1] = end_node->GetY();
		end_node_index[2] = end_node->GetZ();
		path->AddVertex(end_node_index);
		//std::cout << " and end node ID: " << end_node->getID() << std::endl; 

		
		//std::cout << "Added " << path->GetVertexList()->Size() << " nodes to the prune path" << std::endl;

		CheckSegment(cell, end_node, smoothed_tree); //Recursively go check the next segment

		NodeVectorType end_node_children = end_node->GetChildren();

		if (end_node_children.size() != 0)
			continue;				//BASE CASE: We have descendant segments which remain after pruning, so don't prune us

		//Here we are a candidate for pruning because we have no descendant segments
		bool should_prune = ShouldPruneSegment(cell, path);
		//std::cout << "Should we prune this segment? " << should_prune << std::endl;

		if (should_prune)
		{
			Node * node = end_node;
			while (node != start_node)
			{
				//Keep reference to the node to delete
				Node * node_to_delete = node;
				
				//Update the parent node so it doesn't hold onto the node we are deleting
				Node * parent_node = node->GetParent();
				parent_node->RemoveChild(node_to_delete);
				
				//Remove the node to delete from the tree members
				smoothed_tree->RemoveNode(node_to_delete);

				//Update the reference of node
				node = parent_node;

				//Actually delete that node
				delete node_to_delete;
			}
		}
	}

	//BASE CASE: Completely done with all pruning attempts rooted at this subtree
}

bool MicrogliaRegionTracer::ShouldPruneSegment(Cell & cell, PathType::Pointer path)
{
	typedef itk::PathConstIterator< ImageType, PathType > PathIteratorType;
	PathIteratorType path_iter(cell.threshold_label_mask, path);

	itk::uint64_t path_length = 0;
	double num_blank_pixels = 0;

	//std::cout << "VertexList size is: " << path->GetVertexList()->Size() << std::endl;
	
	//std::cout << "End time-point: " << path->EndOfInput() << std::endl;
	//std::cout << "End Index: " << path->EvaluateToIndex(path->EndOfInput()) << std::endl;   

	while (!path_iter.IsAtEnd())
	{
		//std::cout << "Visiting path location: " << path_iter.GetIndex() << std::endl;
		unsigned char pixel_label = path_iter.Get();
		
		if (pixel_label == 1)
			num_blank_pixels += 0.0;
		else if (pixel_label == 0)
			num_blank_pixels += 1.0;
		
		++path_iter;
		++path_length;
	}
	
	if (path_length < 5 )
		return true;

	if (num_blank_pixels >= 2)
		return true;

		
	return false;

}