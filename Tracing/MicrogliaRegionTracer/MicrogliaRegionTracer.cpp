#include "MicrogliaRegionTracer.h"

#define MASK 0
#define PRINT_ALL_IMAGES 0
#define PI (4.0*atan(1.0))

MicrogliaRegionTracer::MicrogliaRegionTracer(const std::string & joint_transforms_filename, const std::string & img_path, const std::string & anchor_filename, const std::string & soma_filename)
{
	this->roi_grabber = new ROIGrabber(joint_transforms_filename, img_path, anchor_filename);
	this->soma_filename = soma_filename;

#ifdef _OPENMP
	std::cerr << "OpenMP is on, setting itkGlobalDefaultNumberOfThreads to 1" << std::endl;
	itk::MultiThreader::SetGlobalDefaultNumberOfThreads( 1 );	//Acquiring threads through OpenMP, so no need for ITK threads
#endif
	aspect_ratio = 3.0; //xy to z ratio in physical space here
}

MicrogliaRegionTracer::~MicrogliaRegionTracer()
{
    delete this->roi_grabber;
}

/*	This function takes in the name of the seed points (where tracing starts) and then grabs a 200x200x100 region with mosaic_roi.
	A new Cell object is created to hold the image and seed point and relevant information. */
void MicrogliaRegionTracer::LoadCellPoints(const std::string & seedpoints_filename)
{
	std::ifstream seed_point_file;
	seed_point_file.open(seedpoints_filename.c_str());

	while (!seed_point_file.eof())
	{	
		itk::uint64_t cellX, cellY, cellZ;
		seed_point_file >> cellX >> cellY >> cellZ;
		//std::cout << "Reading in cell: (" << cellX << ", " << cellY << ", " << cellZ << ")" << std::endl;
		
		Cell* cell = new Cell(cellX, cellY, cellZ);

		cells.push_back(cell);
	}
}

/* This is the main loop where tracing takes place */
void MicrogliaRegionTracer::Trace()
{
	//Create groups of n cells at a time (where n is the number of logical processors)
    //This is needed because OpenMP doesn't play well with ITK
    int num_threads = 1;    //Default 1 thread in a group
#ifdef _OPENMP
    num_threads = omp_get_num_procs();
#endif
    
    int num_groups = std::ceil(cells.size() / (double) num_threads);
	std::cerr << "Number of cells in each group (except the last group): " << num_threads << std::endl;
    std::cerr << "Number of groups of cells to process: " << num_groups << std::endl;

    for (int group_num = 0; group_num < num_groups; group_num++)
    {
        int num_local_threads = (cells.size() / num_threads) ? num_threads : cells.size() % num_threads; //sets num_local_threads to the number of threads if we are not on the last group, otherwise set it to the number of remaining threads
		std::cerr << "number of local threads: " << num_local_threads << std::endl;
		for (int local_thread_num = 0; local_thread_num < num_local_threads; local_thread_num++)
        {
            int global_thread_num = group_num * num_threads + local_thread_num;
            
            Cell* cell = cells[global_thread_num];
            std::cerr << "Reading cell #: " << global_thread_num << std::endl;
                      
            //Setup the size of the initial ROI
            ImageType::SizeType roi_size;
            roi_size[0] = 200;
            roi_size[1] = 200;
            roi_size[2] = 100;
            
            cell->setRequestedSize(roi_size);
            
            //Grab the initial cellimage
            //std::cout << "Grabbing ROI for cell" << std::endl;
            ImageType::IndexType shift_index;
            ImageType::Pointer temp_cell_image = roi_grabber->GetROI(cell, roi_size, shift_index);
            
            //Set the origin of the image to be {0, 0, 0}
            ImageType::PointType origin = temp_cell_image->GetOrigin();
            cell->SetOrigin(origin);
            origin[0] = 0;
            origin[1] = 0;
            origin[2] = 0;
            temp_cell_image->SetOrigin(origin);
            
            //Duplicate the image so that the requests don't propagate back up to roi_grabber
            typedef itk::ImageDuplicator< ImageType > DuplicatorType;
            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(temp_cell_image);
            ftk::TimeStampOverflowSafeUpdate( duplicator.GetPointer() );
            cell->image = duplicator->GetOutput();
            cell->image->DisconnectPipeline();
            
            //Shift_index is just to demonstrate how much the center point of the local point of the image if it got cropped by the left, top, closest point of the image
            cell->setShiftIndex(shift_index);
            //std::cout << cell->getShiftIndex()[0] << " " << cell->getShiftIndex()[1] << " " << cell->getShiftIndex()[2] << std::endl;
            
            //Make the file name of the raw cell image
            std::stringstream cell_filename_stream;
            cell_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".TIF";	//X_Y_Z.TIF
            
            //Write the cell image
            cell->WriteImage(cell_filename_stream.str(), cell->image);
            
            roi_size = cell->image->GetLargestPossibleRegion().GetSize();	//The size of the returned image may not be the size of the image that entered because clipping at edges
            cell->SetSize(roi_size);
        }
        
        //Trace a group of cells
        #pragma omp parallel for
        for (int local_thread_num = 0; local_thread_num < num_local_threads; local_thread_num++)
        {
            int global_thread_num = group_num * num_threads + local_thread_num;
            Cell* cell = cells[global_thread_num];
            
            //Get the mask
        #if MASK	//development only
            cell->GetMask(this->soma_filename);
        #endif
            std::cout << "Calculating candidate pixels for a new cell" << std::endl;
            CalculateCandidatePixels(cell);
            
            std::cout << "Detected " << cell->critical_points_queue.size() << " critical points" << std::endl;
            
            std::cout << "Tree Building" << std::endl;
            BuildTree(cell);
            
            delete cell;
        }
    }
}

/* This function determines the candidate pixels (pixels which we connect to form the tree) */
void MicrogliaRegionTracer::CalculateCandidatePixels(Cell* cell)
{
	CreateIsometricImage(cell);
	
	std::cout << "VesselnessDetection" << std::endl;
	VesselnessDetection(cell);

	std::cout << "Ridge Detection" << std::endl;
	RidgeDetection(cell);
}

/* This function resamples the image and creates an "isometric" image according to the ratio set in the constructor */
void MicrogliaRegionTracer::CreateIsometricImage(Cell* cell)
{
	ImageType::SizeType inputSize = cell->image->GetLargestPossibleRegion().GetSize();
	ImageType::SizeType outputSize = inputSize;
	outputSize[2] = outputSize[2] * 1/aspect_ratio;

	ImageType::SpacingType outputSpacing;
	outputSpacing[0] = 1.0;
	outputSpacing[1] = 1.0;
	outputSpacing[2] = aspect_ratio;

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(cell->image);
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());
	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "resample_filter exception: " << err << std::endl;
		return;
	}

	cell->isometric_image = resample_filter->GetOutput();
	cell->isometric_image->DisconnectPipeline();
#if PRINT_ALL_IMAGES
	cell->WriteImage("isometric_image.mhd", cell->isometric_image);
#endif

	
	ImageType::SpacingType spacing;
	spacing.Fill(1.0);
	cell->isometric_image->SetSpacing(spacing);
}

/* This function does Ridge detection (Laplacian of Gaussian implementation) */
void MicrogliaRegionTracer::RidgeDetection( Cell* cell )
{
	//Add the seed to the critical points
	itk::Index<3> seed_index = {{	cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0], 
									cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1], 
									cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2] }};

	cell->critical_points_queue.push_front(seed_index);
	
	//Calculate the LoG on multiple scales and store into an image
	LoG *log_obj = new LoG();
	std::cout << "Calculating Multiscale LoG" << std::endl;
	LoGImageType::Pointer resampled_multiscale_LoG_image = log_obj->RunMultiScaleLoG(cell);
    delete log_obj;
    
#if PRINT_ALL_IMAGES
	cell->WriteImage("multiscaled_LoG_image.mhd", resampled_multiscale_LoG_image);
#endif

	ImageType::SizeType outputSize = cell->image->GetLargestPossibleRegion().GetSize();

	ImageType::SpacingType outputSpacing;
	outputSpacing[0] = 1.0;
	outputSpacing[1] = 1.0;
	outputSpacing[2] = 1/aspect_ratio;
	//outputSpacing[2] = 1.0;

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< LoGImageType, LoGImageType > ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(resampled_multiscale_LoG_image);
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());

	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "resample_filter exception: " << err << std::endl;
		return;
	}

	cell->multiscale_LoG_image = resample_filter->GetOutput();
    
    
	ImageType::SpacingType spacing;
	spacing.Fill(1.0);
	cell->multiscale_LoG_image->SetSpacing(spacing);
	
	//Make a new image to store the ridge pixels
	ImageType::SizeType size = cell->multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	cell->ridge_image = LoGImageType::New();
	ImageType::IndexType start;
	start.Fill(0);
	ImageType::RegionType region(start, size);
	cell->ridge_image->SetRegions(region);
	cell->ridge_image->Allocate();
	cell->ridge_image->FillBuffer(0);

    //Non-maximal suppression
	//Make a iterator for the image and a neighborhood around the current point we are visiting
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator< LoGImageType > ridge_img_iter(cell->ridge_image, cell->ridge_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > LoG_neighbor_iter(rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());

	while(!LoG_neighbor_iter.IsAtEnd()) 
	{
        LoGImageType::PixelType a = LoG_neighbor_iter.GetPixel(4);
		
		if (a > 0.01) //Ignore pixels less than 1% response
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
	ImageType::SizeType LoG_image_size = cell->multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	cell->critical_point_image = ImageType::New();
	ImageType::IndexType LoG_image_start;
	LoG_image_start.Fill(0);
	ImageType::RegionType LoG_image_region(LoG_image_start, LoG_image_size);
	cell->critical_point_image->SetRegions(LoG_image_region);
	cell->critical_point_image->Allocate();
	cell->critical_point_image->FillBuffer(0);

	//Make a iterator for the image and a neighborhood around the current point we are visiting
	itk::Size<3> LoG_image_rad = {{1,1,1}};
	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell->critical_point_image, cell->critical_point_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > ridge_neighbor_iter(rad, cell->ridge_image, cell->ridge_image->GetLargestPossibleRegion());

	while(!ridge_neighbor_iter.IsAtEnd()) 
	{
		unsigned int neighborhood_size = (LoG_image_rad[0] * 2 + 1) * (LoG_image_rad[1] * 2 + 1) * (LoG_image_rad[2] * 2 + 1);
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
				cell->critical_points_queue.push_back(critical_point_img_iter.GetIndex());
			}
		}

		++critical_point_img_iter;
		++ridge_neighbor_iter;
	}
    
        

	std::ostringstream criticalpointsFileNameStream;
	criticalpointsFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_critical.mhd";
#if PRINT_ALL_IMAGES
	cell->WriteImage(criticalpointsFileNameStream.str(), cell->critical_point_image);
#endif
}

/* This function calculates the Vesselness score */
void MicrogliaRegionTracer::VesselnessDetection(Cell* cell)
{
	typedef itk::Hessian3DToVesselnessMeasureImageFilter< float > VesselnessFilterType;
	VesselnessFilterType::Pointer vesselness_filter = VesselnessFilterType::New();
	vesselness_filter->SetAlpha1(0.5);
	vesselness_filter->SetAlpha2(0.5);

	typedef itk::SymmetricSecondRankTensor< double, 3 >	HessianTensorType;
	typedef itk::Image< HessianTensorType, 3 >			HessianImageType;

	typedef itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, VesselnessImageType > MultiscaleHessianFilterType;

	MultiscaleHessianFilterType::Pointer multiscale_hessian_filter = MultiscaleHessianFilterType::New();
	//multiscale_hessian_filter->SetInput(cell->image);
	multiscale_hessian_filter->SetInput(cell->isometric_image);
	multiscale_hessian_filter->SetSigmaStepMethodToEquispaced();
	multiscale_hessian_filter->SetSigmaMinimum(1.0);
	multiscale_hessian_filter->SetSigmaMaximum(10.0);
	multiscale_hessian_filter->SetNumberOfSigmaSteps(10);
	multiscale_hessian_filter->SetNonNegativeHessianBasedMeasure(true);
	multiscale_hessian_filter->SetHessianToMeasureFilter(vesselness_filter);

	try
	{
		ftk::TimeStampOverflowSafeUpdate( multiscale_hessian_filter.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "multiscale_hessian_filter exception: " << err << std::endl;
	}

	//cell->WriteImage("Resampled_vesselness.mhd", multiscale_hessian_filter->GetOutput());

	//Unsample the image
	ImageType::SizeType outputSize = cell->image->GetLargestPossibleRegion().GetSize();

	ImageType::SpacingType outputSpacing;
	outputSpacing[0] = 1.0;
	outputSpacing[1] = 1.0;
	outputSpacing[2] = 1/aspect_ratio;
	//outputSpacing[2] = 1.0;

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< VesselnessImageType, VesselnessImageType > ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(multiscale_hessian_filter->GetOutput());
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());
	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "resample_filter exception: " << err << std::endl;
		return;
	}

	cell->vesselness_image = resample_filter->GetOutput();
	cell->vesselness_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...
	ImageType::SpacingType spacing;
	spacing.Fill(1.0);
	cell->vesselness_image->SetSpacing(spacing);
	
	std::ostringstream vesselnessFileNameStream;

	vesselnessFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_vesselness.mhd";
#if PRINT_ALL_IMAGES
		cell->WriteImage(vesselnessFileNameStream.str(), cell->vesselness_image);
#endif
}

/* After the candidates pixels are calculated, this function connects all the candidate pixels into a minimum spanning tree based on a cost function */
void MicrogliaRegionTracer::BuildTree(Cell* cell)
{
	std::cerr << "Building Adjacency Graph" << std::endl;
	double** AdjGraph = BuildAdjacencyGraph(cell);

	//for (int k = 0; k < cell->critical_points_queue.size(); k++)
	//	std::cout << AdjGraph[0][k] << std::endl;

	std::cerr << "Building Minimum Spanning Tree" << std::endl;
	Tree* tree = BuildMST1(cell, AdjGraph);

	std::ostringstream swc_filename_stream, swc_filename_stream_local;
	swc_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree.swc";
	swc_filename_stream_local << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree_local.swc";


	WriteTreeToSWCFile(tree, cell, swc_filename_stream.str(), swc_filename_stream_local.str());

	Tree* smoothed_tree = new Tree(*tree);
	SmoothTree(cell, smoothed_tree);

	std::ostringstream swc_filename_stream_smoothed, swc_filename_stream_local_smoothed;
	swc_filename_stream_smoothed << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree_smoothed.swc";
	swc_filename_stream_local_smoothed << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree_local_smoothed.swc";


	WriteTreeToSWCFile(smoothed_tree, cell, swc_filename_stream_smoothed.str(), swc_filename_stream_local_smoothed.str());

	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		delete[] AdjGraph[k];
	delete[] AdjGraph;

	delete tree;
	delete smoothed_tree;
}

/* This function generates a matrix of all possible pairs of candidate pixels */
double** MicrogliaRegionTracer::BuildAdjacencyGraph(Cell* cell)
{
	double** AdjGraph = new double*[cell->critical_points_queue.size()];
	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		AdjGraph[k] = new double[cell->critical_points_queue.size()];

	#pragma omp parallel for
	for (itk::int64_t k = 0; k < (itk::int64_t)cell->critical_points_queue.size(); k++)
	{
		//std::cout << "Calculating distance for node " << k << std::endl;
		for (itk::uint64_t l = 0; l < cell->critical_points_queue.size(); l++)
		{
			AdjGraph[k][l] = CalculateDistance(k, l, cell);
		}
	}

	std::cout << "Done calculating distance graph" << std::endl;
	return AdjGraph;
}

/* This is the distance part of the cost function */
double MicrogliaRegionTracer::CalculateDistance(itk::uint64_t node_from, itk::uint64_t node_to, Cell* cell)
{
	ImageType::IndexType node1 = cell->critical_points_queue[node_from];
	ImageType::IndexType node2 = cell->critical_points_queue[node_to];

	ImageType::IndexType trace_vector;
	trace_vector[0] = node2[0] - node1[0];
	trace_vector[1] = node2[1] - node1[1];
	trace_vector[2] = node2[2] - node1[2];

	double mag_trace_vector = sqrt(pow(trace_vector[0], 2.0) + pow(trace_vector[1], 2.0) + pow(trace_vector[2], 2.0));
	
	if ((node_from == node_to) || (mag_trace_vector < 0.0001))	//Very likely that these two points are equal so we shouldn't try to connect them
		return std::numeric_limits< double >::max();

#if MASK
	//If we are connecting from the root node
	if (node_from == 0)
	{	
		if (cell->soma_label_image->GetPixel(node1) == cell->soma_label_image->GetPixel(node2)) //If it is inside the soma our centroid belongs to, we give it 0 weight
			return 0;
		else if (cell->soma_label_image->GetPixel(node2) != 0)											//We are trying to connect to a pixel in another soma, so we give it infinite weight
			return std::numeric_limits< double >::max();	
	}
#endif

	PathType::Pointer path = PathType::New();
	path->Initialize();

	//std::cout << "Start Index: " << node1 << " " << "End Index: " << node2 << std::endl;

	path->AddVertex(node1);
	path->AddVertex(node2);

	typedef itk::PathConstIterator< VesselnessImageType, PathType > PathIteratorType;
	PathIteratorType path_iter(cell->vesselness_image, path);

	double sum_of_vesselness_values = 0;
	itk::uint64_t path_length = 0;
	itk::uint64_t num_blank_pixels = 0;
	
	while (!path_iter.IsAtEnd())
	{
		//std::cout << "Path iterator position: " << path_iter.GetPathPosition() << std::endl;
		//std::cout << "Path iterator index: " << path_iter.GetIndex() << std::endl;
		if (path_iter.Get() < 1.0)
			++num_blank_pixels;
		++path_iter;
		++path_length;
	}
	
	if (num_blank_pixels >= 3)
		return std::numeric_limits< double >::max();	//Thresholding based on jumping large gaps
	else
		return mag_trace_vector;
}

/* This function generates the minimum spanning tree (Prim's algorithm implementation, the output is a Tree structure */
Tree* MicrogliaRegionTracer::BuildMST1(Cell* cell, double** AdjGraph)
{	
	Tree* tree = new Tree();
	ImageType::IndexType root_index = cell->critical_points_queue.front();
	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], 1));

	//Make a copy of the AdjGraph called TreeAdjGraph which we will use to build the MST (Important that we keep the original AdjGraph to measure tortuosity)
	double** TreeAdjGraph = new double*[cell->critical_points_queue.size()];
	for (int k = 0; k < cell->critical_points_queue.size(); k++)
	{
		TreeAdjGraph[k] = new double[cell->critical_points_queue.size()];
		for (int l = 0; l < cell->critical_points_queue.size(); l++)
			TreeAdjGraph[k][l] = AdjGraph[k][l];
	}	

	//Root node should have infinite weights to connect to
	for (int m = 0; m < cell->critical_points_queue.size(); m++)
		TreeAdjGraph[m][0] = std::numeric_limits<double>::max();

	//Prim's algorithm
	//for each node but the last, find the minimum weight unconnected node and connect it to the tree
	for (itk::uint64_t l = 0; l < cell->critical_points_queue.size() - 1; l++)
	{
		itk::uint64_t minimum_connected_node_id = 0;
		itk::uint64_t minimum_node_index_to_id = 0;
		double minimum_node_cost = std::numeric_limits<double>::max();
		Node* minimum_connected_node = NULL;
		std::vector<Node*> member_nodes = tree->GetMemberNodes();
		std::vector<Node*>::iterator member_nodes_iter;
        
		for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)	//For each node that is already part of the Tree in Prim's algorithm
		{
			Node* connected_node = *member_nodes_iter;

			for (itk::uint64_t k = 0; k < cell->critical_points_queue.size(); k++) //Search through all the unconnected nodes to this connected node
			{
				itk::uint64_t connected_node_id = connected_node->getID() - 1; //Node IDs are 1 greater than than vector indices, so we subtract 1 here
				
                //Weights for tortuosity and distance
				double alpha = 0.5;			//this affects the cost weight for the angle
				double beta = 1.0 - alpha;	//this affects the cost weight for the distance
								
				double node_to_candidate_length = TreeAdjGraph[connected_node_id][k];
				
                double cost;
				if (connected_node->GetParent() != NULL)	//is not the root node (has a parent)
				{
					itk::uint64_t parent_id = connected_node->GetParent()->getID() - 1;

					double parent_to_candidate_length = AdjGraph[parent_id][k];
					double parent_to_node_length = AdjGraph[parent_id][connected_node_id];
					
					//Law of cosines solved for the missing angle. PI then is subtracted from that, since it is deviation from a straight line.
					//Theta ranges between 0 and PI with 0 being the least deviation from a straight-line with the parent trace
					double theta = PI - acos((pow(parent_to_node_length, 2.0) + pow(node_to_candidate_length, 2.0) - pow(parent_to_candidate_length, 2.0)) / (2 * parent_to_node_length * node_to_candidate_length)); 
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
		
		std::cout << "Found new edge from " << minimum_connected_node_id << " to " << minimum_node_index_to_id << " Location: " << cell->critical_points_queue[minimum_connected_node_id][0] << " " << cell->critical_points_queue[minimum_connected_node_id][1] << " " << cell->critical_points_queue[minimum_connected_node_id][2] << " " << cell->critical_points_queue[minimum_node_index_to_id][0] << " " << cell->critical_points_queue[minimum_node_index_to_id][1] << " "  << cell->critical_points_queue[minimum_node_index_to_id][2] << " cost: " << minimum_node_cost << std::endl;

		ImageType::IndexType point_index;
		point_index[0] = cell->critical_points_queue[minimum_node_index_to_id][0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
		point_index[1] = cell->critical_points_queue[minimum_node_index_to_id][1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
		point_index[2] = cell->critical_points_queue[minimum_node_index_to_id][2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];
		
		ImageType::IndexType point_local_index;
		point_local_index[0] = cell->critical_points_queue[minimum_node_index_to_id][0];
		point_local_index[1] = cell->critical_points_queue[minimum_node_index_to_id][1];
		point_local_index[2] = cell->critical_points_queue[minimum_node_index_to_id][2];

		
		for (int m = 0; m < cell->critical_points_queue.size(); m++)
			TreeAdjGraph[m][minimum_node_index_to_id] = std::numeric_limits<double>::max();						 //Node already has parents, we dont want it to have more parents
		TreeAdjGraph[minimum_node_index_to_id][minimum_connected_node_id] = std::numeric_limits<double>::max(); //So the parent doesn't become the child of its child

		//Connect the unconnected point
		Node* new_connected_node = new Node(point_local_index[0], point_local_index[1], point_local_index[2], minimum_node_index_to_id + 1); //vector indices are 1 less than SWC IDs, so we add one here
		new_connected_node->SetParent(minimum_connected_node);
		minimum_connected_node->AddChild(new_connected_node);
		tree->AddNode(new_connected_node, minimum_connected_node);
	}

	//Update next available ID
	cell->next_available_ID = cell->critical_points_queue.size() + 1;	//We used up IDs [1, critical_pts_queue.size()], so here we start at critical_points_queue.size + 1
    
    for (int k = 0; k < cell->critical_points_queue.size(); k++)
		delete[] TreeAdjGraph[k];
	delete[] TreeAdjGraph;
    
	return tree;
}

/*	This function will take a tree structure and write out the "local" and "global" SWC file corresponding to that tree structure.
*	Local here means an SWC file that will fit into the region of interest we are interested in and global means that it will fit into the montage */ 
void MicrogliaRegionTracer::WriteTreeToSWCFile(Tree* tree, Cell* cell, std::string filename, std::string filename_local)
{
	std::cout << "Entering WriteTreeToSWCFile" << std::endl;
	std::ofstream traceFile, traceFile_local;
	
	std::cout << "Opening " << filename << std::endl;
	traceFile.open(filename.c_str());
	std::cout << "Opening " << filename_local << std::endl;
	traceFile_local.open(filename_local.c_str());
	
	Node* root = tree->GetRoot();

	itk::uint64_t tree_depth = 0; //root node is defined as tree depth 0
	WriteLinkToParent(root, tree_depth, cell, traceFile, traceFile_local);	//Recursive function that does the actual tree traversal and writing the SWC lines

	traceFile.close();
	traceFile_local.close();	
}

//This function does a depth-first traversal of the tree, maybe it makes sense to do a breadth-first traversal for some uses?
//CAREFUL: This function may overflow the stack due to recursion pushing parameters onto the stack and will crash if you have a tree deep enough. Rewrite iteratively or increase stack size if this is the case...
void MicrogliaRegionTracer::WriteLinkToParent(Node* node, itk::uint64_t tree_depth, Cell* cell, std::ofstream &traceFile, std::ofstream &traceFile_local)
{
	//Calculate some node indices
	ImageType::PointType node_index, node_index_local;
	node_index_local[0] = node->x;
	node_index_local[1] = node->y;
	node_index_local[2] = node->z;
	node_index[0] = node_index_local[0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
	node_index[1] = node_index_local[1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
	node_index[2] = node_index_local[2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];

	itk::int64_t parent_node_id;
	if (node->GetParent() == NULL) 
		parent_node_id = -1; //This node has no parent, so this is the root node, so parent ID is -1
	else
		parent_node_id = node->GetParent()->getID();

	std::vector< Node* > children = node->GetChildren();
	
	if (tree_depth == 1 && children.size() == 0)
		return;														//BASE CASE: Don't write out a trace if we are at depth one and have no children because we are a trace to the edge of the soma

	//Write out the SWC lines
	traceFile		<< node->getID() << " 3 " << node_index[0]			<< " " << node_index[1]			<< " "  << node_index[2]		<< " 1 " << parent_node_id << std::endl;
	traceFile_local << node->getID() << " 3 " << node_index_local[0]	<< " " << node_index_local[1]	<< " "  << node_index_local[2]	<< " 1 " << parent_node_id << std::endl;
	
	if (children.size() == 0) 
		return;														//BASE CASE: No children, so don't visit them

	std::vector< Node* >::iterator children_iter;
	for (children_iter = children.begin(); children_iter != children.end(); ++children_iter)
	{
		Node* child_node = *children_iter;
		WriteLinkToParent(child_node, tree_depth+1, cell, traceFile, traceFile_local);
	}

	return;															//BASE CASE: Finished visiting the entire subtree that is rooted at this node
}

/* This function smoothes a Tree structure */
void MicrogliaRegionTracer::SmoothTree(Cell* cell, Tree* smoothed_tree)
{
	CreateSpeedImage(cell);

	SmoothSegments(cell, smoothed_tree, smoothed_tree->GetRoot());
	//Tree* smoothed_tree = SmoothSegments2(cell, tree);
}

/* The Tree segments are traversed here and SmoothPath is called on each segment */
void MicrogliaRegionTracer::SmoothSegments(Cell* cell, Tree* smoothed_tree, Node* start_node) //Smooth AND Prune
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
		double segment_length = 0.0;
		
		//Make a new path for each possible segment from this start_node
		PathType::Pointer path = PathType::New();
		path->Initialize();
		path->AddVertex(start_node_index);

		Node* start_node_child = *start_node_children_iter;	//Keep position of start_node_child so we can determine along which branch to smooth later
		Node* child_node = start_node_child;
		std::vector< Node* > child_node_children = start_node_child->GetChildren();
		ImageType::IndexType last_node_index = start_node_index;

		std::cout << "Visiting path: " << start_node->getID() << " [" << start_node->x << ", " << start_node->y << ", " << start_node->z << "] ";
		
        //Keep going down the segment until we hit a branch point or the leaf node
		while (child_node_children.size() < 2 && child_node_children.size() != 0)
		{
			std::cout << child_node->getID() << " [" << child_node->x << ", " << child_node->y << ", " << child_node->z << "] ";
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
		std::cout << end_node->getID() << " [" << end_node->x << ", " << end_node->y << ", " << end_node->z << "] " << std::endl;

		//Break the links at each node between the start node, start node child and the end node
		start_node->RemoveChild(start_node_child);
		Node* current_node = start_node_child;
		while (current_node != end_node)
		{
			assert(current_node->GetChildren().size() > 0);			//Can't break the child off of a leaf node, this means we are not checking for the right end node
			assert(current_node->GetChildren().size() < 2);			//We are breaking a link incorrectly if it has 2+ children
			
			Node* next_node = current_node->GetChildren().front();
			current_node->RemoveChild(next_node);
			current_node = next_node;							//Update current node to point to the next node
		}
		
		//If we have no children and we are short, then prune
		if (child_node_children.size() == 0 && segment_length < 5.0)
			continue;
		

		PathType::Pointer speed_path = SmoothPath(cell, smoothed_tree, start_node, end_node, path);	//Get the smoothed path
		
		ReplaceTreeSegmentWithPath(cell, smoothed_tree, speed_path, start_node, end_node);			//Replace the segment in the tree with the smoothed path

		SmoothSegments(cell, smoothed_tree, end_node);							//Call SmoothSegments on the next segment
	}
}

/* This function takes takes in start_node and end_node and a path in between and returns the smoothed path */
MicrogliaRegionTracer::PathType::Pointer MicrogliaRegionTracer::SmoothPath(Cell* cell, Tree* smoothed_tree, Node* start_node, Node* end_node, PathType::Pointer path )
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
	optimizer->SetMaximumStepLength( 0.1 );
	optimizer->SetMinimumStepLength( 0.1 );
	optimizer->SetRelaxationFactor( 0.00 );

	// Create path filter
	SpeedPathFilterType::Pointer pathFilter = SpeedPathFilterType::New();
	pathFilter->SetInput( cell->speed_image );
	pathFilter->SetCostFunction( cost );
	pathFilter->SetOptimizer( optimizer );
	pathFilter->SetTerminationValue( 10.0 );

	// Add path information
	SpeedPathFilterType::PathInfo info;
	SpeedPathFilterType::PointType start;
	start[0] = start_node->x; start[1] = start_node->y; start[2] = start_node->z;
	
	//std::cout << "Setting smoothing start index to: " << start << std::endl;
	

	const PathType::VertexListType *vertex_list = path->GetVertexList();

	for (unsigned int i = 1; i < vertex_list->Size() - 1; ++i)
	{
		PathType::ContinuousIndexType vertex = vertex_list->GetElement(i);
		//std::cout << "Adding smoothing waypoint index to: " << vertex << std::endl;
		//info.AddWayPoint(vertex);
	}

	SpeedPathFilterType::PointType end;
	end[0] = end_node->x; end[1] = end_node->y; end[2] = end_node->z;
	//std::cout << "Adding smoothing endpoint index to: " << end << std::endl;
		
	info.SetStartPoint( end );
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
}

/* This function modifies the tree passed to it and replaces the path between start_node and end_node with the speed_path */
void MicrogliaRegionTracer::ReplaceTreeSegmentWithPath(Cell* cell, Tree* smoothed_tree, PathType::Pointer speed_path, Node* start_node, Node* end_node)
{
	//Convert path output back to tree format
	if ( speed_path->GetVertexList()->Size() == 0 )
	{
		std::cerr << "WARNING: Path contains no points!" << std::endl;
		start_node->AddChild(end_node);
		end_node->SetParent(start_node);
		return;
	}

	const PathType::VertexListType *smoothed_vertex_list = speed_path->GetVertexList();

	Node* current_node = start_node;
	for (unsigned int i = 0; i < smoothed_vertex_list->Size(); ++i)
	{
		PathType::ContinuousIndexType path_index = smoothed_vertex_list->GetElement(i);

		Node* new_node = new Node(path_index[0], path_index[1], path_index[2], cell->next_available_ID );

		//std::cout << "Adding smoothed path index: " << path_index << std::endl;

		//Update next available ID
		++(cell->next_available_ID);

		//Connect this new node to the current_node
		current_node->AddChild(new_node);
		new_node->SetParent(current_node);
		smoothed_tree->AddNode(new_node, current_node);

		current_node = new_node;
	}
	//std::cout << "Connecting last index: " << end_node->x << " " << end_node->y << " " << end_node->z << std::endl;
	current_node->AddChild(end_node);
	end_node->SetParent(current_node);
}

/* The speed image for SmoothPath is created here */
void MicrogliaRegionTracer::CreateSpeedImage(Cell* cell)
{
	//Normalize vesselness image
	typedef itk::RescaleIntensityImageFilter< DistanceImageType > RescaleIntensityFilterType;
	RescaleIntensityFilterType::Pointer rescale_filter = RescaleIntensityFilterType::New();
	rescale_filter->SetOutputMinimum(0.0);
	rescale_filter->SetOutputMaximum(1.0);
	rescale_filter->SetInput(cell->vesselness_image);

	try
	{
		rescale_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "rescale_filter exception: " << err << std::endl;
		return;
	}

	//cell->speed_image = rescale_filter->GetOutput();

	//Raise normalized vesselness_image to third power
	typedef itk::PowImageFilter< DistanceImageType > PowImageFilterType;
	PowImageFilterType::Pointer pow_image_filter = PowImageFilterType::New();
	pow_image_filter->SetInput1(rescale_filter->GetOutput());
	pow_image_filter->SetConstant2( 0.33 );

	try
	{
		pow_image_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "pow_image_filter exception: " << err << std::endl;
		return;
	}

	cell->speed_image = pow_image_filter->GetOutput();

#if PRINT_ALL_IMAGES
	cell->WriteImage("speed_image.mhd", cell->speed_image);
#endif
}

/* Calculate the Euclidean Distance */
double MicrogliaRegionTracer::CalculateEuclideanDistance(ImageType::IndexType node1, ImageType::IndexType node2)
{
	ImageType::IndexType trace_vector;
	trace_vector[0] = node2[0] - node1[0];
	trace_vector[1] = node2[1] - node1[1];
	trace_vector[2] = node2[2] - node1[2];

	return sqrt(pow(trace_vector[0], 2.0) + pow(trace_vector[1], 2.0) + pow(trace_vector[2], 2.0));
}

/* This function smooths the tree from the root node to the leaf node in one go without taking into account the in between critical points as opposed to the previous method of smoothing between branch points */
Tree* MicrogliaRegionTracer::SmoothSegments2(Cell* cell, Tree* tree)
{
	Tree* smoothed_tree = new Tree();
	
	Node* root_node = tree->GetRoot();
	itk::Index<3> root_node_index;
	root_node_index[0] = root_node->x;
	root_node_index[1] = root_node->y;
	root_node_index[2] = root_node->z;

	smoothed_tree->SetRoot(root_node);

	std::vector<Node*> leaf_nodes;	
	tree->GetLeafNodes(leaf_nodes);

	//Smooth all the paths that is connected to this start_node
	std::vector< Node* >::iterator leaf_nodes_iter;
	for (leaf_nodes_iter = leaf_nodes.begin(); leaf_nodes_iter != leaf_nodes.end(); ++leaf_nodes_iter)
	{
		//Make a new path for each possible segment from this start_node
		PathType::Pointer path = PathType::New();
		path->Initialize();
		path->AddVertex(root_node_index);
		
		Node* leaf_node = *leaf_nodes_iter;
		itk::Index<3> leaf_node_index;
		leaf_node_index[0] = leaf_node->x;
		leaf_node_index[1] = leaf_node->y;
		leaf_node_index[2] = leaf_node->z;
		
		path->AddVertex(leaf_node_index);

		PathType::Pointer speed_path = SmoothPath(cell, tree, root_node, leaf_node, path);

		ReplaceTreeSegmentWithPath(cell, tree, speed_path, root_node, leaf_node);
	}

	return smoothed_tree;
}