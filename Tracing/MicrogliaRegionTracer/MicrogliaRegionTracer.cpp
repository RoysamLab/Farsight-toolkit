#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename)
{
	roi_grabber = new ROIGrabber(joint_transforms_filename, img_path, anchor_filename);
}

void MicrogliaRegionTracer::LoadCellPoints(std::string seedpoints_filename, std::string soma_filename)
{
	std::ifstream seed_point_file;
	seed_point_file.open(seedpoints_filename.c_str());

	while (!seed_point_file.eof())
	{
		
		itk::uint64_t cellX, cellY, cellZ;
		seed_point_file >> cellX >> cellY >> cellZ;
		//std::cout << "Reading in cell: (" << cellX << ", " << cellY << ", " << cellZ << ")" << std::endl;
		
		Cell* cell = new Cell(cellX, cellY, cellZ);

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

		//Duplicate the image so that the requests don't propagate back up to roi_grabber
		typedef itk::ImageDuplicator< ImageType > DuplicatorType;
		DuplicatorType::Pointer duplicator = DuplicatorType::New();
		duplicator->SetInputImage(temp_cell_image);
		duplicator->Update();

		cell->image = duplicator->GetOutput();

		//Need to do this to make the MaskImageFilter work since it expects similar origins
		ImageType::PointType origin = cell->image->GetOrigin();
		cell->SetOrigin(origin);
		origin[0] = 0;
		origin[1] = 0;
		origin[2] = 0;
		cell->image->SetOrigin(origin);


		cell->setShiftIndex(shift_index);
		//std::cout << cell->getShiftIndex()[0] << " " << cell->getShiftIndex()[1] << " " << cell->getShiftIndex()[2] << std::endl;
		
		//Make the file name of the raw cell image
		std::stringstream cell_filename_stream;
		cell_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".TIF";	//X_Y_Z.TIF

		//Write the cell image
		Cell::WriteImage(cell_filename_stream.str(), cell->image);

		roi_size = cell->image->GetLargestPossibleRegion().GetSize();	//The size of the returned image may not be the size of the image that entered because clipping at edges
		cell->SetSize(roi_size);

		//Get the mask
		cell->GetMask(soma_filename);

		//Get the masked image
		cell->ComputeMaskedImage();

		cell->masked_image = cell->image;

		cells.push_back(cell);
	}
}

void MicrogliaRegionTracer::Trace()
{
	//Trace cell by cell
	#pragma omp parallel for
	for (int k = 0; k < cells.size(); k++)
	{
		Cell* cell = cells[k];
		
		std::cout << "Calculating candidate pixels for a new cell" << std::endl;
		CalculateCandidatePixels(cell);

		std::cout << "Detected " << cell->critical_points_queue.size() << " critical points" << std::endl;
		
		std::cout << "Tree Building" << std::endl;
		BuildTree(cell);

		delete cell;
	}
}

void MicrogliaRegionTracer::Trace2()
{
	//Trace cell by cell
	//#pragma omp parallel for
	for (int k = 0; k < cells.size(); k++)
	{	
		Cell* cell = cells[k];

		std::cout << "Thresholding" << std::endl;
		cell->MaximumEntropyThreshold();

		std::cout << "Skeletonize" << std::endl;
		cell->Skeletonize();
		cell->critical_point_image = cell->skeleton_image;
  
		cell->ComputeCriticalPointsVector(cell->skeleton_image);

		std::cout << "Tracing based on skeleton" << std::endl;
		TraceSkeletonImage(cell);

		delete cell;
	}
}

void MicrogliaRegionTracer::CalculateCandidatePixels(Cell* cell)
{
	std::cout << "Ridge Detection" << std::endl;
	RidgeDetection(cell);
	
	std::cout << "VesselnessDetection" << std::endl;
	VesselnessDetection(cell);
}

void MicrogliaRegionTracer::RidgeDetection( Cell* cell )
{
	itk::Index<3> seed_index = {{cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0], cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1], cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2] }};
	cell->critical_points_queue.push_front(seed_index);
	
	LoG *log_obj = new LoG();
	//Calculate the LoG on multiple scales and store into an iamge
	std::cout << "Calculating Multiscale LoG" << std::endl;
	cell->multiscale_LoG_image = log_obj->RunMultiScaleLoG(cell);
	
	//Make a new image to store the critical points	
	ImageType::SizeType size = cell->multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	cell->critical_point_image = ImageType::New();
	ImageType::IndexType start;
	start.Fill(0);

	ImageType::RegionType region(start, size);
	cell->critical_point_image->SetRegions(region);
	cell->critical_point_image->Allocate();
	cell->critical_point_image->FillBuffer(0);

	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell->critical_point_image, cell->critical_point_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > LoG_neighbor_iter(rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());

	while(!LoG_neighbor_iter.IsAtEnd()) 
	{
		unsigned int neighborhood_size = (rad[0] * 2 + 1) * (rad[1] * 2 + 1) * (rad[2] * 2 + 1);
		unsigned int center_pixel_offset_index = neighborhood_size / 2;

		LoGImageType::PixelType center_pixel_intensity = LoG_neighbor_iter.GetPixel(center_pixel_offset_index);
		
		if (center_pixel_intensity >= 0.02)	//Must have greater than 2% response from LoG to even be considered for local maximum
		{	
			bool local_maximum = true;
			
			for (unsigned int neighborhood_index = 0; neighborhood_index < neighborhood_size; ++neighborhood_index)
			{
				if (neighborhood_index != center_pixel_offset_index)
				{
					bool isInBounds; //true if the pixel is in bounds
					LoGImageType::PixelType neighbor_pixel_intensity = LoG_neighbor_iter.GetPixel(neighborhood_index, isInBounds);

					if (isInBounds && neighbor_pixel_intensity > center_pixel_intensity)
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
		++LoG_neighbor_iter;
	}

	std::ostringstream criticalpointsFileNameStream;
	criticalpointsFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_critical.mhd";
	//Cell::WriteImage(criticalpointsFileNameStream.str(), cell->critical_point_image);
}

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
	multiscale_hessian_filter->SetInput(cell->image);
	multiscale_hessian_filter->SetSigmaStepMethodToEquispaced();
	multiscale_hessian_filter->SetSigmaMinimum(2.0);
	multiscale_hessian_filter->SetSigmaMaximum(6.0);
	multiscale_hessian_filter->SetNumberOfSigmaSteps(10);
	multiscale_hessian_filter->SetNonNegativeHessianBasedMeasure(true);
	multiscale_hessian_filter->SetHessianToMeasureFilter(vesselness_filter);

	try
	{
		multiscale_hessian_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "multiscale_hessian_filter exception: " << &err << std::endl;
	}

	cell->vesselness_image = multiscale_hessian_filter->GetOutput();
	
	std::ostringstream vesselnessFileNameStream;

	vesselnessFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_vesselness.mhd";
	//Cell::WriteImage(vesselnessFileNameLocalStream.str(), cell->vesselness_image);
}

void MicrogliaRegionTracer::RidgeDetection2( Cell* cell )
{
	ImageType::SizeType size = cell->image->GetLargestPossibleRegion().GetSize();

	//Make a new image to store the critical points	
	cell->critical_point_image = ImageType::New();
	ImageType::IndexType start;
	start.Fill(0);

	ImageType::RegionType region(start, size);
	cell->critical_point_image->SetRegions(region);
	cell->critical_point_image->Allocate();
	cell->critical_point_image->FillBuffer(0);

	//Make a new image to store the vesselness scores	
	cell->vesselness_image = VesselnessImageType::New();
	VesselnessImageType::IndexType vesselness_start;
	vesselness_start.Fill(0);

	VesselnessImageType::RegionType vesselness_region(vesselness_start, size);
	cell->vesselness_image->SetRegions(vesselness_region);
	cell->vesselness_image->Allocate();
	cell->vesselness_image->FillBuffer(0);

	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell->critical_point_image, cell->critical_point_image->GetLargestPossibleRegion());
	itk::ImageRegionIterator< VesselnessImageType > vesselness_image_iter(cell->vesselness_image, cell->vesselness_image->GetLargestPossibleRegion());

	itk::Size<3> rad = {{1,1,1}};
	itk::Size<3> local_rad = {{1,1,1}};

	itk::NeighborhoodIterator< LoGImageType > neighbor_iter(rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());
	itk::NeighborhoodIterator< LoGImageType > local_maxima_iter(local_rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());

	//Making a size value that removes some room (kind of like reverse padding) so we dont go out of bounds later
	itk::Size<3> unpad_size = cell->multiscale_LoG_image->GetLargestPossibleRegion().GetSize();
	unpad_size[0] = unpad_size[0] - (rad[0] + 1);
	unpad_size[1] = unpad_size[1] - (rad[1] + 1); 
	unpad_size[2] = unpad_size[2] - (rad[2] + 1);

	critical_point_img_iter.GoToBegin();
	neighbor_iter.GoToBegin();
	vesselness_image_iter.GoToBegin();
	local_maxima_iter.GoToBegin();

	//Need the maximum intensity of the multiscale LoG image for Vesselness
	typedef itk::MinimumMaximumImageCalculator< LoGImageType > ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(cell->multiscale_LoG_image);
	imageCalculatorFilter->ComputeMaximum();
	double max_intensity_multiscale = imageCalculatorFilter->GetMaximum();


	while(!neighbor_iter.IsAtEnd()) 
	{
		itk::Index<3> index = neighbor_iter.GetIndex();

		//checking to see if we are at the edge of the image
		if ( (index[0] < rad[0] + 1) || (index[1] < rad[1] + 1) || (index[2] < rad[2] + 1) ||								
			(index[0] >= (int)unpad_size[0]) || (index[1] >= (int)unpad_size[1]) || (index[2] >= (int)unpad_size[2]) )		
		{
			++critical_point_img_iter;
			++neighbor_iter;
			++vesselness_image_iter;
			++local_maxima_iter;
			continue;
		}

		//Non-maximal suppression (only keep the seedpoint with the maximum value in the 5x5x5 neighborhood)
		bool local_maximum = true;
		itk::uint64_t local_maximum_search_size = local_maxima_iter.GetSize(0) * local_maxima_iter.GetSize(1) * local_maxima_iter.GetSize(2);
		for (int i = 0; i < local_maximum_search_size; ++i)
		{
			if (local_maxima_iter.GetPixel(i) > local_maxima_iter.GetPixel(local_maximum_search_size / 2))
				local_maximum = false;
		}

		if (local_maximum)
		{		
			critical_point_img_iter.Set(255);
		}

		float vesselness_score = RunHessian(cell->multiscale_LoG_image, neighbor_iter, max_intensity_multiscale);

		if (vesselness_score > cell->vesselness_image_maximum_intensity)
			cell->vesselness_image_maximum_intensity = vesselness_score;

		vesselness_image_iter.Set(vesselness_score);

		if (vesselness_image_iter.Get() > 0.005 && local_maximum)
			cell->critical_points_queue.push_back(vesselness_image_iter.GetIndex());

		++local_maxima_iter;
		++critical_point_img_iter;
		++neighbor_iter;
		++vesselness_image_iter;
	}

	std::ostringstream criticalpointsFileNameStream;
	std::ostringstream vesselnessFileNameLocalStream;	

	criticalpointsFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_critical.mhd";
	vesselnessFileNameLocalStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_vesselness.mhd";

	Cell::WriteImage(criticalpointsFileNameStream.str(), cell->critical_point_image);
	Cell::WriteImage(vesselnessFileNameLocalStream.str(), cell->vesselness_image);
}

double MicrogliaRegionTracer::RunHessian( LoGImageType::Pointer log_image, itk::NeighborhoodIterator<LoGImageType> neighbor_iter, double max_intensity_multiscale )
{
	// set the diagonal terms in neighborhood iterator, this is the offsets for the diametrically opposing pixels
	itk::Offset<3>
		xp =  {{2 ,  0 ,   0}},
		xn =  {{-2,  0,    0}},
		yp =  {{0,   2,   0}},
		yn =  {{0,  -2,    0}},
		zp =  {{0,   0,    2}},
		zn =  {{0,   0,   -2}};

	//{0, 0, 0} is the center pixel of the 3x3x3 neighborhood. The constants are then the pixel index starting from the top-left corner of the front face.
	//x: left is -1
	//y: up is -1
	//z: out of page is -1
	unsigned int
		//{ x,    y ,  z }
		xy1 =  17, //{ 1 ,   1 ,  0 },
		xy2 =  9,  //{ -1,  -1 ,  0 },
		xy3 =  15, //{ -1,   1 ,  0 },
		xy4 =  11, //{ 1 ,  -1 ,  0 },

		yz1 =  25, //{ 0 ,   1 ,  1 },
		yz2 =  1,  //{ 0 ,  -1 , -1 },
		yz3 =  19, //{ 0 ,  -1 ,  1 },
		yz4 =  7,  //{ 0 ,   1 , -1 },

		xz1 =  23, //{ 1 ,   0 ,  1 },
		xz2 =  3,  //{-1 ,   0 , -1 },
		xz3 =  21, //{-1 ,   0 ,  1 },
		xz4 =  5;  //{ 1 ,   0 , -1 };
	
	typedef itk::SymmetricSecondRankTensor<double,3> TensorType;
	
	TensorType hessian;
	itk::Index<3> index = neighbor_iter.GetIndex();
	hessian[0] = log_image->GetPixel( index + xp ) +
		log_image->GetPixel( index + xn ) -
		2*neighbor_iter.GetPixel( 13 );
	hessian[3] = log_image->GetPixel( index + yp ) +
		log_image->GetPixel( index + yn ) -
		2*neighbor_iter.GetPixel( 13 );
	hessian[5] = log_image->GetPixel( index + zp ) +
		log_image->GetPixel( index + zn ) -
		2*neighbor_iter.GetPixel( 13 );
	hessian[1] = neighbor_iter.GetPixel(xy1) + neighbor_iter.GetPixel(xy2) -
		neighbor_iter.GetPixel(xy3) - neighbor_iter.GetPixel(xy4);
	hessian[2] = neighbor_iter.GetPixel(xz1) + neighbor_iter.GetPixel(xz2) -
		neighbor_iter.GetPixel(xz3) - neighbor_iter.GetPixel(xz4);
	hessian[4] = neighbor_iter.GetPixel(yz1) + neighbor_iter.GetPixel(yz2) -
		neighbor_iter.GetPixel(yz3) - neighbor_iter.GetPixel(yz4);

	typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
	typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;

	EigenValuesArrayType ev;
	EigenVectorMatrixType em;

	hessian.ComputeEigenAnalysis (ev, em);	//Compute Eigenvalues

	double vesselness_score = ComputeVesselness(ev[0], ev[1], ev[2], max_intensity_multiscale);
	
	if (vesselness_score < 0)
		std::cout << "negative vesselness score detected: " << vesselness_score << std::endl;

	return vesselness_score;
}

double MicrogliaRegionTracer::ComputeVesselness( double ev1, double ev2, double ev3, double maximum_intensity )
{
	double lambda1, lambda2, lambda3; //constrain lambda1 <= lambda2 <= lambda3

	double ev1_magnitude = std::abs(ev1);
	double ev2_magnitude = std::abs(ev2);
	double ev3_magnitude = std::abs(ev3);


	///std::cout << ev1_magnitude << " " << ev2_magnitude << " " << ev3_magnitude << std::endl;
	
	if (ev1_magnitude <= ev2_magnitude && ev1_magnitude <= ev3_magnitude) //ev1 is the smallest eigenvalue
	{
		lambda1 = ev1;
		if (ev2_magnitude < ev3_magnitude)
		{
			lambda2 = ev2;
			lambda3 = ev3;
		}
		else
		{
			lambda2 = ev3;
			lambda3 = ev2;
		}
	}
	else if (ev2_magnitude <= ev1_magnitude && ev2_magnitude <= ev3_magnitude) //ev2 is the smallest eigenvalue
	{
		lambda1 = ev2;
		if (ev1_magnitude < ev3_magnitude)
		{
			lambda2 = ev1;
			lambda3 = ev3;
		}
		else
		{
			lambda2 = ev3;
			lambda3 = ev1;
		}
	}
	else
	{
		lambda1 = ev3;
		if (ev1_magnitude < ev2_magnitude)
		{
			lambda2 = ev1;
			lambda3 = ev2;
		}
		else
		{
			lambda2 = ev2;
			lambda3 = ev1;
		}
	}


	double vesselness_score;
	if (lambda3 < 0 && lambda2 < 0)
	{
		lambda1 = std::abs(lambda1);
		lambda2 = std::abs(lambda2);
		lambda3 = std::abs(lambda3);

		double r_a = lambda1 / sqrt(lambda2 * lambda3);
		double r_b = lambda2 / lambda3;
		double s = sqrt(pow(lambda1, 2) + pow(lambda2, 2) + pow(lambda3, 2));

		//std::cout << r_a << " " << r_b << " " << s << std::endl;

		double alpha = 0.5 * maximum_intensity;
		double beta = 0.5 * maximum_intensity;
		double gamma = 0.25 * maximum_intensity;
		
		vesselness_score =			exp(-1 * pow(r_a, 2) / (2 * pow(alpha, 2))) 
								* (1 -	exp(-1 * pow(r_b, 2) / (2 * pow(beta, 2)))) 
								* (1 -	exp(-1 * pow(s, 2) / (2 * pow(gamma, 2))));

		//std::cout << vesselness_score << std::endl;
	}
	else
		vesselness_score = 0;

	return vesselness_score;
}


void MicrogliaRegionTracer::BuildTree(Cell* cell)
{
	ImageType::IndexType cell_index;
	cell_index[0] = cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0];
	cell_index[1] = cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1];
	cell_index[2] = cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2];

	cell->critical_points_queue.push_front(cell_index);	//Add the centroid to the critical points

	std::cout << "Building Adjacency Graph" << std::endl;
	double** AdjGraph = BuildAdjacencyGraph(cell);

	Tree* tree = BuildMST1(cell, AdjGraph);

	std::ostringstream swc_filename_stream, swc_filename_stream_local;
	swc_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree.swc";
	swc_filename_stream_local << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree_local.swc";


	WriteTreeToSWCFile(tree, cell, swc_filename_stream.str(), swc_filename_stream_local.str());

	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		delete[] AdjGraph[k];
	delete[] AdjGraph;

	delete tree;
}

double** MicrogliaRegionTracer::BuildAdjacencyGraph(Cell* cell)
{
	double** AdjGraph = new double*[cell->critical_points_queue.size()];
	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		AdjGraph[k] = new double[cell->critical_points_queue.size()];

	for (itk::uint64_t k = 0; k < cell->critical_points_queue.size(); k++)
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

	//If we are connecting from the root node
	if (node_from == 0)
	{	
		if (cell->soma_label_image->GetPixel(node1) == cell->soma_label_image->GetPixel(node2)) //If it is inside the soma our centroid belongs to, we give it 0 weight
			return 0;
		else
			return std::numeric_limits< double >::max();	
	}
	typedef itk::PolyLineParametricPath< 3 > PathType;
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
	
	//path_iter.GoToBegin();
	while (!path_iter.IsAtEnd())
	{
		//std::cout << "Path iterator position: " << path_iter.GetPathPosition() << std::endl;
		//std::cout << "Path iterator index: " << path_iter.GetIndex() << std::endl;
		if (path_iter.Get() < 1.0)
			++num_blank_pixels;
		++path_iter;
		++path_length;
	}
	
	if (num_blank_pixels >= 5)
		return 100000;
	else
		return mag_trace_vector;
}

Tree* MicrogliaRegionTracer::BuildMST1(Cell* cell, double** AdjGraph)
{	
	Tree* tree = new Tree();
	ImageType::IndexType root_index = cell->critical_points_queue.front();
	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], 1));
	
	//std::ofstream traceFile;
	//std::ofstream traceFile_local;
	//std::ostringstream traceFileNameStream;
	//std::ostringstream traceFileNameLocalStream;
	//traceFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".swc";
	//traceFileNameLocalStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_local.swc";
	//traceFile.open(traceFileNameStream.str().c_str());
	//traceFile_local.open(traceFileNameLocalStream.str().c_str());
	//std::cout << "Opening " << traceFileNameStream.str() << std::endl;

	//Root node should have infinite weights to connect to
	for (int m = 0; m < cell->critical_points_queue.size(); m++)
		AdjGraph[m][0] = std::numeric_limits<double>::max();

	//Write out the root point
	//traceFile << "1 3 " << cell->getX() << " " << cell->getY() << " " << cell->getZ() << " 1 -1" << std::endl;	//global coordinates
	//traceFile_local << "1 3 " << cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0] << " " << cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1] << " " << cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2] << " 1 -1" << std::endl;	//local coordinates

	//for each node but the last, find its child
	for (itk::uint64_t l = 0; l < cell->critical_points_queue.size() - 1; l++)
	{
		//std::cout << "Calculating nearest neighbor for point " << l << std::endl;
		itk::uint64_t minimum_node_index_from_id = 0;
		itk::uint64_t minimum_node_index_to_id = 0;
		double minimum_node_distance = std::numeric_limits<double>::max();
		Node* minimum_parent_node = NULL;
		
		//For each connected point
		std::vector<Node*> member_nodes = tree->GetMemberNodes();
		std::vector<Node*>::iterator member_nodes_iter;
		for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)
		{
			Node* node = *member_nodes_iter;

			//Search through all the nodes and find the minimum distance
			for (itk::uint64_t k = 0; k < cell->critical_points_queue.size(); k++)
			{
				if (AdjGraph[node->getID() - 1][k] < minimum_node_distance) //from l (connected point) to k (unconnected point) if the current distance is less than the minimum distance
				{
					minimum_parent_node = node;
					minimum_node_index_from_id = node->getID() - 1;
					minimum_node_index_to_id = k;
					minimum_node_distance = AdjGraph[minimum_node_index_from_id][minimum_node_index_to_id];
				}
			}	
		}

		if (minimum_node_distance >= 100000)
			break;	//Minimum distance way too far
		
		//std::cout << "Found new edge from " << minimum_node_index_from_id << " to " << minimum_node_index_to_id << " Location: " << cell->critical_points_queue[minimum_node_index_from_id][0] << " " << cell->critical_points_queue[minimum_node_index_from_id][1] << " " << cell->critical_points_queue[minimum_node_index_from_id][2] << " " << cell->critical_points_queue[minimum_node_index_to_id][0] << " " << cell->critical_points_queue[minimum_node_index_to_id][1] << " "  << cell->critical_points_queue[minimum_node_index_to_id][2] << std::endl;
		ImageType::IndexType cell_origin;
		cell_origin[0] = cell->critical_points_queue[minimum_node_index_to_id][0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
		cell_origin[1] = cell->critical_points_queue[minimum_node_index_to_id][1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
		cell_origin[2] = cell->critical_points_queue[minimum_node_index_to_id][2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];
		
		ImageType::IndexType cell_origin_local;
		cell_origin_local[0] = cell->critical_points_queue[minimum_node_index_to_id][0];
		cell_origin_local[1] = cell->critical_points_queue[minimum_node_index_to_id][1];
		cell_origin_local[2] = cell->critical_points_queue[minimum_node_index_to_id][2];

		//traceFile << minimum_node_index_to_id + 1 << " 3 " << cell_origin[0] << " " << cell_origin[1] << " "  << cell_origin[2] << " 1 " << minimum_parent_node->getID() << std::endl;

		//traceFile_local << minimum_node_index_to_id + 1 << " 3 " << cell_origin_local[0] << " " << cell_origin_local[1] << " "  << cell_origin_local[2] << " 1 " << minimum_parent_node->getID() << std::endl;

		//Set the weight of the AdjGraph entries for this minimum edge to infinite so it is not visited again
		for (int m = 0; m < cell->critical_points_queue.size(); m++)
			AdjGraph[m][minimum_node_index_to_id] = std::numeric_limits<double>::max();
		AdjGraph[minimum_node_index_to_id][minimum_node_index_from_id] = std::numeric_limits<double>::max(); //So the parent doesn't become the child of its child

		//Connect the unconnected point
		Node* new_connected_node = new Node(cell_origin_local[0], cell_origin_local[1], cell_origin_local[2], minimum_node_index_to_id + 1);
		new_connected_node->SetParent(minimum_parent_node);
		tree->AddNode(new_connected_node, minimum_parent_node);
	}

	std::cout << cell->getRequestedSize()[0] << " " << cell->getRequestedSize()[1] << " " << cell->getRequestedSize()[2] << std::endl;
	std::cout << "Shift: " << cell->getShiftIndex()[0] << " " << cell->getShiftIndex()[1] << " " << cell->getShiftIndex()[2] << std::endl;
	//std::cout << "Closing " << traceFileNameStream.str() << std::endl;
	//traceFile.close();
	//traceFile_local.close();
	return tree;
}



void MicrogliaRegionTracer::TraceSkeletonImage(Cell* cell)
{
	ImageType::IndexType root_index = FindNearestCriticalPointToCentroid(cell);

	//Make a new image to store the visited points	
	ImageType::Pointer visited_image = ImageType::New();
	ImageType::IndexType start;
	start.Fill(0);

	itk::Size<3> size = cell->thresholded_image->GetLargestPossibleRegion().GetSize();

	ImageType::RegionType region(start, size);
	visited_image->SetRegions(region);
	visited_image->Allocate();
	visited_image->FillBuffer(0);
		
	std::ofstream traceFile;
	std::ofstream traceFile_local;
	std::ostringstream traceFileNameStream;
	std::ostringstream traceFileNameLocalStream;
	traceFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".swc";
	traceFileNameLocalStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_local.swc";
	traceFile.open(traceFileNameStream.str().c_str());
	traceFile_local.open(traceFileNameLocalStream.str().c_str());
	std::cout << "Opening " << traceFileNameStream.str() << std::endl;

	itk::int64_t swc_line_number = 1;
	FollowSkeleton(cell, root_index, visited_image, -1, swc_line_number, traceFile, traceFile_local);

	traceFile.close();
	traceFile_local.close();
}

void MicrogliaRegionTracer::FollowSkeleton( Cell* cell, ImageType::IndexType index, ImageType::Pointer visited_image, itk::int64_t parent_id, itk::int64_t &swc_line_number, std::ofstream &traceFile, std::ofstream &traceFile_local)
{
	itk::int64_t id = swc_line_number;
	traceFile_local << id << " 3 " << index[0] << " " << index[1] << " "  << index[2] << " 1 " << parent_id << std::endl;
	itk::Size<3> rad = {{1,1,1}};

	itk::uint64_t neighborhood_size = (2 * rad[0] + 1) * (2 * rad[1] + 1) * (2 * rad[2] + 1);

	itk::NeighborhoodIterator< ImageType > neighborhood_iter(rad, cell->skeleton_image, cell->skeleton_image->GetLargestPossibleRegion().GetSize());
	neighborhood_iter.SetLocation(index);

	for (int k = 0; k < neighborhood_size; k++)
	{
		bool IsInBounds;
		ImageType::PixelType intensity = neighborhood_iter.GetPixel(k, IsInBounds);
		if (k != neighborhood_size/2 && IsInBounds)
		{
			ImageType::IndexType child_index = neighborhood_iter.GetIndex(k);

			if (visited_image->GetPixel(child_index) == 0)	//Not visited, so we can go on
			{
				if (intensity != 0)
				{
					//std::cout << " Index: " << index << " Child Index: " << child_index << std::endl;
					visited_image->SetPixel(index, std::numeric_limits< ImageType::PixelType >::max());
					FollowSkeleton(cell, child_index, visited_image, id, ++swc_line_number, traceFile, traceFile_local);
				}
			}
		}
	}
}

MicrogliaRegionTracer::ImageType::IndexType MicrogliaRegionTracer::FindNearestCriticalPointToCentroid(Cell* cell)
{
	ImageType::IndexType centroid_location_local, root_location;
	centroid_location_local[0] = cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0];
	centroid_location_local[1] = cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1];
	centroid_location_local[2] = cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2];

	double min_distance = std::numeric_limits< double >::max();

	std::deque< ImageType::IndexType >::iterator critical_points_queue_iter;

	for (critical_points_queue_iter = cell->critical_points_queue.begin(); critical_points_queue_iter != cell->critical_points_queue.end(); ++critical_points_queue_iter)
	{
		ImageType::IndexType critical_point_index = *critical_points_queue_iter;

		double distance = sqrt(
								pow(critical_point_index[0] - centroid_location_local[0], 2.0) +
								pow(critical_point_index[1] - centroid_location_local[1], 2.0) +
								pow(critical_point_index[2] - centroid_location_local[2], 2.0) ) ;

		if (distance < min_distance)
		{
			min_distance = distance;
			root_location = critical_point_index;
		}
	}

	return root_location;
}

void MicrogliaRegionTracer::WriteTreeToSWCFile(Tree* tree, Cell* cell, std::string filename, std::string filename_local)
{
	std::cout << "Entering WriteTreeToSWCFile" << std::endl;
	std::ofstream traceFile, traceFile_local;
	
	std::cout << "Opening " << filename << std::endl;
	traceFile.open(filename.c_str());
	std::cout << "Opening " << filename_local << std::endl;
	traceFile_local.open(filename_local.c_str());
	
	Node* root = tree->getRoot();

	itk::uint64_t tree_depth = 0; //root node is defined as tree depth 0
	WriteLinkToParent(root, tree_depth, cell, traceFile, traceFile_local);

	traceFile.close();
	traceFile_local.close();	
}

void MicrogliaRegionTracer::WriteLinkToParent(Node* node, itk::uint64_t tree_depth, Cell* cell, std::ofstream &traceFile, std::ofstream &traceFile_local)
{
	//Calculate some node indices
	ImageType::IndexType node_index, node_index_local;
	node_index_local[0] = node->x;
	node_index_local[1] = node->y;
	node_index_local[2] = node->z;
	node_index[0] = node_index_local[0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
	node_index[1] = node_index_local[1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
	node_index[2] = node_index_local[2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];

	itk::int64_t parent_node_id;
	if (node->GetParent() == NULL) //This node has no parent, so this is the root node, so parent ID is -1
		parent_node_id = -1;
	else
		parent_node_id = node->GetParent()->getID();

	std::vector< Node* > children = node->GetChildren();
	if (tree_depth == 1 && children.size() == 0)
		return; //Don't write out a trace if we are at depth one and have no children because we are a trace to the edge of the soma

	//Write out the SWC lines
	traceFile		<< node->getID() << " 3 " << node_index[0]			<< " " << node_index[1]			<< " "  << node_index[2]		<< " 1 " << parent_node_id << std::endl;
	traceFile_local << node->getID() << " 3 " << node_index_local[0]	<< " " << node_index_local[1]	<< " "  << node_index_local[2]	<< " 1 " << parent_node_id << std::endl;
	
	if (children.size() == 0) //No children, so don't visit them
		return;

	std::vector< Node* >::iterator children_iter;
	for (children_iter = children.begin(); children_iter != children.end(); ++children_iter)
	{
		Node* child_node = *children_iter;
		WriteLinkToParent(child_node, tree_depth+1, cell, traceFile, traceFile_local);
	}
}
