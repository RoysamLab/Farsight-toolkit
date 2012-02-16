#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename)
{
	roi_grabber = new ROIGrabber(joint_transforms_filename, img_path, anchor_filename);
}

MicrogliaRegionTracer::ImageType::Pointer MicrogliaRegionTracer::GetMaskedImage(MaskedImageType::Pointer mask, ImageType::Pointer image)
{
	typedef itk::MaskNegatedImageFilter< ImageType, MaskedImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetMaskImage(mask);
	maskFilter->SetInput(image);
	try
	{
		maskFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "maskFilter Exception: " << err << std::endl;
	}

	return maskFilter->GetOutput();
}

MicrogliaRegionTracer::ImageType::Pointer MicrogliaRegionTracer::GetMaskedImage(std::string filename, ImageType::Pointer image)
{
	typedef itk::ImageFileReader< MaskedImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}

	return GetMaskedImage(reader->GetOutput(), image);
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
		cell->image = roi_grabber->GetROI(cell, roi_size, shift_index);

		//Need to do this to make the MaskImageFilter work since it expects similar origins
		ImageType::PointType origin = cell->image->GetOrigin();
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
		WriteImage(cell_filename_stream.str(), cell->image);

		roi_size = cell->image->GetLargestPossibleRegion().GetSize();	//The size of the returned image may not be the size of the image that entered because clipping at edges
		cell->setSize(roi_size);

		////Get the masked image
		//cell->masked_image = GetMaskedImage(soma_filename, cell->image);

		////Make the file name of the masked cell image
		//std::stringstream masked_cell_filename_stream;
		//masked_cell_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_masked.TIF";	//X_Y_Z_masked.TIF
		//
		////Write the masked cell image
		//WriteImage(masked_cell_filename_stream.str(), cell->masked_image);

		cells.push_back(cell);
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, itk::Image< unsigned char, 3>::Pointer image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned char, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, itk::Image< unsigned short, 3>::Pointer image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned short, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, itk::Image< float , 3 >::Pointer image)
{
	typedef itk::ImageFileWriter< itk::Image< float, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::Trace()
{
	std::vector<Cell*>::iterator cells_iter;
	
	//Trace cell by cell
	#pragma omp parallel for
	for (int k = 0; k < cells.size(); k++)
	//for (cells_iter = cells.begin(); cells_iter != cells.end(); cells_iter++)
	{
		//Cell* cell = *cells_iter;
		Cell* cell = cells[k];
		
		std::cout << "Calculating candidate pixels for a new cell" << std::endl;
		CalculateCandidatePixels(cell);

		std::cout << "Detected " << cell->critical_points_queue.size() << " critical points" << std::endl;
		std::cout << "Tree Building" << std::endl;
		BuildTree(cell);

		std::cout << "Kappa Sigma Thresholding" << std::endl;
		KappaSigmaThreshold(cell);
		HuangThreshold(cell);
		IntermodesThreshold(cell);
		IsoDataThreshold(cell);
		MaximumEntropyThreshold(cell);
		MomentsThreshold(cell);
		OtsuThreshold(cell);
		OtsuMultipleThreshold(cell);
		RenyiEntropyThreshold(cell);
		ShanbhagThreshold(cell);
		YenThreshold(cell);

		//SmoothPath(cell);
		delete cell;
	}
}

void MicrogliaRegionTracer::CalculateCandidatePixels(Cell* cell)
{
	LoG *log_obj = new LoG();
	//Calculate the LoG on multiple scales and put them into a vector
	std::cout << "Calculating Multiscale LoG" << std::endl;
	cell->multiscale_LoG_image = log_obj->RunMultiScaleLoG(cell);

	std::cout << "Starting ridge detection" << std::endl;
	RidgeDetection(cell);
}

void MicrogliaRegionTracer::RidgeDetection( Cell* cell )
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

	itk::ImageRegionIterator<ImageType> critical_point_img_iter(cell->critical_point_image, cell->critical_point_image->GetLargestPossibleRegion());
	itk::ImageRegionIterator<VesselnessImageType> vesselness_image_iter(cell->vesselness_image, cell->vesselness_image->GetLargestPossibleRegion());
	
	itk::Size<3> rad = {{1,1,1}};
	itk::Size<3> local_rad = {{1,1,1}};

	itk::NeighborhoodIterator<LoGImageType> neighbor_iter(rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());
	itk::NeighborhoodIterator<LoGImageType> local_maxima_iter(local_rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());

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
		vesselness_image_iter.Set(vesselness_score);
		
		if (vesselness_image_iter.Get() > 0.025 && local_maximum)
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

	WriteImage(criticalpointsFileNameStream.str(), cell->critical_point_image);
	WriteImage(vesselnessFileNameLocalStream.str(), cell->vesselness_image);
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

	double** AdjGraph = BuildAdjacencyGraph(cell);

	Tree* tree = BuildMST1(cell, AdjGraph);

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

	#pragma omp parallel for
	for (itk::int64_t k = 0; k < cell->critical_points_queue.size(); k++)
	{
		for (itk::uint64_t l = 0; l < cell->critical_points_queue.size(); l++)
		{
			AdjGraph[k][l] = CalculateDistance(k, l, cell);
		}
	}
	return AdjGraph;
}

double MicrogliaRegionTracer::CalculateDistance(itk::uint64_t k, itk::uint64_t l, Cell* cell)
{
	ImageType::IndexType node1 = cell->critical_points_queue[k];
	ImageType::IndexType node2 = cell->critical_points_queue[l];

	itk::int64_t distance_x = node1[0] - node2[0];
	itk::int64_t distance_y = node1[1] - node2[1];
	itk::int64_t distance_z = node1[2] - node2[2];

	return sqrt(pow(distance_x, 2.0) + pow(distance_y, 2.0) + pow(distance_z, 2.0));
}

Tree* MicrogliaRegionTracer::BuildMST1(Cell* cell, double** AdjGraph)
{	
	Tree* tree = new Tree();
	ImageType::IndexType root_index = cell->critical_points_queue.front();
	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], 0));
	
	std::ofstream traceFile;
	std::ofstream traceFile_local;
	std::ostringstream traceFileNameStream;
	std::ostringstream traceFileNameLocalStream;
	traceFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".swc";
	traceFileNameLocalStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_local.swc";
	traceFile.open(traceFileNameStream.str().c_str());
	traceFile_local.open(traceFileNameLocalStream.str().c_str());
	std::cout << "Opening " << traceFileNameStream.str() << std::endl;

	//Root node should have infinite weights to connect to
	for (int m = 0; m < cell->critical_points_queue.size(); m++)
		AdjGraph[m][0] = std::numeric_limits<double>::max();

	//Write out the root point
	traceFile << "1 3 " << cell->getX() << " " << cell->getY() << " " << cell->getZ() << " 1 -1" << std::endl;	//global coordinates
	traceFile_local << "1 3 " << cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0] << " " << cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1] << " " << cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2] << " 1 -1" << std::endl;	//local coordinates

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
				if (AdjGraph[node->getID()][k] < minimum_node_distance) //from l (connected point) to k (unconnected point) if the current distance is less than the minimum distance
				{
					minimum_parent_node = node;
					minimum_node_index_from_id = node->getID();
					minimum_node_index_to_id = k;
					minimum_node_distance = AdjGraph[node->getID()][k];
				}
			}	
		}

		
		
		//std::cout << "Found new edge from " << minimum_node_index_from_id << " to " << minimum_node_index_to_id << " Location: " << cell->critical_points_queue[minimum_node_index_from_id][0] << " " << cell->critical_points_queue[minimum_node_index_from_id][1] << " " << cell->critical_points_queue[minimum_node_index_from_id][2] << " " << cell->critical_points_queue[minimum_node_index_to_id][0] << " " << cell->critical_points_queue[minimum_node_index_to_id][1] << " "  << cell->critical_points_queue[minimum_node_index_to_id][2] << std::endl;
		ImageType::IndexType cell_origin;
		cell_origin[0] = cell->critical_points_queue[minimum_node_index_to_id][0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
		cell_origin[1] = cell->critical_points_queue[minimum_node_index_to_id][1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
		cell_origin[2] = cell->critical_points_queue[minimum_node_index_to_id][2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];
		
		ImageType::IndexType cell_origin_local;
		cell_origin_local[0] = cell->critical_points_queue[minimum_node_index_to_id][0];
		cell_origin_local[1] = cell->critical_points_queue[minimum_node_index_to_id][1];
		cell_origin_local[2] = cell->critical_points_queue[minimum_node_index_to_id][2];

		traceFile << minimum_node_index_to_id + 1 << " 3 " << cell_origin[0] << " " << cell_origin[1] << " "  << cell_origin[2] << " 1 " << minimum_node_index_from_id + 1 << std::endl;

		traceFile_local << minimum_node_index_to_id + 1 << " 3 " << cell_origin_local[0] << " " << cell_origin_local[1] << " "  << cell_origin_local[2] << " 1 " << minimum_node_index_from_id + 1 << std::endl;

		//Set the weight of the AdjGraph entries for this minimum edge to infinite so it is not visited again
		for (int m = 0; m < cell->critical_points_queue.size(); m++)
			AdjGraph[m][minimum_node_index_to_id] = std::numeric_limits<double>::max();
		AdjGraph[minimum_node_index_to_id][minimum_node_index_from_id] = std::numeric_limits<double>::max(); //So the parent doesnt become the child of its child

		//Connect the unconnected point
		Node* new_connected_node = new Node(cell->critical_points_queue[minimum_node_index_to_id][0], cell->critical_points_queue[minimum_node_index_to_id][1], cell->critical_points_queue[minimum_node_index_to_id][2], minimum_node_index_to_id);
		new_connected_node->SetParent(minimum_parent_node);
		tree->AddNode(new_connected_node, minimum_parent_node);
	}

	std::cout << cell->getRequestedSize()[0] << " " << cell->getRequestedSize()[1] << " " << cell->getRequestedSize()[2] << std::endl;
	std::cout << "Shift: " << cell->getShiftIndex()[0] << " " << cell->getShiftIndex()[1] << " " << cell->getShiftIndex()[2] << std::endl;
	std::cout << "Closing " << traceFileNameStream.str() << std::endl;
	traceFile.close();
	traceFile_local.close();
	return tree;
}

void MicrogliaRegionTracer::SmoothPath(Cell* cell)
{
	typedef itk::DanielssonDistanceMapImageFilter< LoGImageType, DistanceImageType, VoronoiImageType > DistanceMapFilterType;
	DistanceMapFilterType::Pointer dist_map_filter = DistanceMapFilterType::New();
	dist_map_filter->SetInput(cell->multiscale_LoG_image);
	try
	{
		dist_map_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{                                             
		std::cerr << "dist_map_filter exception at " << __FILE__ << " " << __LINE__ << std::endl;
		std::cerr << &err << std::endl;
	}

	VoronoiImageType::Pointer voronoi_image = dist_map_filter->GetVoronoiMap();
	DistanceImageType::Pointer dist_image = dist_map_filter->GetDistanceMap();

	std::ostringstream voronoi_image_filename_stream;
	std::ostringstream dist_image_filename_stream;
	voronoi_image_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_voronoi.mhd";
	dist_image_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_dist.mhd";

	WriteImage(voronoi_image_filename_stream.str(), voronoi_image);
	WriteImage(dist_image_filename_stream.str(), dist_image);
}

void MicrogliaRegionTracer::KappaSigmaThreshold(Cell* cell)
{
	typedef itk::KappaSigmaThresholdImageFilter< ImageType > KappaSigmaFilterType;
	KappaSigmaFilterType::Pointer kappa_sigma_filter = KappaSigmaFilterType::New();
	kappa_sigma_filter->SetInput(cell->image);
	kappa_sigma_filter->SetSigmaFactor(2.0);
	kappa_sigma_filter->SetNumberOfIterations(2);
	kappa_sigma_filter->Update();
	
	std::ostringstream kappa_sigma_filename_stream;
	kappa_sigma_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_kappa_sigma.mhd";
	WriteImage(kappa_sigma_filename_stream.str(), kappa_sigma_filter->GetOutput());
}

void MicrogliaRegionTracer::HuangThreshold(Cell* cell)
{
	typedef itk::HuangThresholdImageFilter< ImageType, ImageType > HuangThresholdImageFilterType;
	HuangThresholdImageFilterType::Pointer huang_threshold_filter = HuangThresholdImageFilterType::New();
	huang_threshold_filter->SetInput(cell->image);
	huang_threshold_filter->Update();

	std::ostringstream huang_filename_stream;
	huang_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_huang.mhd";
	WriteImage(huang_filename_stream.str(), huang_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::IntermodesThreshold(Cell* cell)
{
	typedef itk::IntermodesThresholdImageFilter< ImageType, ImageType > IntermodesThresholdImageFilterType;
	IntermodesThresholdImageFilterType::Pointer Intermodes_threshold_filter = IntermodesThresholdImageFilterType::New();
	Intermodes_threshold_filter->SetInput(cell->image);
	Intermodes_threshold_filter->Update();

	std::ostringstream Intermodes_filename_stream;
	Intermodes_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Intermodes.mhd";
	WriteImage(Intermodes_filename_stream.str(), Intermodes_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::IsoDataThreshold(Cell* cell)
{
	typedef itk::IsoDataThresholdImageFilter< ImageType, ImageType > IsoDataThresholdImageFilterType;
	IsoDataThresholdImageFilterType::Pointer IsoData_threshold_filter = IsoDataThresholdImageFilterType::New();
	IsoData_threshold_filter->SetInput(cell->image);
	IsoData_threshold_filter->Update();

	std::ostringstream IsoData_filename_stream;
	IsoData_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_IsoData.mhd";
	WriteImage(IsoData_filename_stream.str(), IsoData_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::MaximumEntropyThreshold(Cell* cell)
{
	typedef itk::MaximumEntropyThresholdImageFilter< ImageType, ImageType > MaximumEntropyThresholdImageFilterType;
	MaximumEntropyThresholdImageFilterType::Pointer MaximumEntropy_threshold_filter = MaximumEntropyThresholdImageFilterType::New();
	MaximumEntropy_threshold_filter->SetInput(cell->image);
	MaximumEntropy_threshold_filter->Update();

	std::ostringstream MaximumEntropy_filename_stream;
	MaximumEntropy_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_MaximumEntropy.mhd";
	WriteImage(MaximumEntropy_filename_stream.str(), MaximumEntropy_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::MomentsThreshold(Cell* cell)
{
	typedef itk::MomentsThresholdImageFilter< ImageType, ImageType > MomentsThresholdImageFilterType;
	MomentsThresholdImageFilterType::Pointer Moments_threshold_filter = MomentsThresholdImageFilterType::New();
	Moments_threshold_filter->SetInput(cell->image);
	Moments_threshold_filter->Update();

	std::ostringstream Moments_filename_stream;
	Moments_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Moments.mhd";
	WriteImage(Moments_filename_stream.str(), Moments_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::OtsuThreshold(Cell* cell)
{
	typedef itk::OtsuThresholdImageFilter< ImageType, ImageType > OtsuThresholdImageFilterType;
	OtsuThresholdImageFilterType::Pointer Otsu_threshold_filter = OtsuThresholdImageFilterType::New();
	Otsu_threshold_filter->SetInput(cell->image);
	Otsu_threshold_filter->Update();

	std::ostringstream Otsu_filename_stream;
	Otsu_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Otsu.mhd";
	WriteImage(Otsu_filename_stream.str(), Otsu_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::OtsuMultipleThreshold(Cell* cell)
{
	typedef itk::OtsuMultipleThresholdsImageFilter< ImageType, ImageType > OtsuMultipleThresholdImageFilterType;
	OtsuMultipleThresholdImageFilterType::Pointer Otsu_Multiple_thresholds_filter = OtsuMultipleThresholdImageFilterType::New();
	Otsu_Multiple_thresholds_filter->SetInput(cell->image);
	Otsu_Multiple_thresholds_filter->Update();

	std::ostringstream Otsu_Multiple_filename_stream;
	Otsu_Multiple_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Otsu_Multiple.mhd";
	WriteImage(Otsu_Multiple_filename_stream.str(), Otsu_Multiple_thresholds_filter->GetOutput());
}

void MicrogliaRegionTracer::RenyiEntropyThreshold(Cell* cell)
{
	typedef itk::RenyiEntropyThresholdImageFilter< ImageType, ImageType > RenyiEntropyThresholdImageFilterType;
	RenyiEntropyThresholdImageFilterType::Pointer RenyiEntropy_threshold_filter = RenyiEntropyThresholdImageFilterType::New();
	RenyiEntropy_threshold_filter->SetInput(cell->image);
	RenyiEntropy_threshold_filter->Update();

	std::ostringstream RenyiEntropy_filename_stream;
	RenyiEntropy_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_RenyiEntropy.mhd";
	WriteImage(RenyiEntropy_filename_stream.str(), RenyiEntropy_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::ShanbhagThreshold(Cell* cell)
{
	typedef itk::ShanbhagThresholdImageFilter< ImageType, ImageType > ShanbhagThresholdImageFilterType;
	ShanbhagThresholdImageFilterType::Pointer Shanbhag_threshold_filter = ShanbhagThresholdImageFilterType::New();
	Shanbhag_threshold_filter->SetInput(cell->image);
	Shanbhag_threshold_filter->Update();

	std::ostringstream Shanbhag_filename_stream;
	Shanbhag_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Shanbhag.mhd";
	WriteImage(Shanbhag_filename_stream.str(), Shanbhag_threshold_filter->GetOutput());
}

void MicrogliaRegionTracer::YenThreshold(Cell* cell)
{
	typedef itk::YenThresholdImageFilter< ImageType, ImageType > YenThresholdImageFilterType;
	YenThresholdImageFilterType::Pointer Yen_threshold_filter = YenThresholdImageFilterType::New();
	Yen_threshold_filter->SetInput(cell->image);
	Yen_threshold_filter->Update();

	std::ostringstream Yen_filename_stream;
	Yen_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_Yen.mhd";
	WriteImage(Yen_filename_stream.str(), Yen_threshold_filter->GetOutput());
}



