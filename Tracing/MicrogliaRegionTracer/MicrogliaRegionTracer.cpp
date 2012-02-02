#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename)
{
	roi_grabber = new ROIGrabber(joint_transforms_filename, img_path, anchor_filename);
}

void MicrogliaRegionTracer::LoadImage(ImageType::Pointer image)
{
	this->image = image;
}

void MicrogliaRegionTracer::LoadImage(std::string filename)
{
	typedef itk::ImageFileReader< ImageType > ReaderType;
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

	LoadImage(reader->GetOutput());
}

void MicrogliaRegionTracer::LoadCellPoints(std::string filename)
{
	std::ifstream seed_point_file;
	seed_point_file.open(filename.c_str());

	while (!seed_point_file.eof())
	{
		size_t cellX, cellY, cellZ;
		seed_point_file >> cellX >> cellY >> cellZ;

		//std::cout << "Reading in cell: (" << cellX << ", " << cellY << ", " << cellZ << ")" << std::endl;
		cells.push_back(new Cell(cellX, cellY, cellZ));
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, ImageType::Pointer image)
{
	typedef itk::ImageFileWriter< ImageType > WriterType;
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

void MicrogliaRegionTracer::WriteVesselnessImage(std::string filename, VesselnessImageType::Pointer image)
{
	typedef itk::ImageFileWriter< VesselnessImageType > WriterType;
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

void MicrogliaRegionTracer::WriteCellImages()
{
	std::vector<Cell*>::iterator cells_iter;

	if (cells.size() == 0)
	{
		std::cout << "There are no cells" << std::endl;
		exit(0);
	}
	//clock_t startTime = clock();
	for (cells_iter = cells.begin(); cells_iter != cells.end(); cells_iter++)
	{
		Cell* cell = *cells_iter;

		ImageType::SizeType roi_size;
		roi_size[0] = 200;
		roi_size[1] = 200;
		roi_size[2] = 100;
		
		//Grab the cell and its ROI
		ImageType::Pointer cell_image = roi_grabber->GetROI(cell, roi_size);

		//Make the file name of the raw cell image
		std::stringstream cell_filename_stream;
		cell_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << ".TIF";	//X_Y_Z.TIF
		
		//Write the cell image
		WriteImage(cell_filename_stream.str(), cell_image);
	}
	//std::cout << "Grabbed " << cells.size() << " cells in " << (clock() - startTime)/CLOCKS_PER_SEC << " seconds" << std::endl;
}

void MicrogliaRegionTracer::Trace()
{
	std::vector<Cell*>::iterator cells_iter;
	
	//Trace cell by cell
	for (cells_iter = cells.begin(); cells_iter != cells.end(); cells_iter++)
	{
		Cell* cell = *cells_iter;
		std::vector<ImageType::IndexType> critical_points_vector;
		
		std::cout << "Calculating candidate pixels for a new cell" << std::endl;
		CalculateCandidatePixels(cell, critical_points_vector);

		std::cout << "Tree Building" << std::endl;
		BuildTree(cell, critical_points_vector);
	}
}

void MicrogliaRegionTracer::CalculateCandidatePixels(Cell* cell, std::vector<ImageType::IndexType> &critical_points_vector)
{
	ImageType::SizeType roi_size;
	roi_size[0] = 200;
	roi_size[1] = 200;
	roi_size[2] = 100;
	
	//Grab the initial cellimage
	std::cout << "Grabbing ROI for cell" << std::endl;
	ImageType::Pointer cell_image = roi_grabber->GetROI(cell, roi_size);
	
	//WriteImage("cell.tif", cell_image);

	LoG *log_obj = new LoG();
	
	//Calculate the LoG on multiple scales and put them into a vector
	std::cout << "Calculating Multiscale LoG" << std::endl;
	std::vector<LoGImageType::Pointer> log_cellimage_vector = log_obj->RunMultiScaleLoG(cell, cell_image);

//	std::cout << "Starting ridge detection" << std::endl;
//	RidgeDetection(log_cellimage_vector, cell_image->GetBufferedRegion().GetSize(), critical_points_vector);
}

void MicrogliaRegionTracer::RidgeDetection( std::vector<LoGImageType::Pointer> log_cellimage_vector, ImageType::SizeType size, std::vector<ImageType::IndexType> &critical_points_vector )
{
	//Make a new image to store the critical points	
	ImageType::Pointer critical_point_image = ImageType::New();
	ImageType::IndexType start;
	start.Fill(0);

	ImageType::RegionType region(start, size);
	critical_point_image->SetRegions(region);
	critical_point_image->Allocate();
	critical_point_image->FillBuffer(0);

	//Make a new image to store the vesselness scores	
	VesselnessImageType::Pointer vesselness_image = VesselnessImageType::New();
	VesselnessImageType::IndexType vesselness_start;
	vesselness_start.Fill(0);

	VesselnessImageType::RegionType vesselness_region(vesselness_start, size);
	vesselness_image->SetRegions(vesselness_region);
	vesselness_image->Allocate();
	vesselness_image->FillBuffer(0);

	itk::ImageRegionIterator<ImageType> critical_point_img_iter(critical_point_image, critical_point_image->GetBufferedRegion());
	itk::ImageRegionIterator<VesselnessImageType> vesselness_image_iter(vesselness_image, vesselness_image->GetBufferedRegion());
	
	std::vector<LoGImageType::Pointer>::iterator log_cellimage_vector_iter;
	
	//For each log image, get the critical points
	for (log_cellimage_vector_iter = log_cellimage_vector.begin(); log_cellimage_vector_iter != log_cellimage_vector.end(); log_cellimage_vector_iter++)
	{
		std::cout << "Starting ridge detection for a new scale" << std::endl;
		LoGImageType::Pointer log_image = *log_cellimage_vector_iter;

		itk::Size<3> rad = {{1,1,1}};
		itk::NeighborhoodIterator<LoGImageType> neighbor_iter(rad, log_image, log_image->GetBufferedRegion());
		
		//Making a size value that removes some room (kind of like reverse padding) so we dont go out of bounds later
		itk::Size<3> unpad_size = log_image->GetBufferedRegion().GetSize();
		unpad_size[0] = unpad_size[0] - 3;
		unpad_size[1] = unpad_size[1] - 3; 
		unpad_size[2] = unpad_size[2] - 3;
		
		critical_point_img_iter.GoToBegin();
		neighbor_iter.GoToBegin();
		vesselness_image_iter.GoToBegin();

		while(!neighbor_iter.IsAtEnd()) 
		{
			itk::Index<3> index = neighbor_iter.GetIndex();
			
			//checking to see if we are at the edge of the image
			if ( (index[0] < 2) || (index[1] < 2) || (index[2] < 2) ||								
				(index[0] > (int)unpad_size[0]) || (index[1] > (int)unpad_size[1]) || (index[2] > (int)unpad_size[2]) )		
			{
				++critical_point_img_iter;
				++neighbor_iter;
				++vesselness_image_iter;
				continue;
			}

			//Diametrically opposing pairs
			////Summing the average of the greater of the diametrically pairs of pixels
			//float sum_of_greater_of_diametrically_opposing_pairs = 0.0;
			//for (unsigned int i=0; i < 13; ++i)
			//{
			//	sum_of_greater_of_diametrically_opposing_pairs += std::max<float>(neighbor_iter.GetPixel(i), neighbor_iter.GetPixel(26 - i));
			//}
			//
			////Now divide by the number of pairs to get the average
			//float avg_of_greater_of_diametrically_opposing_pairs = sum_of_greater_of_diametrically_opposing_pairs / 13.0f;

			//float center_pixel_intensity = neighbor_iter.GetPixel(13) ;	//center pixel intensity
			//if (center_pixel_intensity - avg_of_greater_of_diametrically_opposing_pairs > 0)
			//{
			//	critical_point_img_iter.Set(255);
			//	
			//	float vesselness_score = RunHessian(log_image, neighbor_iter);
			//	if (vesselness_image_iter.Get() < vesselness_score)
			//		vesselness_image_iter.Set(vesselness_score);
			//}

			//Local maximum
			bool local_maximum = true;
			for (int i = 0; i <= 26; ++i)
			{
				if (neighbor_iter.GetPixel(i) > neighbor_iter.GetPixel(13))
					local_maximum = false;
			}
			
			if (local_maximum)
			{
				critical_point_img_iter.Set(255);
		
				float vesselness_score = RunHessian(log_image, neighbor_iter);
				if (vesselness_image_iter.Get() < vesselness_score)
					vesselness_image_iter.Set(vesselness_score);
			}

			++critical_point_img_iter;
			++neighbor_iter;
			++vesselness_image_iter;
		}
	}

	//WriteImage("critical_points.TIF", critical_point_image);
	//WriteVesselnessImage("vesselness.mhd", vesselness_image);

	vesselness_image_iter.GoToBegin();

	while(!vesselness_image_iter.IsAtEnd())
	{
		if (vesselness_image_iter.Get() > 0.04)
			critical_points_vector.push_back(vesselness_image_iter.GetIndex());
		
		++vesselness_image_iter;
	}
}

double MicrogliaRegionTracer::RunHessian( LoGImageType::Pointer log_image, itk::NeighborhoodIterator<LoGImageType> neighbor_iter )
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

	double vesselness_score = ComputeVesselness(ev[0], ev[1], ev[2]);
	
	if (vesselness_score < 0)
		std::cout << "negative vesselness score detected" << std::endl;

	return vesselness_score;
}

double MicrogliaRegionTracer::ComputeVesselness( double ev1, double ev2, double ev3 )
{
	double ev1_magnitude = std::abs(ev1);
	double ev2_magnitude = std::abs(ev2);
	double ev3_magnitude = std::abs(ev3);
	
	if (ev1_magnitude < ev2_magnitude && ev1_magnitude < ev3_magnitude) //ev1 is the smallest eigenvalue
		return (ev2_magnitude + ev3_magnitude)/2 - ev1_magnitude;
	else if (ev2_magnitude < ev1_magnitude && ev2_magnitude < ev3_magnitude) //ev2 is the smallest eigenvalue
		return (ev1_magnitude + ev3_magnitude)/2 - ev2_magnitude;
	else
		return (ev1_magnitude + ev2_magnitude)/2 - ev3_magnitude;
}

void MicrogliaRegionTracer::BuildTree(Cell* cell, std::vector<ImageType::IndexType> &critical_points_vector)
{
	ImageType::IndexType cell_index;
	cell_index[0] = 100;
	cell_index[1] = 100;
	cell_index[2] = 50;

	critical_points_vector.push_back(cell_index);	//Add the centroid to the critical points

	double** AdjGraph = BuildAdjacencyGraph(critical_points_vector);

	Tree* tree = BuildMST(critical_points_vector, AdjGraph);

	for (int k = 0; k < critical_points_vector.size(); k++)
		delete[] AdjGraph[k];
	delete[] AdjGraph;

	delete tree;
}

double** MicrogliaRegionTracer::BuildAdjacencyGraph(std::vector<ImageType::IndexType> critical_points_vector)
{
	double** AdjGraph = new double*[critical_points_vector.size()];
	for (int k = 0; k < critical_points_vector.size(); k++)
		AdjGraph[k] = new double[critical_points_vector.size()];

	for (itk::uint64_t k = 0; k < critical_points_vector.size(); k++)
	{
		for (itk::uint64_t l = 0; l < critical_points_vector.size(); l++)
		{
			AdjGraph[k][l] = CalculateDistance(k, l, critical_points_vector);
		}
	}
	return AdjGraph;
}

double MicrogliaRegionTracer::CalculateDistance(itk::uint64_t k, itk::uint64_t l, std::vector<ImageType::IndexType> critical_points_vector)
{
	ImageType::IndexType node1 = critical_points_vector[k];
	ImageType::IndexType node2 = critical_points_vector[l];

	itk::int64_t distance_x = node1[0] - node2[0];
	itk::int64_t distance_y = node1[1] - node2[1];
	itk::int64_t distance_z = node1[2] - node2[2];

	return sqrt(pow(distance_x, 2.0) + pow(distance_y, 2.0) + pow(distance_z, 2.0));
}

Tree* MicrogliaRegionTracer::BuildMST(std::vector<MicrogliaRegionTracer::ImageType::IndexType> critical_points_vector, double** AdjGraph)
{	
	Tree* tree = new Tree();
	ImageType::IndexType root_index = critical_points_vector.back();
	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], critical_points_vector.size() - 1));
	
	std::ofstream traceFile;
	std::ostringstream traceFileNameStream;
	traceFileNameStream << root_index[0] << "_" << root_index[1] << "_" << root_index[2] << ".swc" << std::endl;
	traceFile.open(traceFileNameStream.str().c_str());
	//traceFile.open("trace.swc");
	std::cout << "Opening " << traceFileNameStream.str() << std::endl;

	//Root node should have infinite weights to connect to
	for (int m = 0; m < critical_points_vector.size(); m++)
		AdjGraph[m][critical_points_vector.size() - 1] = std::numeric_limits<double>::max();

	traceFile << "1 3 100 100 50 1 -1" << std::endl;
	//do this for each point but the root point
	for (itk::uint64_t l = 0; l < critical_points_vector.size() - 1; l++)
	{
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

			//Search through all the points and find the minimum distance
			for (itk::uint64_t k = 0; k < critical_points_vector.size(); k++)
			{
				if (AdjGraph[node->getID()][k] < minimum_node_distance && node->getID() != k)
				{
					minimum_parent_node = node;
					minimum_node_index_from_id = node->getID();
					minimum_node_index_to_id = k;
					minimum_node_distance = AdjGraph[node->getID()][k];
				}
			}	
		}

		//std::cout << "Found new edge from " << minimum_node_index_from_id << " to " << minimum_node_index_to_id << " Location: " << critical_points_vector[minimum_node_index_from_id][0] << " " << critical_points_vector[minimum_node_index_from_id][1] << " " << critical_points_vector[minimum_node_index_from_id][2] << " " << critical_points_vector[minimum_node_index_to_id][0] << " " << critical_points_vector[minimum_node_index_to_id][1] << " "  << critical_points_vector[minimum_node_index_to_id][2] << std::endl;
		
		if (minimum_node_index_from_id == critical_points_vector.size() - 1)	//parent is root
			traceFile << minimum_node_index_to_id + 2 << " 3 " << critical_points_vector[minimum_node_index_to_id][0] << " " << critical_points_vector[minimum_node_index_to_id][1] << " "  << critical_points_vector[minimum_node_index_to_id][2] << " 1 1" << std::endl;
		else
			traceFile << minimum_node_index_to_id + 2 << " 3 " << critical_points_vector[minimum_node_index_to_id][0] << " " << critical_points_vector[minimum_node_index_to_id][1] << " "  << critical_points_vector[minimum_node_index_to_id][2] << " 1 " << minimum_node_index_from_id + 2 << std::endl;

		//Set the weight of the AdjGraph entry for this minimum edge to infinite so it is not visited again
		for (int m = 0; m < critical_points_vector.size(); m++)
			AdjGraph[m][minimum_node_index_to_id] = std::numeric_limits<double>::max();
		AdjGraph[minimum_node_index_to_id][minimum_node_index_from_id] = std::numeric_limits<double>::max();

		//Connect the unconnected point
		Node* new_connected_node = new Node(critical_points_vector.at(minimum_node_index_to_id)[0], critical_points_vector.at(minimum_node_index_to_id)[1], critical_points_vector.at(minimum_node_index_to_id)[2], minimum_node_index_to_id);
		new_connected_node->SetParent(minimum_parent_node);
		tree->AddNode(new_connected_node, minimum_parent_node);
		//critical_points_vector.erase(critical_points_vector.begin()+minimum_node_index_to_id);
	}
	std::cout << "Closing " << traceFileNameStream.str() << std::endl;
	traceFile.close();
	return tree;
}
