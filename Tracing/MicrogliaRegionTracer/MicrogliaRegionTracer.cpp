#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename, std::string soma_filename)
{
	this->roi_grabber = new ROIGrabber(joint_transforms_filename, img_path, anchor_filename);
	this->soma_filename = soma_filename;

	itk::MultiThreader::SetGlobalDefaultNumberOfThreads( 16 );	//Acquiring threads through OpenMP, so no need for ITK threads
}

void MicrogliaRegionTracer::LoadCellPoints(std::string seedpoints_filename)
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

		//Need to do this to make the MaskImageFilter work since it expects similar origins
		ImageType::PointType origin = temp_cell_image->GetOrigin();
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

		//Get the mask
		//cell->GetMask(this->soma_filename);
		
		std::cout << "Calculating candidate pixels for a new cell" << std::endl;
		CalculateCandidatePixels(cell);

		std::cout << "Detected " << cell->critical_points_queue.size() << " critical points" << std::endl;
		
		std::cout << "Tree Building" << std::endl;
		BuildTree(cell);

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
	//Add the seed to the critical points
	itk::Index<3> seed_index = {{	cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0], 
									cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1], 
									cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2] }};

	cell->critical_points_queue.push_front(seed_index);
	
	//Calculate the LoG on multiple scales and store into an image
	LoG *log_obj = new LoG();
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

	//Make a iterator for the image and a neighborhood around the current point we are visiting
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator< ImageType > critical_point_img_iter(cell->critical_point_image, cell->critical_point_image->GetLargestPossibleRegion());
	itk::ConstNeighborhoodIterator< LoGImageType > LoG_neighbor_iter(rad, cell->multiscale_LoG_image, cell->multiscale_LoG_image->GetLargestPossibleRegion());

	while(!LoG_neighbor_iter.IsAtEnd()) 
	{
		unsigned int neighborhood_size = (rad[0] * 2 + 1) * (rad[1] * 2 + 1) * (rad[2] * 2 + 1);
		unsigned int center_pixel_offset_index = neighborhood_size / 2;

		LoGImageType::PixelType center_pixel_intensity = LoG_neighbor_iter.GetPixel(center_pixel_offset_index);
		
		if (center_pixel_intensity >= 0.03)	//Must have greater than 3% response from LoG to even be considered for local maximum (so we dont choose points that are part of noise)
		{	
			bool local_maximum = true;
			
			for (unsigned int neighborhood_index = 0; neighborhood_index < neighborhood_size; ++neighborhood_index)
			{
				if (neighborhood_index != center_pixel_offset_index)
				{
					bool isInBounds; //true if the pixel is in bounds
					LoGImageType::PixelType neighbor_pixel_log_intensity = LoG_neighbor_iter.GetPixel(neighborhood_index, isInBounds);

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
		++LoG_neighbor_iter;
	}

	std::ostringstream criticalpointsFileNameStream;
	criticalpointsFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_critical.mhd";
	//cell->WriteImage(criticalpointsFileNameStream.str(), cell->critical_point_image);
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
		ftk::TimeStampOverflowSafeUpdate( multiscale_hessian_filter.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "multiscale_hessian_filter exception: " << &err << std::endl;
	}

	cell->vesselness_image = multiscale_hessian_filter->GetOutput();
	cell->vesselness_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...
	
	std::ostringstream vesselnessFileNameStream;

	vesselnessFileNameStream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_vesselness.mhd";
	//cell->(vesselnessFileNameLocalStream.str(), cell->vesselness_image);
}

void MicrogliaRegionTracer::BuildTree(Cell* cell)
{
	//Seed point should be added in Ridge Detection
	//ImageType::IndexType cell_index;
	//cell_index[0] = cell->getRequestedSize()[0]/2 + cell->getShiftIndex()[0];
	//cell_index[1] = cell->getRequestedSize()[1]/2 + cell->getShiftIndex()[1];
	//cell_index[2] = cell->getRequestedSize()[2]/2 + cell->getShiftIndex()[2];

	//cell->critical_points_queue.push_front(cell_index);	//Add the centroid to the critical points

	std::cout << "Building Adjacency Graph" << std::endl;
	double** AdjGraph = BuildAdjacencyGraph(cell);

	/*Tree* tree = BuildMST1(cell, AdjGraph);

	std::ostringstream swc_filename_stream, swc_filename_stream_local;
	swc_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree.swc";
	swc_filename_stream_local << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_tree_local.swc";


	WriteTreeToSWCFile(tree, cell, swc_filename_stream.str(), swc_filename_stream_local.str());*/

	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		delete[] AdjGraph[k];
	delete[] AdjGraph;

	//delete tree;
}

double** MicrogliaRegionTracer::BuildAdjacencyGraph(Cell* cell)
{
	double** AdjGraph = new double*[cell->critical_points_queue.size()];
	for (int k = 0; k < cell->critical_points_queue.size(); k++)
		AdjGraph[k] = new double[cell->critical_points_queue.size()];

	//#pragma omp parallel for
	for (itk::int64_t k = 0; k < cell->critical_points_queue.size(); k++)
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

//This is the cost function to connect 2 points
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

	////If we are connecting from the root node
	//if (node_from == 0)
	//{	
	//	if (cell->soma_label_image->GetPixel(node1) == cell->soma_label_image->GetPixel(node2)) //If it is inside the soma our centroid belongs to, we give it 0 weight
	//		return 0;
	//	else if (cell->soma_label_image->GetPixel(node2) != 0)											//We are trying to connect to a pixel in another soma, so we give it infinite weight
	//		return std::numeric_limits< double >::max();	
	//}
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
		return std::numeric_limits< double >::max();	//Thresholding based on jumping large gasp
	else
		return mag_trace_vector;
}

//Tree* MicrogliaRegionTracer::BuildMST1(Cell* cell, double** AdjGraph)
//{	
//	Tree* tree = new Tree();
//	ImageType::IndexType root_index = cell->critical_points_queue.front();
//	tree->SetRoot(new Node(root_index[0], root_index[1], root_index[2], 1));
//
//	//Root node should have infinite weights to connect to
//	for (int m = 0; m < cell->critical_points_queue.size(); m++)
//		AdjGraph[m][0] = std::numeric_limits<double>::max();
//
//	//Prim's algorithm
//	//for each node but the last, find its child
//	for (itk::uint64_t l = 0; l < cell->critical_points_queue.size() - 1; l++)
//	{
//		//std::cout << "Calculating nearest neighbor for point " << l << std::endl;
//		itk::uint64_t minimum_node_index_from_id = 0;
//		itk::uint64_t minimum_node_index_to_id = 0;
//		double minimum_node_distance = std::numeric_limits<double>::max();
//		Node* minimum_parent_node = NULL;
//		
//		//For each connected point
//		std::vector<Node*> member_nodes = tree->GetMemberNodes();
//		std::vector<Node*>::iterator member_nodes_iter;
//		for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)
//		{
//			Node* node = *member_nodes_iter;
//
//			//Search through all the nodes and find the minimum distance
//			for (itk::uint64_t k = 0; k < cell->critical_points_queue.size(); k++)
//			{
//				itk::uint64_t node_index_from_id = node->getID() - 1; //Node IDs are 1 greater than than vector indices, so we subtract 1 here
//				
//				if (AdjGraph[node_index_from_id][k] < minimum_node_distance) //from l (connected point) to k (unconnected point) if the current distance is less than the minimum distance
//				{
//					minimum_parent_node = node;
//					minimum_node_index_from_id = node_index_from_id;			
//					minimum_node_index_to_id = k;
//					minimum_node_distance = AdjGraph[minimum_node_index_from_id][minimum_node_index_to_id];
//				}
//			}	
//		}
//
//		if (minimum_node_distance >= std::numeric_limits< double >::max() / 2)
//			break;	//Minimum distance way too far
//		
//		//std::cout << "Found new edge from " << minimum_node_index_from_id << " to " << minimum_node_index_to_id << " Location: " << cell->critical_points_queue[minimum_node_index_from_id][0] << " " << cell->critical_points_queue[minimum_node_index_from_id][1] << " " << cell->critical_points_queue[minimum_node_index_from_id][2] << " " << cell->critical_points_queue[minimum_node_index_to_id][0] << " " << cell->critical_points_queue[minimum_node_index_to_id][1] << " "  << cell->critical_points_queue[minimum_node_index_to_id][2] << std::endl;
//		ImageType::IndexType cell_origin;
//		cell_origin[0] = cell->critical_points_queue[minimum_node_index_to_id][0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
//		cell_origin[1] = cell->critical_points_queue[minimum_node_index_to_id][1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
//		cell_origin[2] = cell->critical_points_queue[minimum_node_index_to_id][2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];
//		
//		ImageType::IndexType cell_origin_local;
//		cell_origin_local[0] = cell->critical_points_queue[minimum_node_index_to_id][0];
//		cell_origin_local[1] = cell->critical_points_queue[minimum_node_index_to_id][1];
//		cell_origin_local[2] = cell->critical_points_queue[minimum_node_index_to_id][2];
//
//		
//		for (int m = 0; m < cell->critical_points_queue.size(); m++)
//			AdjGraph[m][minimum_node_index_to_id] = std::numeric_limits<double>::max();						 //Node already has parents, we dont want it to have more parents
//		AdjGraph[minimum_node_index_to_id][minimum_node_index_from_id] = std::numeric_limits<double>::max(); //So the parent doesn't become the child of its child
//
//		//Connect the unconnected point
//		Node* new_connected_node = new Node(cell_origin_local[0], cell_origin_local[1], cell_origin_local[2], minimum_node_index_to_id + 1); //vector indices are 1 less than SWC IDs, so we add one here
//		new_connected_node->SetParent(minimum_parent_node);
//		tree->AddNode(new_connected_node, minimum_parent_node);
//	}
//
//	return tree;
//}

//void MicrogliaRegionTracer::WriteTreeToSWCFile(Tree* tree, Cell* cell, std::string filename, std::string filename_local)
//{
//	std::cout << "Entering WriteTreeToSWCFile" << std::endl;
//	std::ofstream traceFile, traceFile_local;
//	
//	std::cout << "Opening " << filename << std::endl;
//	traceFile.open(filename.c_str());
//	std::cout << "Opening " << filename_local << std::endl;
//	traceFile_local.open(filename_local.c_str());
//	
//	Node* root = tree->getRoot();
//
//	itk::uint64_t tree_depth = 0; //root node is defined as tree depth 0
//	WriteLinkToParent(root, tree_depth, cell, traceFile, traceFile_local);	//Recursive function that does the actual tree traversal and writing the SWC lines
//
//	traceFile.close();
//	traceFile_local.close();	
//}

////This function does a depth-first traversal of the tree, maybe it makes sense to do a breadth-first traversal for some uses?
////CAREFUL: This function may overflow the stack due to recursion pushing parameters onto the stack and will crash if you have a tree deep enough. Rewrite iteratively or increase stack size if this is the case...
//void MicrogliaRegionTracer::WriteLinkToParent(Node* node, itk::uint64_t tree_depth, Cell* cell, std::ofstream &traceFile, std::ofstream &traceFile_local)
//{
//	//Calculate some node indices
//	ImageType::IndexType node_index, node_index_local;
//	node_index_local[0] = node->x;
//	node_index_local[1] = node->y;
//	node_index_local[2] = node->z;
//	node_index[0] = node_index_local[0] + cell->getX() - cell->getRequestedSize()[0]/2 - cell->getShiftIndex()[0];
//	node_index[1] = node_index_local[1] + cell->getY() - cell->getRequestedSize()[1]/2 - cell->getShiftIndex()[1];
//	node_index[2] = node_index_local[2] + cell->getZ() - cell->getRequestedSize()[2]/2 - cell->getShiftIndex()[2];
//
//	itk::int64_t parent_node_id;
//	if (node->GetParent() == NULL) 
//		parent_node_id = -1; //This node has no parent, so this is the root node, so parent ID is -1
//	else
//		parent_node_id = node->GetParent()->getID();
//
//	std::vector< Node* > children = node->GetChildren();
//	if (tree_depth == 1 && children.size() == 0)
//		return;														//BASE CASE: Don't write out a trace if we are at depth one and have no children because we are a trace to the edge of the soma
//
//	//Write out the SWC lines
//	traceFile		<< node->getID() << " 3 " << node_index[0]			<< " " << node_index[1]			<< " "  << node_index[2]		<< " 1 " << parent_node_id << std::endl;
//	traceFile_local << node->getID() << " 3 " << node_index_local[0]	<< " " << node_index_local[1]	<< " "  << node_index_local[2]	<< " 1 " << parent_node_id << std::endl;
//	
//	if (children.size() == 0) 
//		return;														//BASE CASE: No children, so don't visit them
//
//	std::vector< Node* >::iterator children_iter;
//	for (children_iter = children.begin(); children_iter != children.end(); ++children_iter)
//	{
//		Node* child_node = *children_iter;
//		WriteLinkToParent(child_node, tree_depth+1, cell, traceFile, traceFile_local);
//	}
//
//	return;															//BASE CASE: Finished visiting the entire subtree that is rooted at this node
//}
