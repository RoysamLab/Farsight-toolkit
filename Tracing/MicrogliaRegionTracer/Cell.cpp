#include "Cell.h"

#include "itkMaskNegatedImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"

#include "itkImageRegionIterator.h"

Cell::Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z)
{
	this->cell_x = cell_x;
	this->cell_y = cell_y;
	this->cell_z = cell_z;

	this->next_available_ID = 1;
}

itk::uint64_t Cell::getX() const
{
	return cell_x;
}

itk::uint64_t Cell::getY() const
{
	return cell_y;
}

itk::uint64_t Cell::getZ() const
{
	return cell_z;
}

void Cell::SetSize(const ImageType::SizeType & roi_size)
{
	this->roi_size = roi_size;
}

Cell::ImageType::SizeType Cell::GetSize() const
{
	return roi_size;
}

void Cell::SetOrigin(const ImageType::PointType & roi_origin)
{
	this->roi_origin = roi_origin;
}

Cell::ImageType::PointType Cell::GetOrigin() const
{
	return roi_origin;
}

void Cell::SetRequestedSize(const ImageType::SizeType & cell_requested_size)
{
	this->cell_requested_size = cell_requested_size;
}

Cell::ImageType::SizeType Cell::GetRequestedSize() const
{
	return cell_requested_size;
}

void Cell::SetShiftIndex(ImageType::IndexType shift_index)
{
	this->shift_index = shift_index;
}

Cell::ImageType::IndexType Cell::GetShiftIndex() const
{
	return shift_index;
}

void Cell::ComputeCriticalPointsVector(const ImageType::Pointer & critical_points_image)
{
	critical_points_queue.clear();
	
	itk::ImageRegionIterator<ImageType> critical_points_img_iter(critical_points_image, critical_points_image->GetLargestPossibleRegion());
	critical_points_img_iter.GoToBegin();
	while(!critical_points_img_iter.IsAtEnd())
	{
		if (critical_points_img_iter.Get() != 0)
			critical_points_queue.push_back(critical_points_img_iter.GetIndex());

		++critical_points_img_iter;
	}
}

void Cell::GetMask(const std::string & soma_filename)
{
	typedef itk::ImageFileReader< SomaImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(soma_filename);
	
	//PURPOSELY COMMENTED OUT SO YOU DO NOT TRY TO WRITE THIS CODE, THIS KILLS PERFORMANCE BECAUSE IT ATTEMPTS TO READ THE ENTIRE IMAGE FOR EACH CELL INSTEAD OF JUST THE ROI
	//try
	//{
	//	reader->Update();
	//}
	//catch (itk::ExceptionObject &err)
	//{
	//	std::cerr << "reader Exception: " << err << std::endl;
	//}

	typedef itk::RegionOfInterestImageFilter< SomaImageType, MaskImageType > ROIFilterType;
	ROIFilterType::Pointer roi_filter = ROIFilterType::New();
	
	ImageType::IndexType start;
	start[0] = this->roi_origin[0];
	start[1] = this->roi_origin[1];
	start[2] = this->roi_origin[2];

	ImageType::SizeType size = this->roi_size;
	
	ImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	roi_filter->SetRegionOfInterest(desiredRegion);
	roi_filter->SetInput(reader->GetOutput());

	try
	{
		//roi_filter->Update();
		ftk::TimeStampOverflowSafeUpdate(roi_filter.GetPointer());
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "roi_filter Exception: " << err << std::endl;
	}

	this->mask = roi_filter->GetOutput();
	this->mask->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...

	ImageType::PointType origin;   
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	this->mask->SetOrigin(origin);

	//Make the file name of the mask image
	std::stringstream mask_filename_stream;
	mask_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_mask.TIF";	//X_Y_Z_masked.TIF

	//Write the masked cell image
	//WriteImage(mask_filename_stream.str(), this->mask);

	//Get the label image from the binary image
	typedef itk::BinaryImageToLabelMapFilter< MaskImageType > BinaryToLabelFilterType;
	BinaryToLabelFilterType::Pointer labelMapFilter = BinaryToLabelFilterType::New();
	labelMapFilter->SetInput(this->mask);
	
	try
	{
		ftk::TimeStampOverflowSafeUpdate( labelMapFilter.GetPointer() );
		//labelMapFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelMapFilter exception: " << err << std::endl;
		std::cerr << "Mask image: " << this->mask << std::endl;
		std::cerr << this->mask << std::endl;
		std::cerr << labelMapFilter << std::endl;
	}

	BinaryToLabelFilterType::OutputImageType::Pointer label_map_image = labelMapFilter->GetOutput();
	label_map_image->DisconnectPipeline();

	typedef itk::LabelMapToLabelImageFilter< BinaryToLabelFilterType::OutputImageType, LabelImageType > LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelImageFilter = LabelMapToLabelImageFilterType::New();
	labelImageFilter->SetInput(label_map_image);
	try
	{
		ftk::TimeStampOverflowSafeUpdate( labelImageFilter.GetPointer() );
		//labelImageFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelImageFilter exception: " << err << std::endl;
	}

	this->soma_label_image = labelImageFilter->GetOutput();
	this->soma_label_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...
}

void Cell::ComputeMaskedImage()
{
	typedef itk::MaskNegatedImageFilter< ImageType, MaskImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetMaskImage(this->mask);
	maskFilter->SetInput(this->image);
	try
	{
		//ftk::TimeStampOverflowSafeUpdate( maskFilter.GetPointer() );
		maskFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "maskFilter Exception: " << err << std::endl;
	}

	this->masked_image = maskFilter->GetOutput();
	this->masked_image->DisconnectPipeline();		//Disconnect pipeline so we don't propagate...

	//Make the file name of the masked cell image
	std::stringstream masked_cell_filename_stream;
	masked_cell_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_masked.TIF";	//X_Y_Z_masked.TIF

	//Write the masked cell image
	//WriteImage(masked_cell_filename_stream.str(), this->masked_image);
}

void Cell::WriteImage(const std::string & filename, const itk::Image< unsigned char, 3>::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned char, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		//writer->Update();
		ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

void Cell::WriteImage(const std::string & filename, const itk::Image< unsigned short, 3>::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned short, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		//writer->Update();
		ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

void Cell::WriteImage(const std::string & filename, const itk::Image< float , 3 >::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< float, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		//writer->Update();
		ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

/*	This function will take a tree structure and write out the "local" and "global" SWC file corresponding to that tree structure.
 *	Local here means an SWC file that will fit into the region of interest we are interested in and global means that it will fit into the montage */
void Cell::WriteTreeToSWCFile(Tree * tree, std::string filename, std::string filename_local)
{
	std::cout << "Entering WriteTreeToSWCFile" << std::endl;
	std::ofstream traceFile, traceFile_local;
	
	std::cout << "Opening " << filename << std::endl;
	traceFile.open(filename.c_str());
	std::cout << "Opening " << filename_local << std::endl;
	traceFile_local.open(filename_local.c_str());
	
	Node* root = tree->GetRoot();
    
	itk::uint64_t tree_depth = 0; //root node is defined as tree depth 0
	WriteLinkToParent(root, tree_depth, traceFile, traceFile_local);	//Recursive function that does the actual tree traversal and writing the SWC lines
    
	traceFile.close();
	traceFile_local.close();
}

//This function does a depth-first traversal of the tree, maybe it makes sense to do a breadth-first traversal for some uses?
//CAREFUL: This function may overflow the stack due to recursion pushing parameters onto the stack and will crash if you have a tree deep enough. Rewrite iteratively or increase stack size if this is the case...
void Cell::WriteLinkToParent(Node* node, itk::uint64_t tree_depth, std::ofstream &traceFile, std::ofstream &traceFile_local)
{
	//Calculate some node indices
	ImageType::PointType node_index, node_index_local;
	node_index_local[0] = node->x;
	node_index_local[1] = node->y;
	node_index_local[2] = node->z;
	node_index[0] = node_index_local[0] + this->getX() - this->GetRequestedSize()[0]/2 - this->GetShiftIndex()[0];
	node_index[1] = node_index_local[1] + this->getY() - this->GetRequestedSize()[1]/2 - this->GetShiftIndex()[1];
	node_index[2] = node_index_local[2] + this->getZ() - this->GetRequestedSize()[2]/2 - this->GetShiftIndex()[2];
    
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
		WriteLinkToParent(child_node, tree_depth+1, traceFile, traceFile_local);
	}
    
	return;															//BASE CASE: Finished visiting the entire subtree that is rooted at this node
}
