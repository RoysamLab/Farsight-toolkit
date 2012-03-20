#include "Cell.h"

Cell::Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z)
{
	this->cell_x = cell_x;
	this->cell_y = cell_y;
	this->cell_z = cell_z;
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

void Cell::SetSize(ImageType::SizeType roi_size)
{
	this->roi_size = roi_size;
}

Cell::ImageType::SizeType Cell::GetSize()
{
	return roi_size;
}

void Cell::SetOrigin(ImageType::PointType roi_origin)
{
	this->roi_origin = roi_origin;
}

Cell::ImageType::PointType Cell::GetOrigin()
{
	return roi_origin;
}

void Cell::setRequestedSize(ImageType::SizeType cell_requested_size)
{
	this->cell_requested_size = cell_requested_size;
}

Cell::ImageType::SizeType Cell::getRequestedSize()
{
	return cell_requested_size;
}

void Cell::setShiftIndex(ImageType::IndexType shift_index)
{
	this->shift_index = shift_index;
}

Cell::ImageType::IndexType Cell::getShiftIndex()
{
	return shift_index;
}

void Cell::ComputeCriticalPointsVector(ImageType::Pointer critical_points_image)
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

void Cell::GetMask(std::string soma_filename)
{
	typedef itk::ImageFileReader< MaskImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(soma_filename);
	
	//PURPOSELY COMMENTED OUT SO YOU DO NOT TRY TO WRITE THIS CODE, THIS KILLS PERFORMANCE BECAUSE IT ATTEMPTS TO READ THE ENTIRE IMAGE FOR EACH CELL INSTEAD OF THE ROI
	//try
	//{
	//	reader->Update();
	//}
	//catch (itk::ExceptionObject &err)
	//{
	//	std::cerr << "reader Exception: " << err << std::endl;
	//}

	typedef itk::RegionOfInterestImageFilter< SomaImageType, ImageType > ROIFilterType;
	ROIFilterType::Pointer roi_filter = ROIFilterType::New();
	
	ImageType::IndexType start;
	start[0] = roi_origin[0];
	start[1] = roi_origin[1];
	start[2] = roi_origin[2];

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
		ftk::TimeStampOverflowSafeUpdate( maskFilter.GetPointer() );
		//maskFilter->Update();
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

void Cell::WriteImage(std::string filename, itk::Image< unsigned char, 3>::Pointer image)
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

void Cell::WriteImage(std::string filename, itk::Image< unsigned short, 3>::Pointer image)
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

void Cell::WriteImage(std::string filename, itk::Image< float , 3 >::Pointer image)
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