#include "Cell.h"

Cell::Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z)
{
	this->cell_x = cell_x;
	this->cell_y = cell_y;
	this->cell_z = cell_z;

	this->vesselness_image_maximum_intensity = std::numeric_limits<float>::min();
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
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}

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
		roi_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "roi_filter Exception: " << err << std::endl;
	}

	this->mask = roi_filter->GetOutput();

	ImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	this->mask->SetOrigin(origin);

	//Make the file name of the mask image
	std::stringstream mask_filename_stream;
	mask_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_mask.TIF";	//X_Y_Z_masked.TIF

	std::cout << this->mask->GetLargestPossibleRegion().GetSize() << " " << this->mask->GetOrigin() << std::endl;

	//Write the masked cell image
	WriteImage(mask_filename_stream.str(), this->mask);

	//Get the label image from the binary image
	typedef itk::BinaryImageToLabelMapFilter<MaskImageType> BinaryToLabelFilterType;
	BinaryToLabelFilterType::Pointer labelMapFilter = BinaryToLabelFilterType::New();
	labelMapFilter->SetInput(this->mask);
	
	try
	{
		labelMapFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelMapFilter exception: " << err << std::endl;
	}

	typedef itk::LabelMapToLabelImageFilter< BinaryToLabelFilterType::OutputImageType, LabelImageType > LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelImageFilter = LabelMapToLabelImageFilterType::New();
	labelImageFilter->SetInput(labelMapFilter->GetOutput());
	try
	{
		labelImageFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelImageFilter exception: " << err << std::endl;
	}

	this->soma_label_image = labelImageFilter->GetOutput();
}

void Cell::ComputeMaskedImage()   
{
	typedef itk::MaskNegatedImageFilter< ImageType, MaskImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetMaskImage(this->mask);
	maskFilter->SetInput(this->image);
	try
	{
		maskFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "maskFilter Exception: " << err << std::endl;
	}

	this->masked_image = maskFilter->GetOutput();

	//Make the file name of the masked cell image
	std::stringstream masked_cell_filename_stream;
	masked_cell_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_masked.TIF";	//X_Y_Z_masked.TIF

	//Write the masked cell image
	WriteImage(masked_cell_filename_stream.str(), this->masked_image);
}

void Cell::WriteImage(std::string filename, itk::Image< unsigned char, 3>::Pointer image)
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

void Cell::WriteImage(std::string filename, itk::Image< unsigned short, 3>::Pointer image)
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

void Cell::WriteImage(std::string filename, itk::Image< float , 3 >::Pointer image)
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

void Cell::MaximumEntropyThreshold()
{
	typedef itk::MaximumEntropyThresholdImageFilter< ImageType, ImageType > MaximumEntropyThresholdImageFilterType;
	MaximumEntropyThresholdImageFilterType::Pointer MaximumEntropy_threshold_filter = MaximumEntropyThresholdImageFilterType::New();
	MaximumEntropy_threshold_filter->SetInput(this->image);
	MaximumEntropy_threshold_filter->Update();

	//Inverts the image
	typedef itk::ShiftScaleImageFilter< ImageType, ImageType > InvertFilterType;
	InvertFilterType::Pointer invertFilter = InvertFilterType::New();
	invertFilter->SetShift(-1.0 * std::numeric_limits<ImageType::PixelType>::max());
	invertFilter->SetScale(-1.0 );
	invertFilter->SetInput(MaximumEntropy_threshold_filter->GetOutput());

	try
	{
		invertFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "invertFilter Exception: " << err << std::endl;
	}

	this->thresholded_image = invertFilter->GetOutput();

	std::ostringstream MaximumEntropy_filename_stream;
	MaximumEntropy_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_MaximumEntropy.mhd";
	WriteImage(MaximumEntropy_filename_stream.str(), this->thresholded_image);

	
}

void Cell::Skeletonize()
{
	//typedef itk::SignedMaurerDistanceMapImageFilter< ImageType, DistanceImageType > DistanceMapFilterType;
	//DistanceMapFilterType::Pointer dist_map_filter = DistanceMapFilterType::New();
	//dist_map_filter->SetInput(this->thresholded_image);
	//dist_map_filter->SetSquaredDistance(false);
	//dist_map_filter->SetInsideIsPositive(true);

	//try
	//{
	//	dist_map_filter->Update();
	//}
	//catch (itk::ExceptionObject &err)
	//{                                             
	//	std::cerr << "dist_map_filter exception at " << __FILE__ << " " << __LINE__ << std::endl;
	//	std::cerr << &err << std::endl;
	//}

	//std::ostringstream dist_image_filename_stream;
	//dist_image_filename_stream << cell->getX() << "_" << cell->getY() << "_" << cell->getZ() << "_dist.mhd";
	//WriteImage(dist_image_filename_stream.str(), dist_map_filter->GetOutput());
	//this->distance_map_image = dist_map_filter->GetOutput();

	//Generate the skeleton image
	//skeleton_image = ImageType::New();
	//ImageType::IndexType start;
	//start.Fill(0);
	//ImageType::SizeType size = cell->thresholded_image->GetLargestPossibleRegion().GetSize();

	//ImageType::RegionType region(start, size);
	//skeleton_image->SetRegions(region);
	//skeleton_image->Allocate();
	//skeleton_image->FillBuffer(0);

	//itk::Size<3> rad = {{1,1,1}};
	//itk::NeighborhoodIterator<DistanceImageType> dist_neighbor_iter(rad, this->distance_map_image, this->distance_map_image->GetLargestPossibleRegion());
	//itk::ImageRegionIterator<ImageType> skeleton_img_iter(this->skeleton_image, this->skeleton_image->GetLargestPossibleRegion());

	//skeleton_img_iter.GoToBegin();
	//dist_neighbor_iter.GoToBegin();

	//while (!dist_neighbor_iter.IsAtEnd())
	//{
	//	DistanceImageType::PixelType center_pixel_intensity = dist_neighbor_iter.GetCenterPixel();
	//	unsigned int num_greater = 0;

	//	for (int k = 0; k < 26; k++)
	//	{
	//		if (k != 13)
	//		{
	//			bool IsInBounds;
	//			if (dist_neighbor_iter.GetPixel(k, IsInBounds) > center_pixel_intensity && IsInBounds)
	//				++num_greater;
	//		}
	//	}

	//	if (num_greater < 2) //Largest value or the 2nd largest value in the neighborhood
	//		skeleton_img_iter.Set(std::numeric_limits< ImageType::PixelType >::max());
	//	
	//	++skeleton_img_iter;
	//	++dist_neighbor_iter;
	//}

	typedef itk::BinaryFillholeImageFilter< ImageType > FillholeFilterType;
	FillholeFilterType::Pointer fillhole_filter = FillholeFilterType::New();
	fillhole_filter->SetInput(this->thresholded_image);
	fillhole_filter->SetFullyConnected(false);
	fillhole_filter->SetForegroundValue(255);

	typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType > SkeletonFilterType;
	SkeletonFilterType::Pointer skeleton_filter = SkeletonFilterType::New();
	skeleton_filter->SetInput(fillhole_filter->GetOutput());
	skeleton_filter->Update();

	this->skeleton_image = skeleton_filter->GetOutput();

	std::ostringstream skeleton_filename_stream;
	skeleton_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_skeleton.mhd";
	WriteImage(skeleton_filename_stream.str(), this->skeleton_image);
}