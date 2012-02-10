#include "ROIGrabber.h"

ROIGrabber::ROIGrabber(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename)
{
	roi_filter = new fregl_roi< InputPixelType >(joint_transforms_filename, img_path, anchor_filename, true);
}

ROIGrabber::ImageType::Pointer ROIGrabber::GetROI(Cell* cell, ImageType::SizeType roi_size, ImageType::IndexType &shift_index)
{
	ImageType::PointType roi_origin;

	roi_origin[0] = (itk::int64_t)cell->getX() - (itk::int64_t)roi_size[0]/2;
	roi_origin[1] = (itk::int64_t)cell->getY() - (itk::int64_t)roi_size[1]/2;
	roi_origin[2] = (itk::int64_t)cell->getZ() - (itk::int64_t)roi_size[2]/2;

	if ((itk::int64_t)cell->getX() - (itk::int64_t)roi_size[0]/2 < 0)	
		shift_index[0] = (itk::int64_t)cell->getX() - (itk::int64_t)roi_size[0]/2; 
	else 
		shift_index[0] = 0;
	
	if ((itk::int64_t)cell->getY() - (itk::int64_t)roi_size[1]/2 < 0) 
		shift_index[1] = (itk::int64_t)cell->getY() - (itk::int64_t)roi_size[1]/2; 
	else 
		shift_index[1] = 0;
	
	if ((itk::int64_t)cell->getZ() - (itk::int64_t)roi_size[2]/2 < 0) 
		shift_index[2] = (itk::int64_t)cell->getZ() - (itk::int64_t)roi_size[2]/2; 
	else 
		shift_index[2] = 0;

	roi_filter->setROI(roi_origin, roi_size);

	return roi_filter->getROI();
}