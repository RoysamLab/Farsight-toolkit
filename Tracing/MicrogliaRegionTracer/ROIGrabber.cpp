#include "ROIGrabber.h"

//Constructor, important that you dont make a new fregl_roi for each new object you want to get, otherwise you will lose the in-memory cache
ROIGrabber::ROIGrabber(const std::string & joint_transforms_filename, const std::string & img_path, const std::string & anchor_filename)
	: roi_filter(fregl_roi< InputPixelType >(joint_transforms_filename, img_path, anchor_filename, true, 6))
{
}

ROIGrabber::~ROIGrabber()
{
}

ROIGrabber::ImageType::Pointer ROIGrabber::GetROI(const Cell & cell, const ImageType::SizeType & roi_size, ImageType::IndexType & shift_index)
{
	ImageType::PointType roi_origin;

	roi_origin[0] = (itk::int64_t)cell.getX() - (itk::int64_t)roi_size[0]/2;	//X,Y,Z coordinates are middle of the cell, so we have to subtract half the size to get the origin of the ROI
	roi_origin[1] = (itk::int64_t)cell.getY() - (itk::int64_t)roi_size[1]/2;
	roi_origin[2] = (itk::int64_t)cell.getZ() - (itk::int64_t)roi_size[2]/2;

	//Calculate how much of the left, top, and the top-of-stack got cut off by our size
	if (roi_origin[0] < 0)	
		shift_index[0] = roi_origin[0]; 
	else
		shift_index[0] = 0;
	
	if (roi_origin[1] < 0) 
		shift_index[1] = roi_origin[1]; 
	else 
		shift_index[1] = 0;
	
	if (roi_origin[2] < 0) 
		shift_index[2] = roi_origin[2]; 
	else 
		shift_index[2] = 0;
	
	this->roi_filter.setROI(roi_origin, roi_size);	//Set the ROI

	return this->roi_filter.getROI();				//Get the ROI and return as itkImage
}