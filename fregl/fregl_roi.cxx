#include "fregl_roi.h"

fregl_roi::fregl_roi(std::string joint_xforms_xml_file, std::string img_path, std::string anchor_image, bool nearest_neighbor)
{
	this->joint_xforms_xml_file = joint_xforms_xml_file;
	this->img_path = img_path;
	this->anchor_image = anchor_image;
	this->nearest_neighbor = nearest_neighbor;
}

void fregl_roi::setROI(ImageType::PointType roi_origin, ImageType::SizeType roi_size)
{
	this->roi_origin = roi_origin;
	this->roi_size = roi_size;
}

fregl_roi::ImageType::Pointer fregl_roi::getROI()
{
	fregl_image_manager::Pointer region_montage = new fregl_image_manager(joint_xforms_xml_file, img_path, anchor_image, nearest_neighbor);
	region_montage->set_regionofinterest(roi_origin, roi_size);
	
	region_montage->Update();
	return region_montage->GetOutput();
}
