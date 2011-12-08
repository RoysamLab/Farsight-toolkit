#include "fregl_roi.h"

fregl_roi::fregl_roi(std::string xml_file, std::string img_path, std::string anchor_image)
{
	this->xml_file = xml_file;
	this->img_path = img_path;
	this->anchor_image = anchor_image;
}

void fregl_roi::setROI(ImageType::IndexType roi_origin, ImageType::SizeType roi_size)
{
	this->roi_origin = roi_origin;
	this->roi_size = roi_size;
}

void fregl_roi::getMontage()
{
	//fregl_image_manager::Pointer region_montage = new fregl_image_manager(xml_file, img_path, anchor_image, false /* nearest neighbor */);
	/region_montage->set_regionofinterest(roi_origin, roi_size);
}
