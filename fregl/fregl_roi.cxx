#include "fregl_roi.h"

template < class TPixel >
fregl_roi< TPixel >::fregl_roi(std::string joint_xforms_xml_file, std::string img_path, std::string anchor_image, unsigned int cache_buffer_count, bool nearest_neighbor ) {
    this->joint_xforms_xml_file = joint_xforms_xml_file;
    this->img_path = img_path;
    this->anchor_image = anchor_image;
    this->nearest_neighbor = nearest_neighbor;
    this->region_montage = new fregl_image_manager< TPixel >(joint_xforms_xml_file, img_path, anchor_image, nearest_neighbor);
	
	this->region_montage->set_cache_buffer_count(cache_buffer_count);
	this->region_montage->set_file_caching(true);
}

template < class TPixel >
void fregl_roi< TPixel >::setROI(PointType roi_origin, SizeType roi_size) {
    this->roi_origin = roi_origin;
    this->roi_size = roi_size;
	region_montage->set_regionofinterest(roi_origin, roi_size);
}

template < class TPixel >
typename fregl_roi< TPixel >::ImageTypePointer 
fregl_roi< TPixel >::getROI() {
    region_montage->Update();
    return region_montage->GetOutput();
}

template < class TPixel >
typename fregl_image_manager< TPixel >::Pointer 
fregl_roi< TPixel >::getImageManager() {
    return this->region_montage;
}

//Explicit Instantiation
template class fregl_roi< unsigned char >;
template class fregl_roi< unsigned short >;