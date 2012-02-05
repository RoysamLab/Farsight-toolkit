#ifndef fregl_roi_H
#define fregl_roi_H

#include "itkImage.h"
#include "fregl_image_manager.h"

class fregl_roi
{
public:
	typedef itk::Image< unsigned char, 3 >  ImageType;

	fregl_roi(std::string joint_xforms_xml_file, std::string img_path, std::string anchor_image, bool nearest_neighbor = false);
	~fregl_roi();

	void setROI(ImageType::PointType roi_origin, ImageType::SizeType roi_size);
	ImageType::Pointer getROI();
    fregl_image_manager::Pointer getImageManager ();
	
private:
	ImageType::PointType roi_origin;
	ImageType::SizeType roi_size;

	std::string joint_xforms_xml_file;
	std::string img_path;
	std::string anchor_image;
	bool nearest_neighbor;
    fregl_image_manager::Pointer region_montage;
};

#endif

