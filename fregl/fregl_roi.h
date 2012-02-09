#ifndef fregl_roi_H
#define fregl_roi_H

#include "itkImage.h"
#include "fregl_util.h"
#include "fregl_image_manager.h"

template < class TPixel >
class fregl_roi
{
public:
	typedef TPixel InputPixelType;
	typedef itk::Image< InputPixelType, 3 >	ImageType;
	typedef typename ImageType::Pointer		ImageTypePointer;
	typedef typename ImageType::PointType	PointType;
	typedef typename ImageType::SizeType	SizeType;

	fregl_roi(std::string joint_xforms_xml_file, std::string img_path, std::string anchor_image, bool nearest_neighbor = false);
	
	void setROI(PointType roi_origin, SizeType roi_size);
	ImageTypePointer getROI();
    typename fregl_image_manager< TPixel >::Pointer getImageManager();
	
private:
	PointType roi_origin;
	SizeType roi_size;

	std::string joint_xforms_xml_file;
	std::string img_path;
	std::string anchor_image;
	bool nearest_neighbor;
    typename fregl_image_manager< TPixel >::Pointer region_montage;
};

#endif

