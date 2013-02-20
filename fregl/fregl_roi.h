#ifndef fregl_roi_H
#define fregl_roi_H

#include "itkImage.h"
#include "fregl_util.h"
#include "fregl_image_manager.h"

template < class TPixel >
class fregl_roi
{
public:
	typedef TPixel							InputPixelType;			//TPixel determines our input pixel type
	typedef itk::Image< InputPixelType, 3 >	ImageType;			
	typedef typename ImageType::Pointer		ImageTypePointer;
	typedef typename ImageType::PointType	PointType;
	typedef typename ImageType::SizeType	SizeType;

	fregl_roi(	std::string joint_xforms_xml_file,	//Path to joint_transforms file, can take full path or relative path (from working directory) 
				std::string img_path,				//Path to the images (don't forget the ending "/")
				std::string anchor_image,			//The name of the anchor image (only the name is fine, since you specified the path in the last parameter)
				unsigned int cache_buffer_count,	//Number of images to keep in memory
				bool nearest_neighbor = false		//Whether or not to use nearest neighbor interpolation
				);									//TODO: Overload constructor to turn on/off file_caching
													
	
	void setROI(PointType roi_origin, SizeType roi_size);				//Set the ROI (remember this is in the coordinates of the montage)
	ImageTypePointer getROI();											//Returns you back the image  you asked for
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

