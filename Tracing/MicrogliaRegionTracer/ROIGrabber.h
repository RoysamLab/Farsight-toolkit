#ifndef ROIGRABBER_H
#define ROIGRABBER_H

#include "fregl/fregl_roi.h"
#include "Cell.h"

class ROIGrabber
{
public:
	typedef Cell::InputPixelType InputPixelType;
	typedef Cell::ImageType	ImageType;
	

public:
	explicit ROIGrabber(const std::string & joint_transforms_filename, const std::string & img_path, const std::string & anchor_filename);
	~ROIGrabber();
	
	ImageType::Pointer GetROI(const Cell & cell, const ImageType::SizeType & roi_size, ImageType::IndexType & shift_index);

private:
	fregl_roi< InputPixelType > roi_filter;
	ImageType::Pointer image;
};

#endif