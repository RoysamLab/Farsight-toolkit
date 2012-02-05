#ifndef ROIGRABBER_H
#define ROIGRABBER_H

#include "fregl/fregl_roi.h"
#include "Cell.h"

class ROIGrabber
{
public:
	typedef fregl_roi::ImageType ImageType;

public:
	ROIGrabber(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename);
	~ROIGrabber();
	
	ImageType::Pointer GetROI(Cell* cell, ImageType::SizeType roi_size, ImageType::IndexType &shift_index);

private:
	fregl_roi *roi_filter;
	ImageType::Pointer image;
};

#endif