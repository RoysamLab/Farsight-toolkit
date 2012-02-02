#ifndef ROIGRABBER_H
#define ROIGRABBER_H

#include "Cell.h"
#include "fregl/fregl_roi.h"

class ROIGrabber
{
private:
	typedef fregl_roi::ImageType ImageType;

public:
	ROIGrabber(std::string joint_transforms_filename, std::string img_path, std::string anchor_filename);
	~ROIGrabber();
	
	ImageType::Pointer GetROI(Cell* cell, ImageType::SizeType roi_size);

private:
	fregl_roi *roi_filter;
	ImageType::Pointer image;	
};

#endif