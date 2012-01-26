#ifndef ROIGRABBER_H
#define ROIGRABBER_H

#include "Seed.h"
#include "fregl/fregl_roi.h"

class ROIGrabber
{
private:
	typedef fregl_roi::ImageType ImageType;

public:
	ROIGrabber();
	~ROIGrabber();
	
	ImageType::Pointer GetROI(Seed* seed, ImageType::SizeType roi_size);

private:
	fregl_roi *roi_filter;
	ImageType::Pointer image;	
};

#endif