#ifndef MICROGLIAREGIONTRACER_H
#define MICROGLIAREGIONTRACER_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


#include <fstream>
#include <cstring>
#include <vector>


#include "Seed.h"				//Simple class to hold seed coordinates
#include "ROIGrabber.h"
#include "LoG.h"

#include "time.h"

class MicrogliaRegionTracer
{
public:
	typedef fregl_roi::ImageType ImageType;
	typedef LoG::LoGImageType LoGImageType;

private:
	ImageType::Pointer image;
	std::vector<Seed*> seeds;
	ROIGrabber* roi_grabber;

public:
	MicrogliaRegionTracer();
	~MicrogliaRegionTracer();

	void LoadImage(ImageType::Pointer image);
	void LoadImage(std::string filename);

	void LoadSeedPoints(std::string filename);

	void WriteImage(std::string filename, ImageType::Pointer image);	
	void WriteInitialMicrogliaImages();
	
	void Trace();

	void CalculateCandidatePixel(Seed* seed);
	void RidgeDetection(std::vector<LoGImageType::Pointer> log_seedimage_vector);


private:

	
};

#endif