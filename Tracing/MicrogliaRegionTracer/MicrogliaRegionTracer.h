#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"

#include <fstream>
#include <cstring>
#include <vector>

//Simple class to hold seed coordinates
#include "Seed.h"

//Fregl includes
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>
#include <fregl/fregl_image_manager.h>


class MicrogliaRegionTracer
{
public:
	typedef itk::Image<unsigned char, 3> ImageType;
	typedef itk::Image<float, 3> LoGImageType;

private:
	ImageType::Pointer image;
	std::vector<Seed*> seeds;

public:
	MicrogliaRegionTracer();
	~MicrogliaRegionTracer();

	void LoadImage(ImageType::Pointer image);
	void LoadImage(std::string filename);

	void LoadSeedPoints(std::string filename);

	void WriteImage(std::string filename, ImageType::Pointer image);
	void WriteLoGImage(std::string filename, LoGImageType::Pointer image);

	void RunLoG();

private:
	

	
};