#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer()
{
	roi_grabber = new ROIGrabber("E:/Farsight_Images/MicrogliaRegionTracer/DAPI/joint_transforms.xml", "E:/Farsight_Images/MicrogliaRegionTracer/DAPI/", "8bitkt06045_w410DAPIdsu.TIF");
}

void MicrogliaRegionTracer::LoadImage(ImageType::Pointer image)
{
	this->image = image;
}

void MicrogliaRegionTracer::LoadImage(std::string filename)
{
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "reader Exception: " << err << std::endl;
	}

	LoadImage(reader->GetOutput());
}

void MicrogliaRegionTracer::LoadSeedPoints(std::string filename)
{
	std::ifstream seed_point_file;
	seed_point_file.open(filename.c_str());

	while (!seed_point_file.eof())
	{
		size_t seedX, seedY, seedZ;
		seed_point_file >> seedX >> seedY >> seedZ;

		//std::cout << "Reading in seed: (" << seedX << ", " << seedY << ", " << seedZ << ")" << std::endl;
		seeds.push_back(new Seed(seedX, seedY, seedZ));
	}
}

void MicrogliaRegionTracer::WriteImage(std::string filename, ImageType::Pointer image)
{
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}

void MicrogliaRegionTracer::WriteInitialMicrogliaImages()
{
	std::vector<Seed*>::iterator seeds_iter;

	//clock_t startTime = clock();
	for (seeds_iter = seeds.begin(); seeds_iter != seeds.end(); seeds_iter++)
	{
		Seed* seed = *seeds_iter;

		ImageType::SizeType roi_size;
		roi_size[0] = 100;
		roi_size[1] = 100;
		roi_size[2] = 50;
		
		//Grab the seed and its ROI
		ImageType::Pointer seed_image = roi_grabber->GetROI(seed, roi_size);

		//Make the file name of the raw seed image
		std::stringstream seed_filename_stream;
		seed_filename_stream << seed->getX() << "_" << seed->getY() << "_" << seed->getZ() << ".TIF";	//X_Y_Z.TIF
		
		//Write the seed image
		WriteImage(seed_filename_stream.str(), seed_image);

		////Multiscale LoG
		//LoG *log_obj = new LoG();
		//log_obj->RunMultiScaleLoG(seed, seed_image);
	}
	//std::cout << "Grabbed " << seeds.size() << " cells in " << (clock() - startTime)/CLOCKS_PER_SEC << " seconds" << std::endl;
}

void MicrogliaRegionTracer::Trace()
{
	std::vector<Seed*>::iterator seeds_iter;
	
	//Trace seed by seed
	for (seeds_iter = seeds.begin(); seeds_iter != seeds.end(); seeds_iter++)
	{
		Seed* seed = *seeds_iter;

		CalculateCandidatePixel(seed);
	}
}

void MicrogliaRegionTracer::CalculateCandidatePixel(Seed* seed)
{
	ImageType::SizeType roi_size;
	roi_size[0] = 100;
	roi_size[1] = 100;
	roi_size[2] = 50;
	
	//Grab the initial seedimage
	ImageType::Pointer seed_image = roi_grabber->GetROI(seed, roi_size);
	
	LoG *log_obj = new LoG();
	
	//Calculate the LoG on multiple scales and put them into a vector
	std::vector<LoGImageType::Pointer> log_seedimage_vector = log_obj->RunMultiScaleLoG(seed, seed_image);

	RidgeDetection(log_seedimage_vector);
}

void MicrogliaRegionTracer::RidgeDetection(std::vector<LoGImageType::Pointer> log_seedimage_vector)
{
	// set the diagonal terms in neighborhood iterator, this is the offsets for the diametrically opposing pixels
	itk::Offset<3>
		xp =  {{2 ,  0 ,   0}},
		xn =  {{-2,  0,    0}},
		yp =  {{0,   2,   0}},
		yn =  {{0,  -2,    0}},
		zp =  {{0,   0,    2}},
		zn =  {{0,   0,   -2}};

	//Pixel indicies relative to the top left corner of the 3x3x3 neighborhood.
	//{0, 0, 0} is the center pixel of the 3x3x3 neighborhood.
	//x: left is -1
	//y: up is -1
	//z: out of page is -1
	unsigned int
				   //{ x,    y ,  z }
		xy1 =  17, //{ 1 ,   1 ,  0 },
		xy2 =  9,  //{ -1,  -1 ,  0 },
		xy3 =  15, //{ -1,   1 ,  0 },
		xy4 =  11, //{ 1 ,  -1 ,  0 },

		yz1 =  25, //{ 0 ,   1 ,  1 },
		yz2 =  1,  //{ 0 ,  -1 , -1 },
		yz3 =  19, //{ 0 ,  -1 ,  1 },
		yz4 =  7,  //{ 0 ,   1 , -1 },

		xz1 =  23, //{ 1 ,   0 ,  1 },
		xz2 =  3,  //{-1 ,   0 , -1 },
		xz3 =  21, //{-1 ,   0 ,  1 },
		xz4 =  5;  //{ 1 ,   0 , -1 };
	
	std::vector<LoGImageType::Pointer>::iterator log_seedimage_vector_iter;
	
	//For each log image, get the critical points
	for (log_seedimage_vector_iter = log_seedimage_vector.begin(); log_seedimage_vector_iter != log_seedimage_vector.end(); log_seedimage_vector_iter++)
	{
		LoGImageType::Pointer log_image = *log_seedimage_vector_iter;

			
	}

}
