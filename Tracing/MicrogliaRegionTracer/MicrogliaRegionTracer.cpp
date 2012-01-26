#include "MicrogliaRegionTracer.h"

MicrogliaRegionTracer::MicrogliaRegionTracer()
{
	roi_grabber = new ROIGrabber();
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

		//Multiscale LoG
		LoG *log_obj = new LoG();
		log_obj->RunMultiScaleLoG(seed, seed_image);
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
	std::vector<LoGImageType::Pointer> log_seedimage_vector = log_obj->RunMultiScaleLoG(seed, seed_image);

	RidgeDetection(log_seedimage_vector);
}

void MicrogliaRegionTracer::RidgeDetection(std::vector<LoGImageType::Pointer> log_seedimage_vector)
{
	std::vector<LoGImageType::Pointer>::iterator log_seedimage_vector_iter;
	
	for (log_seedimage_vector_iter = log_seedimage_vector.begin(); log_seedimage_vector_iter != log_seedimage_vector.end(); log_seedimage_vector_iter++)
	{
		LoGImageType::Pointer log_image = *log_seedimage_vector_iter;


	}

}
