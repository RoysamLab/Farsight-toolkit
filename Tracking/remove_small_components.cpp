#include "helpers.h"

using namespace helpers;
template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... \n",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	ReaderType::GlobalWarningDisplayOff();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done\n");
	return reader->GetOutput();

}
template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... \n",filename);
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

int main(int argc, char**argv)
{
	printf("Running %s...",argv[0]);
	if(argc < 4)
	{
		printf(" Usage: %s input_file_name minimum_object_size output_file_name\n",argv[0]);
		std::cin.get();
		return -1;
	}

	LabelImageType::Pointer tempsegmented = readImage<LabelImageType>(argv[1]);
	tempsegmented = fillHoles(tempsegmented,2);
	//writeImage<LabelImageType>(tempsegmented,"C:/Users/arun/Research/Farsight/exe/bin/remove_small_components_debug1.tif");
	tempsegmented = getLargeLabels(tempsegmented,atoi(argv[2]));
	writeImage<LabelImageType>(tempsegmented,argv[3]);
	printf("Done\n");
	return 0;
}