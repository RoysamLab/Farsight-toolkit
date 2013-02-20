#include "helpers.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"


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

bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	return false;
}

int main(int argc, char**argv)
{
	typedef itk::GradientAnisotropicDiffusionImageFilter<FloatImageType, FloatImageType> FilterType;

	if(!file_exists(argv[3]))
	{	
	InputImageType::Pointer im = readImage<InputImageType>(argv[1]);

	MedianFilterType::Pointer mfilter = MedianFilterType::New();
	InputImageType::SizeType radius;
	int rx,ry,rz;
	sscanf(argv[2],"%d,%d,%d",&rx,&ry,&rz);
	radius[0] = rx;
	radius[1] = ry;
	radius[2] = rz;
	printf("doing a median filter of size %d, %d, %d, on |%s| to give |%s|\n",rx,ry,rz, argv[1],argv[3]);
	mfilter->SetRadius(radius);
	mfilter->SetInput(im);
	mfilter->Update();
	writeImage<InputImageType>(mfilter->GetOutput(),argv[3]);
	}
	return 0;
}
