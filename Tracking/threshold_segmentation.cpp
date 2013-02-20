#include "helpers.h"

using namespace helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
double start_t,end_t,diff_t;

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
int main(int argc, char* argv[])
{
	if(argc <4)
	{
		std::cout<<"Usage: threshold_segmentation <InputImageFileName> <OutputImageFileName> <threshold,morhp_open_depth>\n";
		return 0;
	}
	InputImageType::Pointer InImage;
	LabelImageType::Pointer LabImage; 
	int inten_threshold, open_depth, approx_size;
	sscanf(argv[3],"%d,%d,%d",&inten_threshold,&open_depth);
	printf("Intensity Threshold: %d\n Morph Open Depth:%d\n ",inten_threshold,open_depth);
	if(!file_exists(argv[2]))
    {
	   InImage =  readImage<InputImageType>(argv[1]);
	   LabImage = GetThreshSegmentation(InImage,inten_threshold,open_depth);
	   LabImage = GetLargestLabel(LabImage);
	   writeImage<LabelImageType>(LabImage,argv[2]);
	}	   

}