#include "helpers.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"
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



#define MAX_TIME 1000
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
	std::string input = argv[1];
	std::string out = argv[3];

	if(!file_exists((char*)out.c_str())){

	InputImageType::Pointer im = readImage<InputImageType>(input.c_str());
	
	InputImageType::SpacingType spacing;
	spacing[0] = 1;
	spacing[1] = 1;
	spacing[2] = 4.0;
	im->SetSpacing(spacing);
	typedef itk::SmoothingRecursiveGaussianImageFilter<InputImageType,InputImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(im);
	filter->SetSigma(atof(argv[2]));
	
	
	filter->Update();

	InputImageType::Pointer outim = filter->GetOutput();

	IteratorType it1(im,im->GetLargestPossibleRegion());
	IteratorType it2(outim,outim->GetLargestPossibleRegion());
	for(it1.GoToBegin(),it2.GoToBegin();!it1.IsAtEnd();++it1,++it2)
	{
		it2.Set(MIN(MAX(it1.Get()-it2.Get(),0),255));
	}
	writeImage<InputImageType>(outim, out.c_str());
	}
	return 0;
}
int main_old(int argc, char** argv)
{
	int num_t = (argc-1)/2;
	InputImageType::Pointer input[MAX_TIME];
	InputImageType::Pointer output = InputImageType::New();

	if(num_t>MAX_TIME)
	{
		printf("Too many files\n");
		return 0;
	}
	FloatImageType::Pointer mean_image = FloatImageType::New();
	typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
	for(int counter=0; counter< num_t; counter++)
	{
		input[counter] = readImage<InputImageType>(argv[counter+1]);
		printf("Finished reading\n");
		ConstIteratorType input_iter(input[counter],input[counter]->GetLargestPossibleRegion());
		printf("entering if\n");
		if(counter==0)
		{
			mean_image->SetRegions(input[counter]->GetLargestPossibleRegion());
			mean_image->Allocate();
			mean_image->FillBuffer(0);
			printf("Allocated\n");
		}
		printf("before copy\n");
		FloatIteratorType iter(mean_image,mean_image->GetLargestPossibleRegion());
		for(iter.GoToBegin(),input_iter.GoToBegin();!iter.IsAtEnd(); ++iter,++input_iter)
		{
			iter.Set(iter.Get()+input_iter.Get()*1.0/num_t);
		}
	}
	

	output->SetRegions(mean_image->GetLargestPossibleRegion());
	output->Allocate();
	
	for(int counter=0; counter<num_t; counter++)
	{

		typedef itk::ImageRegionConstIterator<FloatImageType> ConstFloatIteratorType;
		ConstFloatIteratorType cfiter(mean_image, mean_image->GetLargestPossibleRegion());
		ConstIteratorType input_iter(input[counter], input[counter]->GetLargestPossibleRegion());
		IteratorType iter(output,output->GetLargestPossibleRegion());
		iter.GoToBegin();
		cfiter.GoToBegin();
		input_iter.GoToBegin();
		for(;!iter.IsAtEnd();++iter,++cfiter,++input_iter)
		{
			int value = input_iter.Get() - cfiter.Get();
			if(value <0)
				value = 0;
			iter.Set(value);
		}
		writeImage<InputImageType>(output,argv[counter+num_t+1]);
	}
}