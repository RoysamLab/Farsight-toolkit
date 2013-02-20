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
	InputImageType::Pointer input = readImage<InputImageType>(argv[1]);

	int start_slice = atoi(argv[2]);
	int end_slice = atoi(argv[3]);

	InputImageType::Pointer output = InputImageType::New();
	
	InputImageType::IndexType index;
	index.Fill(0);
	InputImageType::RegionType region;
	InputImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();
	size[2] = end_slice - start_slice + 1;
	region.SetIndex(index);
	region.SetSize(size);
	output->SetRegions(region);
	output->Allocate();
	index  = input->GetLargestPossibleRegion().GetIndex();
	index[2] = start_slice;
	region.SetIndex(index);
	
	IteratorType iter(output,output->GetLargestPossibleRegion());
	ConstIteratorType iiter(input,region);
	for(iiter.GoToBegin(),iter.GoToBegin(); !iter.IsAtEnd(); ++iter,++iiter)
	{
		iter.Set(iiter.Get());
	}
	writeImage<InputImageType>(output,argv[4]);
	return 0;
}