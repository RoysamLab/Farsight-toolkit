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
	std::string outfname = argv[3];

	//if(!file_exists((char*)outfname.c_str()))
	{

	Input2DImageType16::Pointer input_image = readImage<Input2DImageType16>(input.c_str());
	
	Input2DImageType16::SpacingType spacing;
	spacing[0] = 1;
	spacing[1] = 1;
	input_image->SetSpacing(spacing);
	typedef itk::SmoothingRecursiveGaussianImageFilter<Input2DImageType16,Input2DImageType16> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(input_image);
	filter->SetSigma(atof(argv[2]));
	try
	{
		filter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}

	Input2DImageType16::Pointer backg_image = filter->GetOutput();


	typedef itk::SubtractImageFilter <Input2DImageType16, Input2DImageType16, Float2DImageType> SubtractImageFilterType;
	SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
	subtractFilter->SetInput1(input_image);
	subtractFilter->SetInput2(backg_image);
	try
	{
		subtractFilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	Float2DImageType::Pointer diff_image = subtractFilter->GetOutput();

	//writeImage<Float2DImageType>(diff_image, "C:/Lab/AminFiles/Debug/SegmentationNavin/difference.nrrd");

	Input2DImageType16::Pointer out_image = Input2DImageType16::New();	
	Input2DImageType16::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;    
	out_image->SetOrigin( origin );
	Input2DImageType16::IndexType start;
	start[0] = 0;  // first index on X
	start[1] = 0;  // first index on Y    
	Input2DImageType16::SizeType  size = input_image->GetLargestPossibleRegion().GetSize();
	Input2DImageType16::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	out_image->SetRegions( region );
	out_image->Allocate();
	out_image->FillBuffer(0);
	try
	{
		out_image->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	
	
	Float2DIteratorType diff_iter(diff_image,diff_image->GetLargestPossibleRegion());
	float min_value = std::numeric_limits<float>::max();
	for(diff_iter.GoToBegin();!diff_iter.IsAtEnd();++diff_iter)
	{
		min_value = MIN(min_value,diff_iter.Get());
	}
	std::cout<<"Found minimum value:\t"<<min_value<<std::endl;
	Iterator2DType16 out_iter(out_image,out_image->GetLargestPossibleRegion());
	if(min_value<0.0)
	{
		//std::cout<<"I am in negative loop\n"<<std::endl;
		for(diff_iter.GoToBegin(),out_iter.GoToBegin();!diff_iter.IsAtEnd(),!out_iter.IsAtEnd();++diff_iter,++out_iter)
		{
			out_iter.Set((InputPixelType16)(diff_iter.Get()-min_value));
			//std::cout<<(InputPixelType16)floor(diff_iter.Get()-min_value)<<std::endl;
		}		
	}
	else
	{
		//std::cout<<"I am in positvie  loop\n"<<std::endl;
		for(diff_iter.GoToBegin(),out_iter.GoToBegin();!diff_iter.IsAtEnd(),!out_iter.IsAtEnd();++diff_iter,++out_iter)
		{
			out_iter.Set((InputPixelType16)(diff_iter.Get()+min_value));
		}
	}
		
	writeImage<Input2DImageType16>(out_image, outfname.c_str());
	}
	return 0;
}