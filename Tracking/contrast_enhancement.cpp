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

bool file_exists(char *filename)
{
	FILE *fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	else
		return false;
	
}

int main(int argc, char**argv)
{


	//if(!file_exists(argv[3]))
	//{
		typedef unsigned char InPixType;
	    typedef itk::Image<InPixType ,3> InImageType;
		InImageType::Pointer imdata = readImage<InImageType>(argv[1]);

		InPixType min_intens;
		InPixType max_intens;
		sscanf(argv[2],"%u,%u",&min_intens,&max_intens);
		printf("min_intens %u max_intens %u\n",min_intens,max_intens);
		scanf("%*d",&min_intens);

		// set up the threshold ratio:
		InPixType old_range = max_intens-min_intens;
		InPixType new_range = std::numeric_limits<InPixType>::max();
		float range_ratio = ((float)new_range)/((float)old_range);

		typedef itk::ImageRegionIterator<InImageType> InImageRegIterator;
		InImageRegIterator imageIterator(imdata,imdata->GetLargestPossibleRegion());
		for(imageIterator.GoToBegin();!imageIterator.IsAtEnd(); ++imageIterator)
		{
			InPixType pixval = imageIterator.Get();
			if(pixval<min_intens)
			{
				pixval = min_intens;
			}
			else if(pixval>max_intens)
			{
				pixval  = max_intens;
			}
			pixval = (InPixType)((pixval-min_intens)*range_ratio);
			pixval = ((pixval < 0)?0:((pixval>new_range)?new_range:pixval));
			imageIterator.Set(pixval);
		}

		writeImage<InImageType>(imdata,argv[3]);
	//}

	return 0;
}