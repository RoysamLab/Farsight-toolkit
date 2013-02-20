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
	typedef itk::GradientAnisotropicDiffusionImageFilter<FloatImageType, FloatImageType> FilterType;

	if(!file_exists(argv[3]))
	{
		InputImageType::Pointer im = readImage<InputImageType>(argv[1]);

		typedef itk::CastImageFilter<InputImageType,FloatImageType> CastFilter;
		typedef itk::CastImageFilter<FloatImageType,InputImageType> ReverseCastFilter;
		FilterType::Pointer filter = FilterType::New();


		CastFilter::Pointer cfilter = CastFilter::New();
		cfilter->SetInput(im);
		filter->SetInput(cfilter->GetOutput());

		int num_i;
		float time_step;
		float conductance;
		sscanf(argv[2],"%d,%f,%f",&num_i,&time_step,&conductance);
		printf("time_step %f num_i %d conductance %f\n",time_step,num_i,conductance);
		filter->SetTimeStep(time_step);
		filter->SetNumberOfIterations(num_i);
		filter->SetConductanceParameter(conductance);
		filter->Update();

		FloatImageType::Pointer imf = filter->GetOutput();
		typedef itk::ImageRegionIterator<FloatImageType> IRI;
		IRI imageIterator(imf,imf->GetLargestPossibleRegion());
		for(imageIterator.GoToBegin();!imageIterator.IsAtEnd(); ++imageIterator)
		{
			float value = imageIterator.Get();
			value = ((value < 0)?0:((value>255)?255:value));
			imageIterator.Set(value);
		}

		ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
		rfilter->SetInput(filter->GetOutput());
		rfilter->Update();
		writeImage<InputImageType>(rfilter->GetOutput(),argv[3]);
	}

	return 0;
}