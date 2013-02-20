#include <iostream>
#include <math.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 3 > USImageType;
typedef float FloatPixelType;

int main(int argc, char *argv[])
{
	if( argc < 3 ){
		std::cout << "Usage: " <<(argv[0]) << " ";
		std::cerr << "ROI_Image\t Image_To_Be_Processed\t Out_Image\t Num_Thresholds\t Num_In_Fg \n";
		return EXIT_FAILURE;
	}

	typedef itk::ImageFileReader< USImageType > ReaderType;
	typedef itk::ImageFileWriter< USImageType > WriterType;
	typedef itk::ImageRegionIterator< USImageType > IteratorType; 
	typedef itk::ImageRegionConstIterator< USImageType > ConstIteratorType;
	typedef itk::Statistics::Histogram< FloatPixelType > HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	int num_bin_levs, num_in_fg;
	if( argc > 4 )
	{
		num_bin_levs = atoi(argv[4]);
		num_in_fg    = atoi(argv[5]);
	}
	else
	{
		num_bin_levs = 1;
		num_in_fg    = 1;
	}

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );
	reader->Update();

	USImageType::Pointer input_image=reader->GetOutput();

	std::cout<<"Starting threshold computation\n";
	//Create a temporary histogram container:
	const int numBins = 256;
	double tempHist[numBins];
	for(int i=0; i<numBins; ++i)
	{
		tempHist[i] = 0;
	}

	//Populate the histogram (assume pixel type is actually uchar):
	ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		USPixelType pix = it.Get();
		if(pix <= 255)
		{
			++tempHist[pix];
		}
	}

	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<unsigned short>::max();
	double max = 0.0;
	for(int i=0; i<numBins; ++i)
	{
		if( tempHist[i] > max )
			max = tempHist[i];
	}

	double scaleFactor = 1;
	if(max >= floatIntegerMax)
	{
		scaleFactor = floatIntegerMax / max;
	}

	HistogramType::Pointer histogram = HistogramType::New() ;
	// initialize histogram
	HistogramType::SizeType size;
	HistogramType::MeasurementVectorType lowerBound ;
	HistogramType::MeasurementVectorType upperBound ;

	lowerBound[0] = 0.0;
	upperBound[0] = (float)255.0;
	size.Fill(numBins);

	histogram->Initialize(size, lowerBound, upperBound ) ;

	int i=0;
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter )
	{
		float norm_freq = (float)(tempHist[i] * scaleFactor);
		iter.SetFrequency(norm_freq);
		++i;
	}

	std::cout<<"Histogram computed\n";

	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bin_levs );
	calculator->SetInputHistogram( histogram );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	FloatPixelType thresh;

	for(int i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<FloatPixelType>(*itNum));

	std::cout<<"Threshold computed: "<<thresh<<std::endl;

	int threshold = (int)(thresh+0.5);

	ReaderType::Pointer reader1 = ReaderType::New();
	reader1->SetFileName( argv[2] );
	reader1->Update();

	USImageType::Pointer input_image1 = reader1->GetOutput();

	IteratorType iterator1(input_image, input_image->GetRequestedRegion());
	IteratorType iterator2(input_image1, input_image1->GetRequestedRegion());
	iterator1.GoToBegin(); iterator2.GoToBegin();

	while( !iterator1.IsAtEnd() || !iterator2.IsAtEnd() )
	{
		if( iterator1.Get() > threshold )
			iterator2.Set( 0 );
		++iterator1;
		++iterator2;
	}

	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( input_image1 );
	writer->SetFileName( argv[3] );

	 try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	return EXIT_SUCCESS;
}
