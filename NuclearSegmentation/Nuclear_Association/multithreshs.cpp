#ifndef _MULTITHRESHS_CPP_
#define _MULTITHRESHS_CPP_

#include "itkImage.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

typedef unsigned short USPixelType;
typedef itk::Image< USPixelType, 2 > USImageType;

unsigned short returnthresh( USImageType::Pointer input_image, int num_bins, int num_in_fg ){

//Instantiate the different image and filter types that will be used
	typedef itk::ImageRegionIteratorWithIndex< USImageType > IteratorType;
	typedef itk::ImageRegionConstIterator< USImageType > ConstIteratorType;
	typedef itk::Statistics::ScalarImageToHistogramGenerator< USImageType > ScalarImageToHistogramGeneratorType;
	typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
	CalculatorType::Pointer calculator = CalculatorType::New();
	scalarImageToHistogramGenerator->SetNumberOfBins( 128 );
	calculator->SetNumberOfThresholds( num_bins );
	scalarImageToHistogramGenerator->SetInput( input_image);
	scalarImageToHistogramGenerator->Compute();
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	USPixelType thresh;

	for(int i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<USPixelType>(*itNum));

	return thresh;

}

#endif _MULTITHRESHS_CPP_