#ifndef _CSVRegressionTestcompareTestBaseline_h
#define _CSVRegressionTestcompareTestBaseline_h

#include <fstream>
#include <vector>
#include <utility>

#include "itkIntTypes.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "CSVRegressionTestArgs.h"

namespace CSVRegressionTest
{

/** Compare the test file to the baseline file. */
class CompareTestBaseline
{
public:
	/** Run the comparison of the test data to the baseline data.
 	* \param testCSV stream for the test CSV input.
	* \param baselineCSV stream for the baseline CSV input.
	* \returns Whether comparison passes or fails. */
	bool DoComparison( std::istream & testCSV, std::istream & baselineCSV );

	/** Get the message that describes the result a comparison made with
	 * DoComparison. */
	const std::string & GetComparisonMessage() const;

	/** Set the arguments that parameterize the comparison. */
	void SetArgs( const Args * args );

	CompareTestBaseline();

private:
	/** Write the given error type, one of FractionalError, AbsoluteError,
	 * and StringError. */
	template< class TErrorImage >
	void WriteErrorImage( const std::string & error );

	std::string  ComparisonMessage;
	const Args * Parameters;

	typedef itk::IndexValueType IndexValueType;
	typedef itk::SizeValueType  SizeValueType;

	// For each string entry, we just identify when an error occurred
	// instead of what error, occurred.
	typedef std::vector< IndexValueType >                StringErrorRowType;
	typedef std::vector< StringErrorRowType >            StringErrorType;
	typedef itk::Image< unsigned char, 2 >               StringErrorImageType;

	StringErrorType     StringError;

	typedef std::pair< IndexValueType, double >             NumericalErrorPairType;
	typedef std::vector< NumericalErrorPairType >           NumericalErrorRowType;
	typedef std::vector< NumericalErrorRowType >            NumericalErrorType; 
	typedef itk::Image< double, 2 >                         NumericalErrorImageType;

	NumericalErrorType  AbsoluteError;
	NumericalErrorType  FractionalError;


	typedef StringErrorImageType                          			PNGImageToWriteType;
	typedef itk::ImageFileWriter< PNGImageToWriteType > 			PNGImageWriterType;
	typedef itk::MinimumMaximumImageCalculator< NumericalErrorImageType >   MinMaxCalculatorType;
	typedef itk::IntensityWindowingImageFilter< NumericalErrorImageType, PNGImageToWriteType >
		IntensityWindowingFilterType;
	typedef itk::ResampleImageFilter< PNGImageToWriteType, PNGImageToWriteType >
		ResampleFilterType;
	typedef itk::NearestNeighborInterpolateImageFunction< PNGImageToWriteType >
		ResampleInterpolatorType;
	typedef itk::ImageFileWriter< NumericalErrorImageType >                 NumericalErrorWriterType;

	PNGImageWriterType::Pointer PNGImageWriter;
	MinMaxCalculatorType::Pointer MinMaxCalculator;
	IntensityWindowingFilterType::Pointer IntensityWindowingFilter;
	ResampleFilterType::Pointer ResampleFilter;
	ResampleInterpolatorType::Pointer ResampleInterpolator;
	NumericalErrorWriterType::Pointer NumericalErrorWriter;

	IndexValueType RowCount;
	IndexValueType MaxColumnCount;

	SizeValueType  StringErrorCount;
	SizeValueType  AbsoluteErrorCount;
	SizeValueType  FractionalErrorCount;
};

} // end namespace CSVRegressionTest

#endif
