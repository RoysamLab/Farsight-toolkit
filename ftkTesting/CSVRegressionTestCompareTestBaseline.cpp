#include "CSVRegressionTestCompareTestBaseline.h"

// todo: remove me
#include <iostream>

#include <cmath>
#include <sstream>

#include "itkImageLinearIteratorWithIndex.h"


namespace CSVRegressionTest
{

CompareTestBaseline
::CompareTestBaseline()
{
	this->PNGImageWriter = PNGImageWriterType::New();
	this->MinMaxCalculator = MinMaxCalculatorType::New();
	this->IntensityWindowingFilter = IntensityWindowingFilterType::New();
	this->ResampleFilter = ResampleFilterType::New();
	this->ResampleInterpolator = ResampleInterpolatorType::New();
	this->ResampleFilter->SetInterpolator( this->ResampleInterpolator );
	this->NumericalErrorWriter = NumericalErrorWriterType::New();
}

const std::string &
CompareTestBaseline
::GetComparisonMessage() const
{
	return this->ComparisonMessage;
}

void
CompareTestBaseline
::SetArgs( const Args * args )
{
	this->Parameters = args;
}

bool
CompareTestBaseline
::DoComparison( std::istream & testCSV,
	std::istream & baselineCSV )
{
	this->StringError.clear();
	this->AbsoluteError.clear();
	this->FractionalError.clear();
	this->RowCount = 0;
	this->MaxColumnCount = 1;
	this->StringErrorCount = 0;
	this->AbsoluteErrorCount = 0;
	this->FractionalErrorCount = 0;

	// the stream gets partitioned into lines, then tokens, then type
	// converted to a double or white-space removed string.
	std::string testLine;
	std::string baselineLine;
	std::istringstream testLineStream;
	std::istringstream baselineLineStream;
	std::string testToken;
	std::string baselineToken;
	std::istringstream testTokenStream;
	std::istringstream baselineTokenStream;
	double testTokenAsDouble;
	double baselineTokenAsDouble;
	std::string testTokenAsString;
	std::string baselineTokenAsString;
	bool fractionalErrorOccured = false;
	bool absoluteErrorOccurred = false;
	bool stringErrorOccurred = false;

	// for every line
	while( testCSV.good() )
		{
		++this->RowCount;
		StringErrorRowType stringErrorRow;
		NumericalErrorRowType fractionalErrorRow;
		NumericalErrorRowType absoluteErrorRow;

		getline( testCSV, testLine );
		getline( baselineCSV, baselineLine );
		if( testCSV.good() && !baselineCSV.good() || !testCSV.good() && baselineCSV.good() )
			{
			std::ostringstream ostrm;
			ostrm << "The Baseline does not have the same number of rows"
			      << " as the Test.";
			this->ComparisonMessage = ostrm.str();
			return false;
			}
		testLineStream.str( testLine );
		baselineLineStream.str( baselineLine );

		IndexValueType columnCount = 0;
		// for every token in a line
		while( testLineStream.good() )
			{
			++columnCount;
			// get the token
			getline( testLineStream, testToken, this->Parameters->Delimiter );
			getline( baselineLineStream, baselineToken, this->Parameters->Delimiter );
			if( testLineStream.good() && !baselineLineStream.good() || !testLineStream.good() && baselineLineStream.good() )
				{
				std::ostringstream ostrm;
				ostrm << "The Baseline does not have the same number of columns"
				      << " as the Test in row " << this->RowCount << ".";
				this->ComparisonMessage = ostrm.str();
				return false;
				}

			testTokenStream.clear();
			baselineTokenStream.clear();
			testTokenStream.str( testToken );
			baselineTokenStream.str( baselineToken );

			baselineTokenStream >> baselineTokenAsDouble;
			if( !baselineTokenStream.fail() )
				{
				testTokenStream >> testTokenAsDouble;
				if( testTokenStream.fail() )
					{
					std::ostringstream ostrm;
					ostrm << "Test entry was not a number when Baseline was a number"
					      << " in row " << this->RowCount << " column " << columnCount << ".";
					this->ComparisonMessage = ostrm.str();
					return false;
					}
				double fractionalError = fabs( testTokenAsDouble - baselineTokenAsDouble ) / fabs( baselineTokenAsDouble );
				if( fractionalError >= this->Parameters->FractionalTolerance )
					{
					NumericalErrorPairType error( columnCount - 1, fractionalError );
					fractionalErrorRow.push_back( error );
					fractionalErrorOccured = true;
					++this->FractionalErrorCount;
					}
				double absoluteError = fabs( testTokenAsDouble - baselineTokenAsDouble ); 
				if( absoluteError >= this->Parameters->AbsoluteTolerance )
					{
					NumericalErrorPairType error( columnCount - 1, absoluteError );
					absoluteErrorRow.push_back( error );
					absoluteErrorOccurred = true;
					++this->AbsoluteErrorCount;
					}
				}
			else // not a number, so treat it as a string
				{
				baselineTokenStream.seekg( 0 );
				baselineTokenStream.clear();
				testTokenAsString = "";
				baselineTokenAsString = "";
				// Initialize with non-empty string so we can
				// skip trailing whitespace in the while loop.
				std::string testWord;
				std::string baselineWord;
				while( !baselineTokenStream.eof() )
					{
					testTokenStream >> testWord;
					baselineTokenStream >> baselineWord;
					// trailing whitespace
					if( baselineTokenStream.fail() )
						{
						continue;
						}
					testTokenAsString += testWord;
					baselineTokenAsString += baselineWord;
					}
				if( baselineTokenAsString.compare( testTokenAsString ) )
					{
					stringErrorRow.push_back( columnCount -1 );
					stringErrorOccurred = true;
					++this->StringErrorCount;
					}
				}
			}

		if( columnCount > this->MaxColumnCount )
			{
			this->MaxColumnCount = columnCount;
			}
		this->StringError.push_back( stringErrorRow );
	       	this->FractionalError.push_back( fractionalErrorRow );
		this->AbsoluteError.push_back( absoluteErrorRow );
		testLineStream.clear();
		baselineLineStream.clear();
		}
	--this->RowCount;

	std::ostringstream ostrm;
	if( fractionalErrorOccured )
		{
		ostrm << "All numerical values were not within the specified fractional tolerance: " 
			<< this->Parameters->FractionalTolerance << std::endl;
		this->WriteErrorImage< NumericalErrorImageType >( "FractionalError" );
		}
	if( absoluteErrorOccurred )
		{
		ostrm << "All numerical values were not within the specified absolute tolerance: " <<
		       	this->Parameters->AbsoluteTolerance << std::endl;
		this->WriteErrorImage< NumericalErrorImageType >( "AbsoluteError" );
		}
	if( stringErrorOccurred )
		{
		ostrm << "Not all string entries were the same." << std::endl;
		this->WriteErrorImage< StringErrorImageType >( "StringError" );
		}


	if( fractionalErrorOccured || absoluteErrorOccurred || stringErrorOccurred )
		{
		this->ComparisonMessage = ostrm.str();
		return false;
		}
	this->ComparisonMessage = "The Test and Baseline entries were the same within the tolerances provided.";
	return true;
}

template< class TErrorImage >
void
CompareTestBaseline
::WriteErrorImage( const std::string & error )
{
	std::string filename = this->Parameters->OutputFilePrefix;
	PNGImageToWriteType::Pointer errorImageToWrite;
	
	if( !error.compare( "StringError" ) )
		{
		// Create the image.
		typedef StringErrorImageType::RegionType RegionType;
		RegionType::SizeType size;
		size[0] = this->MaxColumnCount;
		size[1] = this->RowCount;
		StringErrorImageType::RegionType regions;
		regions.SetSize( size );
		errorImageToWrite = StringErrorImageType::New();
		errorImageToWrite->SetRegions( regions );
		errorImageToWrite->Allocate();

		typedef itk::ImageLinearIteratorWithIndex< StringErrorImageType > ImageIteratorType;
		ImageIteratorType imageIt( errorImageToWrite, regions );
		StringErrorType::const_iterator stringErrorIt;
		StringErrorRowType::const_iterator stringErrorRowIt;
		for( imageIt.GoToBegin(), stringErrorIt = this->StringError.begin();
			stringErrorIt != this->StringError.end();
			++stringErrorIt )
			{
			imageIt.GoToBeginOfLine();
			stringErrorRowIt = (*stringErrorIt).begin();
			while( stringErrorRowIt != (*stringErrorIt).end() )
				{
				const IndexValueType nextErrorIndex = *stringErrorRowIt;
				while( imageIt.GetIndex()[0] != nextErrorIndex )
					{
					imageIt.Set( 0 );
					++imageIt;
					}
				imageIt.Set( 255 );
				++imageIt;
				++stringErrorRowIt;
				}
			while( !imageIt.IsAtEndOfLine() )
				{
				imageIt.Set( 0 );
				++imageIt;
				}
			imageIt.NextLine();
			}
		filename += "_StringError";
		this->PNGImageWriter->SetInput( errorImageToWrite );
		this->PNGImageWriter->SetFileName( filename + ".png" );
		this->PNGImageWriter->Update();
		// This gets it deposited on the CDash dashboard.
		std::cout << "<DartMeasurement name=\"StringErrorCount\" type=\"numeric/double\">";
		std::cout << this->StringErrorCount;
		std::cout << "</DartMeasurement>" << std::endl;
		}
	else if( !error.compare( "FractionalError" ) )
		{
		filename += "_FractionalError";
		std::cout << "<DartMeasurement name=\"FractionalErrorCount\" type=\"numeric/double\">";
		std::cout << this->FractionalErrorCount;
		std::cout << "</DartMeasurement>" << std::endl;
		}
	else if( !error.compare( "AbsoluteError" ) )
		{
		filename += "_AbsoluteError";
		std::cout << "<DartMeasurement name=\"AbsoluteErrorCount\" type=\"numeric/double\">";
		std::cout << this->AbsoluteErrorCount;
		std::cout << "</DartMeasurement>" << std::endl;
		}
	else
		{
		throw std::logic_error( "WriteErrorImage: Unknown error type." );
		}

	if( !error.compare( "AbsoluteError" ) || !error.compare( "FractionalError" ) )
		{
		// Create the image.
		typedef NumericalErrorImageType::RegionType RegionType;
		RegionType::SizeType size;
		size[0] = this->MaxColumnCount;
		size[1] = this->RowCount;
		NumericalErrorImageType::RegionType regions;
		regions.SetSize( size );
		NumericalErrorImageType::Pointer numericalErrorImage = NumericalErrorImageType::New();
		numericalErrorImage->SetRegions( regions );
		numericalErrorImage->Allocate();

		typedef itk::ImageLinearIteratorWithIndex< NumericalErrorImageType > ImageIteratorType;
		ImageIteratorType imageIt( numericalErrorImage, regions );
		NumericalErrorType::const_iterator numericalErrorIt;
		NumericalErrorType::const_iterator numericalErrorItEnd;
		NumericalErrorRowType::const_iterator numericalErrorRowIt;
		if( !error.compare( "AbsoluteError" ) )
			{
			numericalErrorIt = this->AbsoluteError.begin();
			numericalErrorItEnd = this->AbsoluteError.end();
			}
		else
			{
			numericalErrorIt = this->FractionalError.begin();
			numericalErrorItEnd = this->FractionalError.end();
			}
		for( imageIt.GoToBegin(); 
			numericalErrorIt != numericalErrorItEnd;
			++numericalErrorIt )
			{
			imageIt.GoToBeginOfLine();
			numericalErrorRowIt = (*numericalErrorIt).begin();
			while( numericalErrorRowIt != (*numericalErrorIt).end() )
				{
				const NumericalErrorPairType nextErrorPair = *numericalErrorRowIt;
				while( imageIt.GetIndex()[0] != nextErrorPair.first )
					{
					imageIt.Set( 0.0 );
					++imageIt;
					}
				imageIt.Set( nextErrorPair.second );
				++imageIt;
				++numericalErrorRowIt;
				}
			while( !imageIt.IsAtEndOfLine() )
				{
				imageIt.Set( 0.0 );
				++imageIt;
				}
			imageIt.NextLine();
			}

		this->NumericalErrorWriter->SetInput( numericalErrorImage );
		this->NumericalErrorWriter->SetFileName( filename + ".mha" );
		this->NumericalErrorWriter->Update();

		this->MinMaxCalculator->SetImage( numericalErrorImage );
		this->MinMaxCalculator->Compute();

		this->IntensityWindowingFilter->SetInput( numericalErrorImage );
		this->IntensityWindowingFilter->SetWindowMinimum( this->MinMaxCalculator->GetMinimum() );
		this->IntensityWindowingFilter->SetWindowMaximum( this->MinMaxCalculator->GetMaximum() );
		this->IntensityWindowingFilter->SetOutputMinimum( 0 );
		this->IntensityWindowingFilter->SetOutputMaximum( 255 );
		this->IntensityWindowingFilter->Update();
		errorImageToWrite = this->IntensityWindowingFilter->GetOutput();
		}

	// Write a PNG version resampled to 640x480, so it is easy to visualize
	// out of the box.  The location of the error can easily be identified
	// in the other images using the index of the erroroneous pixel.
	this->ResampleFilter->SetInput( errorImageToWrite );
	typedef PNGImageToWriteType::RegionType  RegionType;
	typedef RegionType::SizeType 		SizeType;
	SizeType outputSize;
	outputSize[0] = 640;
	outputSize[1] = 480;
	this->ResampleFilter->SetSize( outputSize );
	typedef PNGImageToWriteType::SpacingType SpacingType;
	SpacingType outputSpacing;
	const SizeType inputSize = errorImageToWrite->GetLargestPossibleRegion().GetSize();
	outputSpacing[0] = static_cast< double >( inputSize[0] - 1 ) / static_cast< double >( outputSize[0] - 1 );
	outputSpacing[1] = static_cast< double >( inputSize[1] - 1 ) / static_cast< double >( outputSize[1] - 1 );
	this->ResampleFilter->SetOutputSpacing( outputSpacing );
	this->PNGImageWriter->SetInput( this->ResampleFilter->GetOutput() );
	this->PNGImageWriter->SetFileName( filename + "_Resampled.png" );
	this->PNGImageWriter->Update();
	std::string errorImage = error + "Image";
	std::cout << "<DartMeasurementFile name=\"" << errorImage << "\" type=\"image/png\">";
	std::cout << filename + "_Resampled.png";
	std::cout << "</DartMeasurementFile>" << std::endl;
}

} // end namespace CSVRegressionTest
