#include "CSVRegressionTestcompareTestBaseline.h"

// todo: remove me
#include <iostream>

#include <cmath>
#include <sstream>

namespace CSVRegressionTest
{

bool
compareTestBaseline( const Args & args,
       	std::istream & testCSV,
       	std::istream & baselineCSV,
       	std::string & comparisonMessage )
{
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

	unsigned int rowCount = 1;
	// for every line
	while( testCSV.good() )
		{
		getline( testCSV, testLine );
		getline( baselineCSV, baselineLine );
		if( testCSV.good() && !baselineCSV.good() || !testCSV.good() && baselineCSV.good() )
			{
			std::ostringstream ostrm;
			ostrm << "The Baseline does not have the same number of rows"
			      << " as the Test.";
			comparisonMessage = ostrm.str();
			return false;
			}
		testLineStream.str( testLine );
		baselineLineStream.str( baselineLine );

		unsigned int columnCount = 1;
		// for every token in a line
		while( testLineStream.good() )
			{
			// get the token
			getline( testLineStream, testToken, ',' );
			getline( baselineLineStream, baselineToken, ',' );
			if( testLineStream.good() && !baselineLineStream.good() || !testLineStream.good() && baselineLineStream.good() )
				{
				std::ostringstream ostrm;
				ostrm << "The Baseline does not have the same number of columns"
				      << " as the Test in row " << rowCount << ".";
				comparisonMessage = ostrm.str();
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
					      << " in row " << rowCount << " column " << columnCount << ".";
					comparisonMessage = ostrm.str();
					return false;
					}
				if( fabs( testTokenAsDouble - baselineTokenAsDouble ) / fabs( baselineTokenAsDouble ) >= args.FractionalTolerance )
					{
					std::ostringstream ostrm;
					ostrm <<        "Baseline entry: " << baselineTokenAsDouble << "\n"
					      << "       and Test entry: " << testTokenAsDouble << "\n"
					      << "was not within the specified fractional tolerance: " << args.FractionalTolerance << "\n"
					      << "in row " << rowCount << " column " << columnCount << ".";
					comparisonMessage = ostrm.str();
					return false;
					}
				else if( fabs( testTokenAsDouble - baselineTokenAsDouble ) >= args.AbsoluteTolerance )
					{
					std::ostringstream ostrm;
					ostrm <<        "Baseline entry: " << baselineTokenAsDouble << "\n"
					      << "       and Test entry: " << testTokenAsDouble << "\n"
					      << "was not within the specified absolute tolerance: " << args.AbsoluteTolerance << "\n"
					      << "in row " << rowCount << " column " << columnCount << ".";
					comparisonMessage = ostrm.str();
					return false;
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
					std::ostringstream ostrm;
					ostrm <<        "Baseline entry: " << baselineTokenAsString << "\n"
					      << "       and Test entry: " << testTokenAsString << "\n"
					      << "was not the same\n"
					      << "in row " << rowCount << " column " << columnCount << ".";
					comparisonMessage = ostrm.str();
					return false;
					}
				}
			++columnCount;
			}

		testLineStream.clear();
		baselineLineStream.clear();
		++rowCount;
		}

	comparisonMessage = "The Test and Baseline entries were the same within the tolerances provided.";
	return true;
}

} // end namespace CSVRegressionTest
