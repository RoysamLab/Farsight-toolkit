#include <cstdlib>
#include <fstream>
#include <iostream>

#include "CSVRegressionTestArgs.h"
#include "CSVRegressionTestparseArgs.h"
#include "CSVRegressionTestcompareTestBaseline.h"

// Holds the command line arguments to the CSVRegressionTest program.
int main( int argc, char *argv[] )
{
	// parse the command line arguments
	CSVRegressionTest::Args args;
	bool argsParsedValidly = parseArgs( argc, argv, args );
	if( !argsParsedValidly )
		{
		std::cerr << "Error: command line arguments were not parsed correctly." << std::endl;
		return EXIT_FAILURE;
		}
	// -h, -v
	else if( argsParsedValidly && args.TestCSVFile.empty() )
		{
	 	return EXIT_SUCCESS;
		}

	// open our files.
	std::ifstream testCSVFile( args.TestCSVFile.c_str() );
	if( !testCSVFile.is_open() )
		{
		std::cerr << "Error: Could not open the TestCSVFile." << std::endl;
		return EXIT_FAILURE;
		}
	std::ifstream baselineCSVFile( args.BaselineCSVFile.c_str() );
	if( !baselineCSVFile.is_open() )
		{
		std::cerr << "Error: Could not open the BaselineCSVFile." << std::endl;
		return EXIT_FAILURE;
		}
	std::string comparisonMessage;
	// do the comparison
	bool comparisonResult = compareTestBaseline( args, testCSVFile, baselineCSVFile, comparisonMessage );
	// close the files
	testCSVFile.close();
	baselineCSVFile.close();
	// output the result
	if( !comparisonResult )
	  {
	  std::cerr << "Error: " << comparisonMessage << std::endl;
	  return EXIT_FAILURE;
	  }
	std::cout << comparisonMessage << std::endl;
	return EXIT_SUCCESS;
}
