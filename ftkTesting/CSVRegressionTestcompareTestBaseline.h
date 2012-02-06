#ifndef _CSVRegressionTestcompareTestBaseline_h
#define _CSVRegressionTestcompareTestBaseline_h

#include <fstream>

#include "CSVRegressionTestArgs.h"

namespace CSVRegressionTest
{

/** Compare the test file to the baseline file. 
 * \param args Command line arguments that provide the input parameters for the
 * comparison.
 * \param testCSV stream for the test CSV input.
 * \param baselineCSV stream for the baseline CSV input.
 * \param comparisonMessage Output string that describes the results of the
 * comparison.
 * \returns Whether comparison passes or fails.*/
bool
compareTestBaseline( const Args & args,
  std::istream & testCSV,
  std::istream & baselineCSV,
  std::string & comparisonMessage );

} // end namespace CSVRegressionTest

#endif
