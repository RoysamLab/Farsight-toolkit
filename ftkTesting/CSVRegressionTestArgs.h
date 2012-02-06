#ifndef _CSVRegressionTestArgs_h
#define _CSVRegressionTestArgs_h

#include <string>

namespace CSVRegressionTest
{

// Hold the command line argument parameters.
struct Args
{
  std::string TestCSVFile;
  std::string BaselineCSVFile;

  double AbsoluteTolerance;
  double FractionalTolerance;
};

} // end namespace CSVRegressionTest

#endif
