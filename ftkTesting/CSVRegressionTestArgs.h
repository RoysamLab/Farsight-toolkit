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

  std::string OutputFilePrefix;

  double AbsoluteTolerance;
  double FractionalTolerance;

  char Delimiter;
};

} // end namespace CSVRegressionTest

#endif
