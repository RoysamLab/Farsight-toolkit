#ifndef _CSVRegressionTestparseArgs_h
#define _CSVRegressionTestparseArgs_h

#include "CSVRegressionTestArgs.h"

namespace CSVRegressionTest
{

/** Parse the command line arguments. Return the same value as
 * MetaCommand::Parse(). */
bool
parseArgs( int argc, char * argv[], Args & args );

} // end namespace CSVRegressionTest

#endif
