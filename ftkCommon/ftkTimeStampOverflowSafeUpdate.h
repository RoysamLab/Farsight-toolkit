#ifndef FTK_TIMESTAMP_OVERFLOW_SAFE_UPDATE_H
#define FTK_TIMESTAMP_OVERFLOW_SAFE_UPDATE_H

namespace ftk
{

/** \brief This function calls Update() on a filter in a safe way under
 * long-running, threaded circumstances.
 *
 * ITK Objects keep track of a ModifiedTime.  The ModifiedTime is updated
 * whenever the object changes.  The value is taken from a global integer
 * counter (see itk::TimeStamp), which is incremented on every access so the
 * ModifiedTime is increased when the next Object has Modified() called.  When
 * pipeline processing is initiated with an Update() call, the ModifiedTime's of
 * the Objects are used to determine which parts of a processing pipeline are
 * out-of-date and need to be calculated.
 *
 * This system may encounter problems when the global integer counter overflows
 * when running an very extended analysis.  When only a single filter is being
 * executed, the glober counter rollover may simply cause re-execution of the
 * pipeline.  When multiple pipelines are being executed simultaneously in
 * different threads, however, and those pipelines contain multiple filters or
 * filters with an internal mini-pipeline, errors may result, such as
 *
 *    Requested region is (at least partially) outside the largest possible
 *    region.
 *
 * As a workaround, this function will catch these exceptions and repeatedly try
 * to call Update() for a prescribed number of times so the filter has a chance
 * to Update() with meaningful ModifiedTimes.  If this fails, the exception is
 * re-thrown.
 * 
 * This function returns whether or not the overflow induced exception happened.
 */
template< class TFilter >
bool TimeStampOverflowSafeUpdate( TFilter * filter, const unsigned int numberOfRetrys = 10 );

} // end namespace ftk

#include "ftkTimeStampOverflowSafeUpdate.hxx"

#endif
