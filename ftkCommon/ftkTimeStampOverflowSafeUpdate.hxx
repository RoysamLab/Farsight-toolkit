#ifndef FTK_TIMESTAMP_OVERFLOW_SAFE_UPDATE_HXX
#define FTK_TIMESTAMP_OVERFLOW_SAFE_UPDATE_HXX

#include "ftkTimeStampOverflowSafeUpdate.h"

#include "itkMacro.h"

namespace ftk
{
	template< class TFilter >
	bool TimeStampOverflowSafeUpdate( TFilter * filter, const unsigned int numberOfRetrys )
	{
		bool overflow_exception_happened = false;

		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			if( numberOfRetrys == 0 )
			{
				std::cerr << "Number of retries in TimeStampOverflowSafeUpdate exceeded.  Re-throwing exception." << std::endl;
				throw e;
			}
			TimeStampOverflowSafeUpdate( filter, numberOfRetrys - 1 );
			overflow_exception_happened = true;
		}

		return overflow_exception_happened;
	}

} // end namespace ftk

#endif
