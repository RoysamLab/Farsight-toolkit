
SET(cellTracker_SRCS 
MultiFrameCellTracker.cpp
helpers.cpp
)

SET(cellTracker_HDRS 
helpers.h 
ftkTrackFeatures.h 
MultiFrameCellTracker.h
itkBinaryThinningImageFilter3D.h
itkBinaryThinningImageFilter3D.txx
itkBinaryThinningImageFilter3D.h
)

FIND_PACKAGE(CPLEX)
FIND_PACKAGE(CONCERT)
IF(CPLEX_FOUND)
	INCLUDE_DIRECTORIES(${CPLEX_INCLUDE_DIRS})
ENDIF(CPLEX_FOUND)
IF(CONCERT_FOUND)
	INCLUDE_DIRECTORIES(${CONCERT_INCLUDE_DIRS})
ENDIF(CONCERT_FOUND)
ADD_DEFINITIONS(-DIL_STD)
find_package( Boost)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} )

ADD_LIBRARY(CellTrackerLib 
			${cellTracker_SRCS}
			${cellTracker_HDRS}
)
IF(CPLEX_FOUND)
	TARGET_LINK_LIBRARIES(CellTrackerLib ITKCommon ITKIO  ITKStatistics ftkImage ${CPLEX_LIBRARIES} ${CONCERT_LIBRARIES} wsock32 libKPLS)
ENDIF(CPLEX_FOUND)