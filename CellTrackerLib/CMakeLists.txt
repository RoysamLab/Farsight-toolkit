
SET(cellTracker_SRCS 
MultiFrameCellTracker.cpp
helpers.cpp
ftkTrackFeatures.cpp
#TrackVoting.cpp
)

SET(cellTracker_HDRS 
helpers.h 
ftkTrackFeatures.h 
MultiFrameCellTracker.h
itkBinaryThinningImageFilter3D.h
itkBinaryThinningImageFilter3D.txx
itkBinaryThinningImageFilter3D.h
#TrackVoting.h
)
IF (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
 ENDIF (MSVC)
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

IF(WIN32)
	SET(WSOCK_WIN wsock32)
ELSE(WIN32)
	SET(WSOCK_WIN )
ENDIF(WIN32)

FIND_PACKAGE(VXL)
IF(VXL_FOUND)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ENDIF(VXL_FOUND)	

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} )

ADD_LIBRARY(CellTrackerLib 
			${cellTracker_SRCS}
			${cellTracker_HDRS}
)
IF(CPLEX_FOUND AND VXL_FOUND)
#	TARGET_LINK_LIBRARIES(CellTrackerLib ${ITK_USE_LIBRARIES} ftkImage ftkCommon ${CONCERT_LIBRARIES} ${CPLEX_LIBRARIES} ${WSOCK_WIN} ${LIBKPLS} rgtl rtvl )
	TARGET_LINK_LIBRARIES(CellTrackerLib ${ITK_USE_LIBRARIES} ftkImage ftkCommon ${CONCERT_LIBRARIES} ${CPLEX_LIBRARIES} ${WSOCK_WIN} ${LIBKPLS})
ELSE(CPLEX_FOUND AND VXL_FOUND)
	 MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
ENDIF(CPLEX_FOUND AND VXL_FOUND)
