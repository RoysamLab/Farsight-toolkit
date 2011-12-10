CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK")
ENDIF(VTK_FOUND)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})


################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_SRCS

main.cpp
#ftkNuclearSegmentationNic.cpp
 
)

################################################################################
################################################################################


SET( FTKNUCLEARSEGMENTATIONNIC_HDRS

ftkNuclearSegmentationNic.h
ftkNuclearSegmentationNic.hxx
 
)

################################################################################
################################################################################


ADD_LIBRARY( libftkNuclearSegmentationNic ${FTKNUCLEARSEGMENTATIONNIC_HDRS} ${FTKNUCLEARSEGMENTATIONNIC_SRCS})

ADD_EXECUTABLE( ftkNuclearSegmentationNic
		main.cpp
)

################################################################################
################################################################################


TARGET_LINK_LIBRARIES ( ftkNuclearSegmentationNic ${VTK_LIBRARIES} ${ITK_LIBRARIES} libftkNuclearSegmentationNic )
