cmake_minimum_required(VERSION 2.8)
PROJECT(SPH)

FIND_PACKAGE(VXL)
IF(VXL_FOUND)
    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ELSEIF(VXL_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build without VXL.")
ENDIF(VXL_FOUND)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})

SET( SPH_SRCS 
	Spherical_Harmonic_transform.cpp
	)

SET( SPH_HDRS 
	Spherical_Harmonic_transform.h
	)

ADD_LIBRARY(libSPH ${SPH_SRCS} ${SPH_HDRS})
ADD_EXECUTABLE( SPH SPH_test.cpp )
TARGET_LINK_LIBRARIES(libSPH vnl vnl_algo vnl_io mbl libSPH)
  
TARGET_LINK_LIBRARIES( SPH libSPH)