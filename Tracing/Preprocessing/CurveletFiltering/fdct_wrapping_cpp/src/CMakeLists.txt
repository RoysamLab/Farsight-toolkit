PROJECT(Curvelets)

cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

INCLUDE(FindFFTW215.cmake)

IF(FFTW_FOUND)
	INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
	ADD_LIBRARY(fdct_wrapping fdct_wrapping.cpp ifdct_wrapping.cpp fdct_wrapping_param.cpp)
	TARGET_LINK_LIBRARIES(fdct_wrapping ${FFTW_LIBRARIES})
	IF(WIN32)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
	ELSE(WIN32)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -omp")
	ENDIF(WIN32)
	ADD_EXECUTABLE(curvelets curvelet_preprocessing.cpp)
	TARGET_LINK_LIBRARIES(curvelets fdct_wrapping ITKCommon ITKIO ITKStatistics)
	INSTALL(TARGETS curvelets RUNTIME DESTINATION bin)
	INSTALL(TARGETS fdct_wrapping ARCHIVE DESTINATION lib)
ENDIF(FFTW_FOUND)
