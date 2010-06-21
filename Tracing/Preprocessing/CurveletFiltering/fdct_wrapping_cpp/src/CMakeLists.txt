PROJECT(Curvelets)

    if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
    endif(COMMAND cmake_policy)

cmake_minimum_required(VERSION 2.0)

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
	ADD_EXECUTABLE(test test.cpp)
	ADD_EXECUTABLE(curvelets image_denoising.cpp)
	TARGET_LINK_LIBRARIES(test fdct_wrapping)
	TARGET_LINK_LIBRARIES(curvelets fdct_wrapping ITKCommon ITKIO ITKStatistics)
ENDIF(FFTW_FOUND)
