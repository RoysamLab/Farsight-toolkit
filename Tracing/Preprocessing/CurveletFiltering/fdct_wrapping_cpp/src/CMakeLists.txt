
cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

include(CTest)
include(ExternalProject)
set(EXTERNAL_BASE "${CMAKE_BINARY_DIR}/ExternalProjects")
set_property(DIRECTORY PROPERTY EP_BASE ${EXTERNAL_BASE})

if(CMAKE_EXTRA_GENERATOR)
  set(EXTERNAL_GEN "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else(CMAKE_EXTRA_GENERATOR)
  set(EXTERNAL_GEN "${CMAKE_GENERATOR}")
endif(CMAKE_EXTRA_GENERATOR)

# Set the default build type---this will affect all libraries and
# applications
#
set(BUILD_TYPE "Release")
if(CMAKE_BUILD_TYPE)
  set(BUILD_TYPE "${CMAKE_BUILD_TYPE}")
endif(CMAKE_BUILD_TYPE)

set(INSTALL_DIR "${EXTERNAL_BASE}/Install")

ExternalProject_Add(fftw-2.1.5
  URL http://dl.dropbox.com/u/26629462/fftw-2.1.5.tar.gz
  URL_MD5 1a89d7034071875ff703144fda121e9a
  CMAKE_GENERATOR ${EXTERNAL_GEN}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}/fftw-2.1.5
    -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}
)

SET( FFTW215_INSTALL_DIR "${INSTALL_DIR}/fftw-2.1.5" )

PROJECT(Curvelets)

INCLUDE("${CMAKE_SOURCE_DIR}/Tracing/Preprocessing/CurveletFiltering/fdct_wrapping_cpp/src/FindFFTW215.cmake")

SET( DEPENDS fftw-2.1.5 )

ADD_LIBRARY(CurveletClass
		Curvelet.h	Curvelet.cpp
		)
		
		QT4_WRAP_CPP( GUI_MOC_SRCS CurveletGUI.h)
ADD_LIBRARY(CurveletClassGUI
		CurveletGUI.h	CurveletGUI.cpp
  ${GUI_MOC_SRCS} 
		)
IF(FFTW_FOUND)
	INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
	ADD_LIBRARY(fdct_wrapping fdct_wrapping.cpp ifdct_wrapping.cpp fdct_wrapping_param.cpp)
	TARGET_LINK_LIBRARIES(fdct_wrapping ${FFTW_LIBRARIES})
	IF(WIN32)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
	ELSE(WIN32)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	ENDIF(WIN32)
	ADD_EXECUTABLE(curvelets curvelet_preprocessing.cpp)
	TARGET_LINK_LIBRARIES(curvelets CurveletClass fdct_wrapping ${ITK_LIBRARIES} )
	INSTALL(TARGETS curvelets RUNTIME DESTINATION bin)
	
	ADD_EXECUTABLE(CurveletGUI CurveletGUIDriver.cpp)
	TARGET_LINK_LIBRARIES(CurveletGUI CurveletClassGUI CurveletClass fdct_wrapping ${QT_LIBRARIES} ${ITK_LIBRARIES} )
	INSTALL(TARGETS CurveletGUI RUNTIME DESTINATION bin)
	
	INSTALL(TARGETS fdct_wrapping ARCHIVE DESTINATION lib)
ENDIF(FFTW_FOUND)
