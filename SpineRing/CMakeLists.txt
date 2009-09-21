# CMakeLists.txt for spinereing/detrection
# Author: Hussein Sharafeddin

cmake_minimum_required( VERSION 2.4)
PROJECT ( SpineSegmentation )

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
             INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
             MESSAGE(FATAL_ERROR
             "Cannot build SETracing without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
             INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
             MESSAGE(FATAL_ERROR "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)


IF(ITK_USE_SYSTEM_VXL)
  FIND_PACKAGE(VXL)
  IF(VXL_FOUND)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
    SET(VXL_NUMERICS_LIBRARIES vnl vnl_algo netlib vcl)
  ELSE(VXL_FOUND)
    MESSAGE(FATAL_ERROR "Must set VXL_DIR for ITK_USE_SYSTEM_VXL.")
  ENDIF(VXL_FOUND)
ELSE(ITK_USE_SYSTEM_VXL)
  SET(VXL_NETLIB_INCLUDE_DIR ${ITK_SOURCE_DIR}/Utilities/vxl/v3p/netlib)
  SET(VXL_NUMERICS_LIBRARIES itkvnl itkvnl_algo itkvcl)
ENDIF(ITK_USE_SYSTEM_VXL)


SET (XML_sources
	tinystr.cpp
	tinyxml.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp
	tinystr.h
	tinyxml.h
)	



SET ( SpineSegmentation_sources
      RingSE.cxx
      SpineRing.cxx
      TraceSegNode.cxx
      ImageProc.cxx
      SpinePaths.cxx
      FMPath.cxx
      SpeedFunc.cxx
)

SET ( SpineSegmentation_headers 
      itkSingleImageCostFunction.h
      itkArrivalFunctionToPathFilter.h
      itkSpeedFunctionToPathFilter.h
      SpineRing.h
      SpineConsts.h
      SpineUtils.h
      TraceSegNode.h
      CommonTypeDefs.h
      SpineTypedefs.h
      ImageProc.h
)      


INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${SOURCE_PATH}
  ${VXL_INCLUDE_DIRS}
 )

ADD_EXECUTABLE ( SpineSegmentation ${SpineSegmentation_sources} ${SpineSegmentation_headers}  ${XML_sources} SpineMain.cxx )
TARGET_LINK_LIBRARIES( SpineSegmentation ITKCommon ITKIO ITKNumerics)


