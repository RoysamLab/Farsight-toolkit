cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE( ${ITK_USE_FILE} )
ELSEIF(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE( ${VTK_USE_FILE} )
ELSEIF(VTK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without VTK.")
ENDIF(VTK_FOUND)

SET(SRCS
  vtkPlotEdges.h                vtkPlotEdges.cxx 
  vtkConvertMDLTracesToTrees.h  vtkConvertMDLTracesToTrees.cxx 
  main.cxx)

ADD_EXECUTABLE(connectTraces ${SRCS})
TARGET_LINK_LIBRARIES(connectTraces
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  )

#ADD_EXECUTABLE(test test.cxx vtkPlotEdges.cxx)
#TARGET_LINK_LIBRARIES(test vtkCommon vtkIO vtkGraphics)
