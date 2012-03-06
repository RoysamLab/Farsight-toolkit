CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET( FTKGRAPH_SRCS
  FTKgraph.cpp
  kNearestObjects.txx
)

SET( FTKGRAPH_HDRS
  FTKgraph.h
  kNearestObjects.h
)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_LIBRARY(ftkGraphs ${FTKGRAPH_SRCS} ${FTKGRAPH_HDRS})

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK")
ENDIF(VTK_FOUND)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})


TARGET_LINK_LIBRARIES( ftkGraphs vtkCommon vtkInfovis vtkViews vtkRendering vtkFiltering vtkGraphics   ftkGUI ftkCommon ftkFeatures     ${ITK_LIBRARIES} )
