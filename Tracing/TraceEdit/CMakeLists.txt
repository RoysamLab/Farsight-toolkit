cmake_minimum_required(VERSION 2.8)
PROJECT(TraceEdit)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

SET( Libraries
    tiff
    ${ITK_LIBRARIES}
    vtkRendering
    vtkIO
    vtkHybrid
    vtkWidgets
    vtkCommon
    vtksys
    vtkFiltering
    vtkVolumeRendering
    QVTK
    ftkGUI
    ftkCommon
    ftkGraphs
    ${QT_LIBRARIES} )

QT4_WRAP_CPP( GUI_MOC_SRCS View3D.h MergeModel.h TraceModel.h CellTraceModel.h cellexport.h screenshot.h)

ADD_LIBRARY(Trace
  TraceBit.h        TraceBit.cxx
  TraceLine.h       TraceLine.cxx
  TraceGap.h        TraceGap.cxx
  TraceObject.h     TraceObject.cxx
  TraceModel.h      TraceModel.cxx
  MergeModel.h      MergeModel.cxx
  branchPT.h        BranchPT.cxx
  ImageActors.h     ImageActors.cxx
  CellTrace.h       CellTrace.cxx
  CellTraceModel.h  CellTraceModel.cxx
  vtkPlotEdges.h    vtkPlotEdges.cxx 
  View3D.h          View3D.cxx 
  cellexport.h	    cellexport.cxx
  screenshot.h	    screenshot.cxx
  ${GUI_MOC_SRCS} 
  )
OPTION( USE_GPUREN "Use GPU for Rendering" OFF )
IF(USE_GPUREN)
	ADD_DEFINITIONS(-DUSE_GPUREN)
ENDIF(USE_GPUREN)

TARGET_LINK_LIBRARIES(Trace ${Libraries})

ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE} vtkTraceViewerDriver.cxx )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} Trace )

ADD_EXECUTABLE(merge merge_main.cxx)
TARGET_LINK_LIBRARIES(merge Trace )

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif(BUILD_TESTING)

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
endif(APPLE)

