cmake_minimum_required(VERSION 2.8)
PROJECT(TraceEdit)

QT4_WRAP_CPP( GUI_MOC_SRCS TraceView3D.h MergeModel.h TraceModel.h CellTraceModel.h cellexport.h screenshot.h)

SET( TRACE_SRCS
  TraceBit.cxx
  TraceLine.cxx
  TraceGap.cxx
  TraceObject.cxx
  TraceModel.cxx
  MergeModel.cxx
  BranchPT.cxx
  ImageActors.cxx
  CellTrace.cxx
  CellTraceModel.cxx
  vtkPlotEdges.cxx
  TraceView3D.cxx
  cellexport.cxx
  screenshot.cxx
  GridActors.cxx
  VolumeOfInterest.cxx
  ${GUI_MOC_SRCS} )

SET( TRACE_LIBS
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
    MCLR
    ${QT_LIBRARIES} )

IF (BUILD_SPD)
  LIST( APPEND TRACE_LIBS SPD )
ENDIF (BUILD_SPD)

IF(BUILD_CLUSCLUS)
  LIST( APPEND TRACE_LIBS ClusClus_lib )
ENDIF (BUILD_CLUSCLUS)

IF(QtTestingFound)
  LINK_DIRECTORIES( ${QtTesting_LIBRARY_DIR} )
  ADD_DEFINITIONS(-DUSE_QT_TESTING)
ENDIF()

OPTION( USE_GPUREN "Use GPU for Rendering" OFF )
IF(USE_GPUREN)
	ADD_DEFINITIONS(-DUSE_GPUREN)
ENDIF(USE_GPUREN)

ADD_LIBRARY( Trace STATIC ${TRACE_SRCS} )
TARGET_LINK_LIBRARIES( Trace ${TRACE_LIBS} )

ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE} vtkTraceViewerDriver.cxx )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} Trace)

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif(BUILD_TESTING)

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
endif(APPLE)

