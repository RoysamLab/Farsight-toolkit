PROJECT(TraceEdit)

SET( Libraries
    tiff
    ITKCommon
    ITKNumerics
    ITKStatistics
    ITKIO
    ITKBasicFilters
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
    ${QT_LIBRARIES} )

QT4_WRAP_CPP( GUI_MOC_SRCS View3D.h MergeModel.h TraceModel.h CellTraceModel.h)

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
  ${GUI_MOC_SRCS} 
  )

TARGET_LINK_LIBRARIES(Trace ${Libraries})

ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE} vtkTraceViewerDriver.cxx )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} Trace )

ADD_EXECUTABLE(merge merge_main.cxx)
TARGET_LINK_LIBRARIES(merge Trace )

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif(BUILD_TESTING)

SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME};Trace Editor")

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
endif(APPLE)

