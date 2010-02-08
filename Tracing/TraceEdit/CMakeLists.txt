PROJECT(TraceEdit)

ADD_LIBRARY(Trace
  TraceBit.h      TraceBit.cxx
  TraceLine.h     TraceLine.cxx
  TraceGap.h      TraceGap.cxx
  TraceObject.h   TraceObject.cxx
  TraceModel.h    TraceModel.cxx
  MergeModel.h    MergeModel.cxx
  UndoBuffer.h	  UndoBuffer.cxx
  undobuff.h	    undobuff.cxx
  branchPT.h	    BranchPT.cxx
  ImageActors.h		ImageActors.cxx
  vtkPlotEdges.h	vtkPlotEdges.cxx 
  TraceProjectManager.h TraceProjectManager.cxx)

SET( Sources View3D.h View3D.cxx vtkTraceViewerDriver.cxx )

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
    TinyXML
    ftkGUI
    Trace
    ${QT_LIBRARIES} )

QT4_WRAP_CPP( GUI_MOC_SRCS View3D.h MergeModel.h TraceModel.h)

ADD_EXECUTABLE(${PROJECT_NAME} ${EXE_TYPE}
   ${Sources} ${GUI_MOC_SRCS} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${Libraries} )

SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME};Trace Editor")

if(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS ${PROJECT_NAME} DESTINATION bin )
endif(APPLE)

