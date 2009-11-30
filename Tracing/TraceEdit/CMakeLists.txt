PROJECT(TraceEdit)
SET( Sources
     View3D.h View3D.cxx 
)

#Define (additional) ITK and VTK libraries that are needed for vtkTraceViewer to be built
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
)

QT4_WRAP_CPP( GUI_MOC_SRCS View3D.h MergeModel.h TraceModel.h)

ADD_LIBRARY(Trace
  TraceBit.h    TraceBit.cxx
  TraceLine.h   TraceLine.cxx
  TraceGap.h    TraceGap.cxx
  TraceObject.h TraceObject.cxx
  TraceModel.h  TraceModel.cxx
  MergeModel.h  MergeModel.cxx
  UndoBuffer.h	UndoBuffer.cxx
  undobuff.h	undobuff.cxx
  branchPT.h	BranchPT.cxx
  vtkPlotEdges.h	vtkPlotEdges.cxx
)
ADD_EXECUTABLE( trace_editor  vtkTraceViewerDriver.cxx ${Sources} ${GUI_MOC_SRCS} )
TARGET_LINK_LIBRARIES( trace_editor ${Libraries} Trace TinyXML ftkGUI)

INSTALL(TARGETS trace_editor DESTINATION bin/Tracing)
SET(CPACK_PACKAGE_EXECUTABLES "trace_editor;Trace Editor")

