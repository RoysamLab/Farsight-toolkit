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

QT4_WRAP_CPP( GUI_MOC_SRCS View3D.h MergeModel.h)

ADD_LIBRARY(Trace
  TraceBit.h    TraceBit.cxx
  TraceLine.h   TraceLine.cxx
  TraceGap.h    TraceGap.cxx
  TraceObject.h TraceObject.cxx
  MergeModel.h  MergeModel.cxx
  UndoBuffer.h UndoBuffer.cxx
  undobuff.h undobuff.cxx
)
ADD_EXECUTABLE( trace_editor  vtkTraceViewerDriver.cxx ${Sources} ${GUI_MOC_SRCS} )
TARGET_LINK_LIBRARIES( trace_editor ${Libraries} Trace TinyXML ftkGUI)
#ADD_EXECUTABLE( seed_editor vtkSeedViewerDriver.cxx ${Sources} )
#TARGET_LINK_LIBRARIES( seed_editor ${Libraries} Trace)

INSTALL(TARGETS trace_editor DESTINATION bin)
#INSTALL(TARGETS seed_editor DESTINATION bin)
