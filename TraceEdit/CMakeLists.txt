PROJECT(TraceEdit)
SET( Sources
     View3D.h View3D.cxx
    vtkTraceViewerDriver.cxx

)

#Define (additional) ITK and VTK libraries that are needed for vtkTraceViewer to be built
SET( Libraries
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
    QVTKWidgetPlugin
    TinyXML
)

ADD_LIBRARY(Trace Trace.h Trace.cpp)
ADD_EXECUTABLE(trace_editor  ${Sources})
TARGET_LINK_LIBRARIES( trace_editor ${Libraries} Trace TinyXML)

INSTALL(TARGETS trace_editor DESTINATION bin)
