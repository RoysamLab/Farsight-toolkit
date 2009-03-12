SET( Sources
     View3D.h View3D.cxx
    vtkTraceViewerDriver.cxx Trace.h Trace.cpp 

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

ADD_EXECUTABLE( vtkTraceViewer  ${Sources})
TARGET_LINK_LIBRARIES( vtkTraceViewer ${Libraries})

