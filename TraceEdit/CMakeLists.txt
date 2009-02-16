# set project's name
PROJECT( vtkTraceViewer )
cmake_minimum_required(VERSION 2.6)

#Check if ITK is found
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build vtkTraceViewer without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build vtkTraceViewer without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)

#  FIND_PACKAGE(Qt4)
#  IF(QT4_FOUND)
 #   INCLUDE( ${QT_USE_FILE} )
#    IF (QT_QTXML_FOUND)
#      INCLUDE_DIRECTORIES(${QT_QTXML_INCLUDE_DIR})
#    ELSE(QT_QTXML_FOUND)
#      MESSAGE(FATAL ERROR "Qt XML libraries missing.")
#    ENDIF(QT_QTXML_FOUND)
#  ELSE(QT4_FOUND)
#    MESSAGE(FATAL_ERROR "Cannot build without Qt4")
#  ENDIF(QT4_FOUND)
  

FIND_PACKAGE(LIBXML2)
IF(LIBXML2_FOUND)
      INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
ELSE(LIBXML2_FOUND)
     MESSAGE(FATAL_ERROR
      "&*&*&*& Cannot build without LIBXML2.  Please set LIBXML2_DIR.")
ENDIF(LIBXML2_FOUND)

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
 #   QVTK
 #   QVTKWidgetPlugin
)


ADD_EXECUTABLE( vtkTraceViewer  ${Sources})

TARGET_LINK_LIBRARIES( vtkTraceViewer ${LIBXML2_LIBRARIES} ${Libraries})

