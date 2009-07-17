# the next line sets up include and link directories and defines some variables that we will use.
# you can modify the behavior by setting some variables, e.g.
#   set(QT_USE_OPENGL TRUE)
# -> this will cause cmake to include and link against the OpenGL module

find_package(Qt4 REQUIRED) # find and setup Qt4 for this project

include(${QT_USE_FILE})

# the variable "qtproject_SRCS" contains all .cpp files of this project
set(qtproject_SRCS
    #QT1.cpp
    Seed3D.h Seed3D.cxx
)


#Define (additional) ITK and VTK libraries that are needed for vtkTraceViewer to be builtSET( Libraries    ITKCommon    ITKNumerics    ITKStatistics    ITKIO    ITKBasicFilters    vtkRendering    vtkIO    vtkHybrid    vtkWidgets    vtkCommon    vtksys    vtkFiltering    vtkVolumeRendering    QVTK    QVTKWidgetPlugin   # TinyXML/home/raghav/macecutter)

#Check if ITK is found
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build view3d without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build view 3d without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)


# tell cmake to create .moc files for all files in the variable qtproject_SRCS that require such a file.
# note: this assumes that you use #include "header.moc" in your files
#qt4_automoc(${qtproject_SRCS})

QT4_WRAP_CPP( GUI_MOC_SRCS Seed3D.h )
#ADD_LIBRARY(Trace Trace.h Trace.cpp)ADD_EXECUTABLE( Seededitor  vtkSeeds3D.cxx ${Sources} ${GUI_MOC_SRCS} )

TARGET_LINK_LIBRARIES( Seededitor ${Libraries})

#INSTALL(TARGETS trace_editor DESTINATION bin)	


# create an executable file named "qtproject" from the source files in the variable "qtproject_SRCS".
#add_executable(qtproject ${qtproject_SRCS})

# link the "qtproject" target against the Qt libraries. which libraries exactly, is defined by the "include(${QT_USE_FILE})" line above, which sets up this variable.
#target_link_libraries(qtproject ${QT_LIBRARIES})


