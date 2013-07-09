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


#Define (additional) ITK and VTK libraries that are needed for vtkTraceViewer to be built
SET( Libraries
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
   # TinyXML/home/raghav/macecutter
)

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

ADD_EXECUTABLE( SeedEditor  vtkSeeds3D.cxx Seed3D.cxx ${Sources} ${GUI_MOC_SRCS} )

TARGET_LINK_LIBRARIES( SeedEditor ${Libraries})


#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
SET_TARGET_PROPERTIES( SeedEditor  PROPERTIES WIN32_EXECUTABLE 1)
INSTALL( TARGETS SeedEditor DESTINATION ${INSTALL_BIN_DIR} )

