# This project is designed to be built outside the Insight source tree.
cmake_minimum_required(VERSION 2.6)

#PROJECT(TracingSystem)

#find_package(Qt4 REQUIRED) # find and setup Qt4 for this project
#include(${QT_USE_FILE})

SET( TracingSysLibraries
    ${ITK_LIBRARIES}
    rtvl
    rgtl
    ${VTK_LIBRARIES}
    QVTK
)

SET( GUI_RCS
  icons.qrc
)

#Check if ITK is found
#FIND_PACKAGE(ITK)
#IF(ITK_FOUND)
#    INCLUDE(${ITK_USE_FILE})
#ELSE(ITK_FOUND)
#    MESSAGE(FATAL_ERROR
#    "Cannot build TrackerGUI without ITK.  Please set ITK_DIR.")
#ENDIF(ITK_FOUND)


#FIND_PACKAGE(VTK)
#FIND_PACKAGE(VXL)
#IF(VTK_FOUND AND VXL_FOUND)
#    INCLUDE( ${VTK_USE_FILE} )
#    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
#    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
#ELSEIF(VTK_FOUND AND VXL_FOUND)
#    MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
#ENDIF(VTK_FOUND AND VXL_FOUND)

QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )
QT4_WRAP_CPP(GUI_MOC_SRCS TracingGUI.h TracingView.h ObjectSelection.h ParametersGroup.h SlidersGroup.h OpenSnakeTracer.h dialogs.h dialogs_montage.h ScatterView.h dialogs_scatter.h)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )
#INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_INCLUDE})
#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tinyxml )
#ADD_SUBDIRECTORY( tinyxml )

SET(TracingCore
	TracingCore/PointOperation.h
    TracingCore/PointOperation.cpp
    TracingCore/ImageOperation.h
    TracingCore/ImageOperation.cpp
    TracingCore/SnakeOperation.h
    TracingCore/SnakeOperation.cpp
    TracingCore/Montage/Register.h
    TracingCore/Montage/Register.cpp
    )

ADD_LIBRARY( TraceCoreLibrary ${TracingCore} )

SET(TracingSystem_SOURCES
  main.cpp
  SlidersGroup.cpp
  ParametersGroup.cpp
  dialogs.cpp
  dialogs_montage.cpp
  dialogs_scatter.cpp
  TracingGUI.cpp
  TracingView.cpp
  ObjectSelection.cpp
  OpenSnakeTracer.cpp
  ScatterView.cpp
  Rendering/ImageActors.h
  Rendering/ImageActors.cxx
  ${GUI_MOC_SRCS}
  ${GUI_RC_SRCS}
  )

SET(TracingSystem_LIBRARIES
  ${QT_LIBRARIES} TraceCoreLibrary ftkCommon ${TracingSysLibraries} TinyXML)

add_farsight_executable(TracingSystem "${TracingSystem_SOURCES}" "${TracingSystem_LIBRARIES}")
