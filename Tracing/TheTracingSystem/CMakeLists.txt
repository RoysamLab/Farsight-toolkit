# This project is designed to be built outside the Insight source tree.
cmake_minimum_required(VERSION 2.6)

#PROJECT(TracingSystem)

#find_package(Qt4 REQUIRED) # find and setup Qt4 for this project
#include(${QT_USE_FILE})

SET( TracingSysLibraries
    ${ITK_LIBRARIES}
    rtvl
    rgtl
    vtkRendering
    vtkIO
    vtkHybrid
    vtkWidgets
    vtkCommon
    vtksys
    vtkFiltering
    vtkVolumeRendering
    QVTK
)

SET( GUI_RCS
  icons.qrc
)

#Check if ITK is found
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
    "Cannot build TrackerGUI without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


FIND_PACKAGE(VTK)
FIND_PACKAGE(VXL)
IF(VTK_FOUND AND VXL_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ELSEIF(VTK_FOUND AND VXL_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
ENDIF(VTK_FOUND AND VXL_FOUND)

QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )
QT4_WRAP_CPP(GUI_MOC_SRCS TracingGUI.h TracingView.h ObjectSelection.h ParametersGroup.h SlidersGroup.h OpenSnakeTracer.h dialogs.h dialogs_montage.h ScatterView.h dialogs_scatter.h)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )
#INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_INCLUDE})
#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tinyxml )
#SUBDIRS( tinyxml )

ADD_EXECUTABLE(TracingSystem 
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
      Rendering/itkImageToVTKImageFilter.h
      Rendering/itkImageToVTKImageFilter.txx

      TracingCore/PointOperation.h
      TracingCore/PointOperation.cpp
      TracingCore/ImageOperation.h
      TracingCore/ImageOperation.cpp
      TracingCore/SnakeOperation.h
      TracingCore/SnakeOperation.cpp

      #TracingCore/Filters/itkLabelGeometryImageFilter.h
      #TracingCore/Filters/itkLabelGeometryImageFilter.txx
      #TracingCore/Filters/itkBinaryThinningImageFilter3D.h
      #TracingCore/Filters/itkBinaryThinningImageFilter3D.txx
      #TracingCore/Filters/itkMaximumProjectionImageFilter.h
      #TracingCore/Filters/itkProjectionImageFilter.h
      #TracingCore/Filters/itkProjectionImageFilter.txx
 
      TracingCore/Montage/Register.h
      TracingCore/Montage/Register.cpp

#TracingCore/Vesselness/itkHessianSmoothed3DToVesselnessMeasureImageFilter.h
#      TracingCore/Vesselness/itkHessianSmoothed3DToVesselnessMeasureImageFilter.txx

#TracingCore/Vesselness/itkMultiScaleHessianBasedMeasureImageFilter.h
#      TracingCore/Vesselness/itkMultiScaleHessianBasedMeasureImageFilter.txx

#TracingCore/Vesselness/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h
#      TracingCore/Vesselness/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.txx

     ${GUI_MOC_SRCS}
     ${GUI_RC_SRCS})


TARGET_LINK_LIBRARIES(TracingSystem ${QT_LIBRARIES} ${TracingSysLibraries} TinyXML)
 
#MAKE A GUI EXECUTABLE, NOT A CONSOLE APPLICATION
#SET_TARGET_PROPERTIES( TracingSystem PROPERTIES WIN32_EXECUTABLE 1)
INSTALL( TARGETS TracingSystem DESTINATION bin )

#####################################################################################
# CPACK SETUP FOR BUILDING INSTALLERS
######################################################################################
SET(CPACK_PACKAGE_VERSION_MAJOR 0)
SET(CPACK_PACKAGE_VERSION_MINOR 2)
SET(CPACK_PACKAGE_VERSION_PATCH 1)
SET(CPACK_PACKAGE_VENDOR "Rensselaer Polytechnic Institute")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${TracingSystem_SOURCE_DIR}/farsight-description.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolkit of image analysis methods")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "TracingSystem ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_MAJOR="${CPACK_PACKAGE_VERSION_MAJOR}")
ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_MINOR="${CPACK_PACKAGE_VERSION_MINOR}")
ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_PATCH="${CPACK_PACKAGE_VERSION_PATCH}")

#first screen.  make longer
SET(CPACK_RESOURCE_FILE_WELCOME "${TracingSystem_SOURCE_DIR}/welcome.txt")
#second screen?
SET(CPACK_RESOURCE_FILE_README "${TracingSystem_SOURCE_DIR}/readme.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${TracingSystem_SOURCE_DIR}/license.txt") 

SET(CPACK_NSIS_MENU_LINKS "http://www.farsight-toolkit.org" "Farsight Online Help")
SET(CPACK_PACKAGE_EXECUTABLES "TracingSystem;Tracing System")  

SET(CPACK_PACKAGING_INSTALL_PREFIX "/")

INCLUDE(CPack)

#Qt libraries 
IF(WIN32)
  INSTALL(FILES ${QT_LIBRARY_DIR}/QtCore4.dll ${QT_LIBRARY_DIR}/QtGui4.dll DESTINATION bin)
ENDIF(WIN32)
IF(APPLE)
  INSTALL(DIRECTORY ${QT_QTCORE_LIBRARY_RELEASE} ${QT_QTGUI_LIBRARY_RELEASE} ${QT_QTSQL_LIBRARY_RELEASE} ${QT_QTNETWORK_LIBRARY_RELEASE} DESTINATION bin)
ENDIF(APPLE)

INCLUDE(InstallRequiredSystemLibraries)
