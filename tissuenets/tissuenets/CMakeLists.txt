CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( tissuenets )

find_package( Boost 1.37 COMPONENTS regex graph )

  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
  ELSEIF(VTK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK.")
  ENDIF(VTK_FOUND)


link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} C:/tinyxml)

ADD_EXECUTABLE(tissuenets tissuenets.cpp tissuenets.h vtkSurfaceRepresentation_.cxx vtkSurfaceRepresentation_.h vtkConvertSelection.cxx vtkConvertSelection.h vtkGlyph3D_.cxx vtkGlyph3D_.h )
TARGET_LINK_LIBRARIES( tissuenets vtkCommon vtkIO vtkInfovis vtkViews vtkRendering TinyXML)

