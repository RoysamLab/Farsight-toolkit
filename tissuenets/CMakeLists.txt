CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( bionet )

find_package( Boost 1.37 COMPONENTS regex graph )

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS}  c:/tinyXML)

ADD_EXECUTABLE(bionet bionet.cpp)
TARGET_LINK_LIBRARIES( bionet vtkCommon vtkIO vtkInfovis vtkViews vtkRendering TinyXML)

  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
  ELSEIF(VTK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK.")
  ENDIF(VTK_FOUND)
