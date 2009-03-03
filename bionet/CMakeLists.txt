CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( bionet )

find_package( Boost 1.37 COMPONENTS regex graph )

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )

ADD_EXECUTABLE(bionet bionet.cpp)
TARGET_LINK_LIBRARIES( bionet vtkCommon vtkIO vtkInfovis vtkViews vtkRendering TinyXML)
