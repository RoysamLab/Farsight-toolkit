CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#PROJECT( tissuenets )

find_package( Boost 1.37 COMPONENTS regex graph )

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS})

ADD_EXECUTABLE(tissuenets tissuenets.cpp)
TARGET_LINK_LIBRARIES( tissuenets vtkCommon vtkIO vtkInfovis vtkViews vtkRendering TinyXML)

INSTALL( TARGETS tissuenets DESTINATION bin )