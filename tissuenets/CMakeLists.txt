CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#PROJECT( bionet )

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )

ADD_EXECUTABLE(bionet bionet.cpp)
TARGET_LINK_LIBRARIES( bionet tiff vtkCommon vtkIO vtkInfovis vtkViews vtkRendering TinyXML)

INSTALL( TARGETS bionet DESTINATION bin ) 
INSTALL( FILES tissuenets_demo.py DESTINATION python )
INSTALL( FILES saxon9.jar DESTINATION bin )


SUBDIRS( histogram )
SUBDIRS( tissuenets )

#  FIND_PACKAGE(VTK)
#  IF(VTK_FOUND)
#    INCLUDE( ${VTK_USE_FILE} )
#  ELSEIF(VTK_FOUND)
#    MESSAGE(FATAL_ERROR "Cannot build without VTK.")
#  ENDIF(VTK_FOUND)
