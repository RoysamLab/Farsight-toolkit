CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( read_cell_types )

include_directories( ${Boost_INCLUDE_DIRS} C:/tinyxml)

ADD_EXECUTABLE(read_cell_types read_cell_types.cpp read_cell_types.h)
TARGET_LINK_LIBRARIES( read_cell_types TinyXML)

################################################################################