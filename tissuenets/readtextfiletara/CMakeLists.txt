CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( readtextfiletara )

ADD_EXECUTABLE(readtextfiletara readtextfiletara.cpp)
TARGET_LINK_LIBRARIES(readtextfiletara)
#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

################################################################################
################################################################################
################################################################################
################################################################################
# MACROS (functions)
################################################################################
#MACRO(MACRO_INCLUDE_COMMON_LIBS)


link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )

