cmake_minimum_required (VERSION 2.6)
project (neighborfile)
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

SET (neighborfile_SRCS window.cxx neighborfile.cxx) 
SET (neighborfile_HDRS window.h)
SET (neighborfile_MOC_HDRS window.h)

QT4_WRAP_CPP( neighborfile_MOC_SRCS ${neighborfile_MOC_HDRS} )

add_executable(neighborfile ${neighborfile_SRCS} ${neighborfile_HDRS} ${neighborfile_MOC_SRCS})
TARGET_LINK_LIBRARIES(neighborfile ${QT_LIBRARIES})