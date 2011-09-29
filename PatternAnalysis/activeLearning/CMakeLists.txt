CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( MCLR)


################################################################################
# THESE ARE ITK, VTK, VXL, QT, and CUDA
# WE MAKE SURE THEY ARE PRESENT HERE
################################################################################

FIND_PACKAGE(VTK)
FIND_PACKAGE(VXL)
IF(VTK_FOUND AND VXL_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ELSEIF(VTK_FOUND AND VXL_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
ENDIF(VTK_FOUND AND VXL_FOUND)



SET( MCLR_SRCS
  mclr.cpp
  mclr_SM.cpp
)

SET( MCLR_HDRS
  mclr.h
  mclr_SM.h
)


ADD_LIBRARY(MCLR ${MCLR_SRCS} ${MCLR_HDRS})
#ADD_EXECUTABLE(mclr_main mclr_main.cpp) 
TARGET_LINK_LIBRARIES(MCLR vtkCommon vtkIO vtkCommon vtkHybrid vnl vnl_algo vnl_io mbl)

#TARGET_LINK_LIBRARIES(mclr_main MCLR vtkCommon vtkIO vtkCommon vtkHybrid vnl vnl_algo vnl_io mbl)

