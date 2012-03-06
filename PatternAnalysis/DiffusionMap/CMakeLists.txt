CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( DIFFUSIONMAP )


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



SET( DIFFUSION_SRCS
  DiffusionMap.cpp
)

SET( DIFFUSION_HDRS
  DiffusionMap.h
)


ADD_LIBRARY(diffusion_map ${DIFFUSION_SRCS} ${DIFFUSION_HDRS})

TARGET_LINK_LIBRARIES(diffusion_map ${VTK_LIBRARIES} ftkGraphs vnl vnl_algo vnl_io mbl)



