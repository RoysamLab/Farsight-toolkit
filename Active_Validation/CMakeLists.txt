cmake_minimum_required(VERSION 2.8)

PROJECT(ActiveValidation)

################################################################################
# THESE ARE ITK, VTK, VXL, QT, and 
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

SET(ActiveValidation_SRCS
  ActiveValidation.cpp
  ActiveValidationWizard.cpp
)

SET(ActiveValidation_HDRS
  ActiveValidation.h
  ActiveValidationWizard.h
)
QT4_WRAP_CPP( ActiveValidation_MOC_REQD_HDRS ActiveValidationWizard.h)

#ADD_LIBRARY(libActiveValidation ${ActiveValidation_SRCS} ${ActiveValidation_HDRS})
ADD_EXECUTABLE(ActiveValidation Active_Validation_main.cpp ${ActiveValidation_SRCS} ${ActiveValidation_HDRS} ${ActiveValidation_MOC_REQD_HDRS})
TARGET_LINK_LIBRARIES(ActiveValidation #libActiveValidation
  ${VTK_LIBRARIES}
  vtkCommon vtkIO vtkCommon vtkHybrid vnl vnl_algo vnl_io mbl ClusClus_lib)

