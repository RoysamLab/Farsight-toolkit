cmake_minimum_required(VERSION 2.8)

PROJECT(ActiveValidation)

################################################################################
# THESE ARE ITK, VTK, VXL, QT, and 
# WE MAKE SURE THEY ARE PRESENT HERE
################################################################################
FIND_PACKAGE(VTK)
FIND_PACKAGE(VXL)

SET(ActiveValidation_SRCS
  ActiveValidation.cpp
  ActiveValidationWizard.cpp
)

SET(ActiveValidation_HDRS
  ActiveValidation.h
  ActiveValidationWizard.h
)
QT4_WRAP_CPP( ActiveValidation_MOC_REQD_HDRS ActiveValidationWizard.h)

ADD_LIBRARY(libActiveValidation ${ActiveValidation_SRCS} ${ActiveValidation_HDRS} ${ActiveValidation_MOC_REQD_HDRS})
ADD_EXECUTABLE(ActiveValidation Active_Validation_main.cpp )
TARGET_LINK_LIBRARIES(libActiveValidation 
  ${VTK_LIBRARIES}
  vtkCommon vtkIO vtkCommon vtkHybrid vnl vnl_algo vnl_io mbl ClusClus_lib)
  
TARGET_LINK_LIBRARIES( ActiveValidation libActiveValidation)

