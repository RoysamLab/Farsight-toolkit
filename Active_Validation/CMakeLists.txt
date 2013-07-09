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
  ${VTK_LIBRARIES} vnl vnl_algo vnl_io mbl ClusClus_lib)
  
TARGET_LINK_LIBRARIES( ActiveValidation libActiveValidation)

