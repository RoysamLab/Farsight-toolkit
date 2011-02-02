PROJECT( focus )
cmake_minimum_required(VERSION 2.6)

SET( NSRCS	
    Focus.h
    Focus.cpp
   
    main.cpp
)

ADD_Executable( focus ${NSRCS} )
ADD_Executable( focus_grayscale Focus.h Focus.cpp main_grayscale.cpp)

TARGET_LINK_LIBRARIES( focus ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )
TARGET_LINK_LIBRARIES( focus_grayscale ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )

INSTALL(TARGETS focus RUNTIME DESTINATION bin)
INSTALL(TARGETS focus_grayscale RUNTIME DESTINATION bin)
