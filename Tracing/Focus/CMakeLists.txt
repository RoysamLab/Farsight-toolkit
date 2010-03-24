PROJECT( focus )
cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

SET( NSRCS	
    Focus.h
    Focus.cpp
   
    main.cpp
)

ADD_Executable( focus ${NSRCS} )

TARGET_LINK_LIBRARIES( focus ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )
