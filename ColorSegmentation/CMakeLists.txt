PROJECT( colorbin )
cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

SET( NSRCS	
    dhColors.h
    dhColors.cpp

	dhHistogram.h
    dhHistogram.cpp
    
    dhEvalState.h
    dhEvalState.cpp
    
    dhClassifiers.h
    
    dhSeedGrid.h
    dhSeedGrid.cpp
    
    ColorSegmentation.h
    ColorSegmentation.cpp
   
    main.cpp
)

ADD_Executable( colorbin ${NSRCS} )

TARGET_LINK_LIBRARIES( colorbin ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )
