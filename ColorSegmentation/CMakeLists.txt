PROJECT( colorbin )
cmake_minimum_required(VERSION 2.6)

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

TARGET_LINK_LIBRARIES( colorbin       ${ITK_LIBRARIES} )
