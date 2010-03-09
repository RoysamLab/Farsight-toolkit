SET( NSRCS
	extrema.cpp
	dhArray3D.h
	
    dhColors.h
    dhColors.cpp
    dhImage.h
    dhImage.cpp
    Color_Image.h
    Color_Image.cpp
    PPM_Image.h
    PPM_Image.cpp
	
	Col_Slice.h
	RGB_Slice.h
	RLI_Slice.h
	
	RGB_Atype.h
	RGB_Atype.cpp
	
	dhHistogram.h
    dhHistogram.cpp
    
    dhClassifiers.h
    
    dhEvalState.h
    dhEvalState.cpp
    
    ColorSegmentation.h
    ColorSegmentation.cpp
   
    main.cpp
)

ADD_Executable( colorbin ${NSRCS} )

TARGET_LINK_LIBRARIES( colorbin ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )
