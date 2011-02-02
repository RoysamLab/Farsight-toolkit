SET( NSRCS
	main.cpp
)


ADD_EXECUTABLE( SubtractChannel ${NSRCS} )
TARGET_LINK_LIBRARIES( SubtractChannel ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )