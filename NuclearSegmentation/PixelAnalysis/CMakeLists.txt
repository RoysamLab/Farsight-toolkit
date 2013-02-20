SET( NSRCS
	ftkPixelLevelAnalysis.cpp
)

SET( NHDRS
	ftkPixelLevelAnalysis.h
)

ADD_LIBRARY( PixelLevelAnalysis ${NSRCS} ${NHDRS} )
TARGET_LINK_LIBRARIES( PixelLevelAnalysis ftkCommon Nuclear_Association          ${ITK_LIBRARIES} )