SET( NSRCS
	main.cpp
)


ADD_EXECUTABLE( SubtractChannel ${NSRCS} )
TARGET_LINK_LIBRARIES( SubtractChannel       ${ITK_LIBRARIES} )