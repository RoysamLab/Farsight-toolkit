SET( NSRCS
	CytoplasmSegmentation.cpp
	whole_cell.cpp
	gradientweighteddistancemap.cpp
)

SET( NHDRS
	CytoplasmSegmentation.h
	whole_cell.h
	itkMorphologicalWatershedFromMarkersImageFilter.h
	itkMorphologicalWatershedFromMarkersImageFilter.txx
)

ADD_LIBRARY( CytoplasmSegmentation ${NSRCS} ${NHDRS} )
TARGET_LINK_LIBRARIES( CytoplasmSegmentation Yousef_Nucleus_Seg ftkCommon          ${ITK_LIBRARIES} )