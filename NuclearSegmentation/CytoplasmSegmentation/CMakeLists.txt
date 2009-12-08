SET( NSRCS
	CytoplasmSegmentation.cpp
	whole_cell.cpp
	gradientweighteddistancemap.cpp
)

SET( NHDRS
	CytoplasmSegmentation.h
	whole_cell.h
	itkMorphologicalWatershedFromMarkersImageFilter.h
	itkLabelGeometryImageFilter.h
	itkMorphologicalWatershedFromMarkersImageFilter.txx
	itkLabelGeometryImageFilter.txx
)

ADD_LIBRARY( CytoplasmSegmentation ${NSRCS} ${NHDRS} )
TARGET_LINK_LIBRARIES( CytoplasmSegmentation Yousef_Nucleus_Seg ITKCommon ITKBasicFilters ITKAlgorithms ITKIO )