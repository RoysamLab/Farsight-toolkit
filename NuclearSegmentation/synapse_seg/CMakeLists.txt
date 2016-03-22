SET( NSRCS
	Synapse_Seg.cpp
)

SET( NHDRS
	Synapse_Seg.h
	Synapse_Seg.txx
)

ADD_LIBRARY( Synapse_Seg ${NSRCS} ${NHDRS} )
TARGET_LINK_LIBRARIES( Synapse_Seg ftkCommon ${ITK_LIBRARIES} )