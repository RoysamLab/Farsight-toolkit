SET( NSRCS
	Voting_Seg.cpp
)

SET( NHDRS
	GeneralIO.h
	Voting_Seg.h
)

ADD_LIBRARY( Voting_Seg ${NSRCS} ${NHDRS} )
TARGET_LINK_LIBRARIES( Voting_Seg LIBftkVoting Yousef_Nucleus_Seg ftkCommon ${ITK_LIBRARIES} )