SET( SRCS
	ftkNuclearAssociationRules.h
	ftkNuclearAssociationRules.cpp
)

ADD_LIBRARY( Nuclear_Association ${SRCS} )
TARGET_LINK_LIBRARIES ( Nuclear_Association ftkCommon ITKBasicFilters ITKCommon ITKIO)
#INSTALL( TARGETS Nuclear_Association DESTINATION . )