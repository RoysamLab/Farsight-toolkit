SET( SRCS
	ftkNuclearAssociationRules.h
	ftkNuclearAssociationRules.cpp
	AssociationAuxFns.cpp
)

ADD_LIBRARY( Nuclear_Association ${SRCS} )
TARGET_LINK_LIBRARIES ( Nuclear_Association CytoplasmSegmentation ftkCommon
  ftkFeatures ${ITK_LIBRARIES})
#INSTALL( TARGETS Nuclear_Association DESTINATION . )
