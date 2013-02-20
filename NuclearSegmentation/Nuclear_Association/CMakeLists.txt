SET( SRCS
	ftkNuclearAssociationRules.h
	ftkNuclearAssociationRules.cpp
	AssociationAuxFns.cpp
	VolumeOfInterest.cxx
	VolumeOfInterest.h
)

ADD_LIBRARY( Nuclear_Association ${SRCS} )
TARGET_LINK_LIBRARIES ( Nuclear_Association CytoplasmSegmentation ftkCommon
  ftkFeatures ${ITK_LIBRARIES} ${VTK_LIBRARIES} vnl vnl_algo )
#INSTALL( TARGETS Nuclear_Association DESTINATION . )
