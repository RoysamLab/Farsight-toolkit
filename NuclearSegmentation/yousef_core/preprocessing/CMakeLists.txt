SET( SRCS
	GrAnisDiff.h
	GrAnisDiff.cpp
)

ADD_LIBRARY( preprocessing ${SRCS} )
TARGET_LINK_LIBRARIES ( preprocessing ITKBasicFilters ITKCommon ITKIO)
#INSTALL( TARGETS preprocessing DESTINATION . )
