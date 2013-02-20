SET( SRCS
	GrAnisDiff.h
	GrAnisDiff.cpp
)

ADD_LIBRARY( preprocessing ${SRCS} )
TARGET_LINK_LIBRARIES ( preprocessing        ${ITK_LIBRARIES} )
#INSTALL( TARGETS preprocessing DESTINATION . )
