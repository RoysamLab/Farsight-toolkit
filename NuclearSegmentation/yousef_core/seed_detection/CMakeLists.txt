SET( SRCS
	seedsDetection_2D_Cpp.cxx
	seedsdetection.h
)

ADD_LIBRARY( seedsDetection_2D ${SRCS} )


TARGET_LINK_LIBRARIES ( seedsDetection_2D        ${ITK_LIBRARIES} )
#INSTALL( TARGETS seedsDetection_2D DESTINATION . )

