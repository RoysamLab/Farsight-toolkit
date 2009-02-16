SET( SRCS
	seedsDetection_3D.cxx
	seedsdetection_3D.h
	itklaplacianrecursivegaussianimagefilternew.txx
	itklaplacianrecursivegaussianimagefilternew.h
)

ADD_LIBRARY( seedsDetection_3D ${SRCS} )
TARGET_LINK_LIBRARIES ( seedsDetection_3D ITKBasicFilters ITKCommon ITKIO)
#INSTALL( TARGETS seedsDetection_3D DESTINATION . )

