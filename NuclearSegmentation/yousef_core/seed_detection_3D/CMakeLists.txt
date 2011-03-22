SET( SRCS
	seedsDetection_3D.cxx	seedsdetection_3D.h
 itklaplacianrecursivegaussianimagefilternew.txx
	itklaplacianrecursivegaussianimagefilternew.h )


#########################
If CUDA is found, then add CUDA file and add to library
###########################

if (CUDA_FOUND)
    SET (SRCS ${SRCS} LocalMaxima3D.cu )
    CUDA_ADD_LIBRARY( seedsDetection_3D ${SRCS} )

else (CUDA_FOUND)
    ADD_LIBRARY( seedsDetection_3D ${SRCS} )

endif (CUDA_FOUND)

TARGET_LINK_LIBRARIES ( seedsDetection_3D ITKBasicFilters ITKCommon ITKIO)

#INSTALL( TARGETS seedsDetection_3D DESTINATION . )

