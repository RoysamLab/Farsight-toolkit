SET( SRCS
	seedsDetection_3D.cxx	seedsdetection_3D.h
 itklaplacianrecursivegaussianimagefilternew.txx
	itklaplacianrecursivegaussianimagefilternew.h )


#########################
If CUDA is found, then add CUDA file and add to library
###########################

if (BUILD_CUDA)
    SET (SRCS ${SRCS} LocalMaxima3D.cu )
    CUDA_ADD_LIBRARY( seedsDetection_3D ${SRCS} )

else (BUILD_CUDA)
    ADD_LIBRARY( seedsDetection_3D ${SRCS} )

endif (BUILD_CUDA)

TARGET_LINK_LIBRARIES ( seedsDetection_3D        ${ITK_LIBRARIES} )

#INSTALL( TARGETS seedsDetection_3D DESTINATION . )

