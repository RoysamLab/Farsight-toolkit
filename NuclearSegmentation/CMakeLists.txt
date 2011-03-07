SET( NUCLEARSEGMENTATION_SRCS
  ftkNuclearSegmentation.cpp
)

SET( NUCLEARSEGMENTATION_HDRS
  ftkNuclearSegmentation.h
)


ADD_LIBRARY( NuclearSegmentation
            ${NUCLEARSEGMENTATION_SRCS}
            ${NUCLEARSEGMENTATION_HDRS} 
)

TARGET_LINK_LIBRARIES( NuclearSegmentation ftkCommon ftkFeatures ftkImage Yousef_Nucleus_Seg Nuclear_Association PixelLevelAnalysis )
#INSTALL( TARGETS NuclearSegmentation DESTINATION . )

#Remove SegmentationViewer (segren) because it does not compile on Windows 7 machine.
#SUBDIRS( yousef_core Nuclear_Association exe CytoplasmSegmentation PixelAnalysis NucleusEditor SegmentationViewer )
SUBDIRS( yousef_core Nuclear_Association exe CytoplasmSegmentation PixelAnalysis NucleusEditor SubtractChannel )

IF (BUILD_OPENCL)
    IF (MSVC)
        INCLUDE_DIRECTORIES( $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/inc/ )
	INCLUDE_DIRECTORIES( $ENV{NVSDKCOMPUTE_ROOT}/shared/inc )
	LINK_DIRECTORIES( $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/lib )
	LINK_DIRECTORIES( $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/lib/x64 )
	LINK_DIRECTORIES( $ENV{NVSDKCOMPUTE_ROOT}/shared/lib )
	TARGET_LINK_LIBRARIES(NuclearSegmentation debug oclUtils64D.lib optimized oclUtils64.lib )
	TARGET_LINK_LIBRARIES(NuclearSegmentation debug shrUtils64D.lib optimized shrUtils64.lib )
	TARGET_LINK_LIBRARIES(NuclearSegmentation OpenCL.lib)
    ENDIF(MSVC)
ENDIF(BUILD_OPENCL)