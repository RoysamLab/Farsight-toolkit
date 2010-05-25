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
SUBDIRS( yousef_core Nuclear_Association exe CytoplasmSegmentation PixelAnalysis NucleusEditor )