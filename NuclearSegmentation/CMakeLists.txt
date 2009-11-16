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

TARGET_LINK_LIBRARIES( NuclearSegmentation ftkCommon ftkImage Yousef_Nucleus_Seg Nuclear_Association)
#INSTALL( TARGETS NuclearSegmentation DESTINATION . )

SUBDIRS( yousef_core  Nuclear_Association exe NucleusEditor Histopathology SegmentationViewer )

