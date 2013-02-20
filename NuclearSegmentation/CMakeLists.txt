SET( NUCLEARSEGMENTATION_SRCS
  ftkNuclearSegmentation.cpp
)

SET( NUCLEARSEGMENTATION_HDRS
  ftkNuclearSegmentation.h
)

INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/SQLite/NESqlite )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

IF(USE_TRACKING)
  ADD_DEFINITIONS(-DUSE_TRACKING)
  INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/CellTrackerLib )
  SET(TRACK_LIB CellTrackerLib)
  IF (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:msvcrt.lib")
  ENDIF (MSVC)
ELSE(USE_TRACKING)
  SET(TRACK_LIB )
ENDIF(USE_TRACKING)

ADD_LIBRARY( NuclearSegmentation
            ${NUCLEARSEGMENTATION_SRCS}
            ${NUCLEARSEGMENTATION_HDRS} 
)

TARGET_LINK_LIBRARIES( NuclearSegmentation ftkCommon ftkFeatures ftkImage Yousef_Nucleus_Seg Nuclear_Association PixelLevelAnalysis ${TRACK_LIB} )
#INSTALL( TARGETS NuclearSegmentation DESTINATION . )



#Remove SegmentationViewer (segren) because it does not compile on Windows 7 machine.
#ADD_SUBDIRECTORY( yousef_core Nuclear_Association exe CytoplasmSegmentation PixelAnalysis NucleusEditor SegmentationViewer )
ADD_SUBDIRECTORY( yousef_core )
ADD_SUBDIRECTORY( Nuclear_Association )
ADD_SUBDIRECTORY( exe )
ADD_SUBDIRECTORY( CytoplasmSegmentation )
ADD_SUBDIRECTORY( PixelAnalysis )
ADD_SUBDIRECTORY( NucleusEditor )
ADD_SUBDIRECTORY( SubtractChannel )
# Statistical_Model_Segmentation )

IF (BUILD_OPENCL)
    IF (MSVC)
	INCLUDE_DIRECTORIES( $ENV{CUDA_INC_PATH} )
	LINK_DIRECTORIES( $ENV{CUDA_LIB_PATH} )
	TARGET_LINK_LIBRARIES(NuclearSegmentation OpenCL.lib)
    ENDIF(MSVC)
ENDIF(BUILD_OPENCL)

IF (BUILD_MODEL_SEG)
ADD_SUBDIRECTORY(Statistical_Model_Segmentation)
ENDIF(BUILD_MODEL_SEG)
