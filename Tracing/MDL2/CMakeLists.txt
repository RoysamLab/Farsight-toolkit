PROJECT(MDL2)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL2 ${Farsight_SOURCE_DIR}/Tracing/MDL2 )

SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)

ADD_LIBRARY(mdlVolumeProcess mdlVolumeProcess.h mdlVolumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKIO ITKStatistics ITKBasicFilters ITKCommon)

#ADD_LIBRARY(mdlAnisoDiffuse mdlAnisoDiffuse.h mdlAnisoDiffuse.cpp)
#TARGET_LINK_LIBRARIES(mdlAnisoDiffuse ITKIO ITKStatistics)

ADD_EXECUTABLE(pp main.cpp)
TARGET_LINK_LIBRARIES(pp mdlVolumeProcess)