PROJECT(MDL2)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL2 ${Farsight_SOURCE_DIR}/Tracing/MDL2 )

SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)

#CHECK TO MAKE SURE ITK LIBRARIES ARE NEEDED FOR THESE:

ADD_LIBRARY(mdlVolumeProcess mdlVolumeProcess.h mdlVolumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKIO ITKStatistics ITKBasicFilters ITKCommon)

ADD_LIBRARY(mdlIntegratedSkeleton mdlIntegratedSkeleton.h mdlIntegratedSkeleton.cpp)
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKBasicFilters ITKCommon)

ADD_LIBRARY(mdlUtils mdlUtils.h mdlUtils.cpp)

ADD_LIBRARY(mdlMST mdlMST.h mdlMST.cpp WeightedMahalsnobisDistance.cxx)
TARGET_LINK_LIBRARIES(mdlMST mdlUtils)

ADD_LIBRARY(mdlBSplineFitting mdlBSplineFitting.cpp mdlBSplineFitting.h)

ADD_EXECUTABLE(volumeProcess2 main.cpp)
TARGET_LINK_LIBRARIES(volumeProcess2 mdlVolumeProcess)

ADD_EXECUTABLE(mdlpipe pipeline_main.cpp)
TARGET_LINK_LIBRARIES(mdlpipe mdlVolumeProcess mdlIntegratedSkeleton mdlMST mdlBSplineFitting  )