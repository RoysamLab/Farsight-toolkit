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

ADD_LIBRARY(mdlUtils mdlUtils.h mdlUtils.cpp)

ADD_LIBRARY(mdlVolumeProcess mdlVolumeProcess.h mdlVolumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(mdlVolumeProcess ITKIO ITKStatistics ITKBasicFilters ITKCommon)

#ADD_EXECUTABLE(mdlVolumeProcess2 volproc_main.cpp)
#TARGET_LINK_LIBRARIES(volumeProcess2 mdlVolumeProcess)

ADD_LIBRARY(mdlIntegratedSkeleton mdlIntegratedSkeleton.h mdlIntegratedSkeleton.cpp)
TARGET_LINK_LIBRARIES(mdlIntegratedSkeleton mdlUtils ITKBasicFilters ITKCommon ITKIO ITKStatistics )

ADD_EXECUTABLE(mdl_skeletonize skeleton_main.cpp)
TARGET_LINK_LIBRARIES(mdl_skeletonize mdlIntegratedSkeleton)

ADD_LIBRARY(mdlMST mdlMST.h mdlMST.cpp WeightedMahalsnobisDistance.cxx)
TARGET_LINK_LIBRARIES(mdlMST mdlUtils mdlVolumeProcess)

ADD_EXECUTABLE(mdl_mst mst_main.cpp)
TARGET_LINK_LIBRARIES(mdl_mst mdlMST)

ADD_EXECUTABLE(new_mst_try new_mst_try.cpp)
TARGET_LINK_LIBRARIES(new_mst_try mdlMST)

ADD_LIBRARY(mdlBSplineFitting mdlBSplineFitting.cpp mdlBSplineFitting.h)

#ADD_EXECUTABLE(mdlpipe pipeline_main.cpp)
#TARGET_LINK_LIBRARIES(mdlpipe mdlVolumeProcess mdlIntegratedSkeleton mdlMST mdlBSplineFitting  )

#ADD_EXECUTABLE(mdlpipe2 pipeline_main2.cpp)
#TARGET_LINK_LIBRARIES(mdlpipe2 mdlVolumeProcess mdlIntegratedSkeleton mdlMST mdlBSplineFitting  )

#ADD_EXECUTABLE(PreDistanceMAP PredistanceMAP.cpp)
#TARGET_LINK_LIBRARIES(PreDistanceMAP mdlVolumeProcess)

SUBDIRS( wizard )

