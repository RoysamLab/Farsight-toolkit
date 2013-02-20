PROJECT(MDL2)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL2 ${Farsight_SOURCE_DIR}/Tracing/MDL2 )

SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)

ADD_LIBRARY(mdlUtils mdlUtils.h mdlUtils.cpp)

ADD_LIBRARY(mdlVolumeProcess mdlVolumeProcess.h mdlVolumeProcess.cpp ${cbin_SRCS})
TARGET_LINK_LIBRARIES(mdlVolumeProcess         ${ITK_LIBRARIES} )

ADD_LIBRARY(mdlIntegratedSkeleton mdlIntegratedSkeleton.h mdlIntegratedSkeleton.cpp)
TARGET_LINK_LIBRARIES(mdlIntegratedSkeleton mdlUtils          ${ITK_LIBRARIES} )

ADD_EXECUTABLE(mdl_skeletonize skeleton_main.cpp)
TARGET_LINK_LIBRARIES(mdl_skeletonize mdlIntegratedSkeleton)

ADD_LIBRARY(mdlMST mdlMST.h mdlMST.cpp WeightedMahalsnobisDistance.cxx)
TARGET_LINK_LIBRARIES(mdlMST mdlUtils mdlVolumeProcess)

ADD_EXECUTABLE(mdl_mst mst_main.cpp)
TARGET_LINK_LIBRARIES(mdl_mst mdlMST)

ADD_LIBRARY(mdlBSplineFitting mdlBSplineFitting.cpp mdlBSplineFitting.h)
TARGET_LINK_LIBRARIES(mdlBSplineFitting mdlUtils      ${ITK_LIBRARIES} )

ADD_SUBDIRECTORY( wizard )

