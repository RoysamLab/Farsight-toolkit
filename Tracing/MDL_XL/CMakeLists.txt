PROJECT(MDL_XL)
FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} ${Farsight_BINARY_DIR}/Tracing/MDL_XL ${Farsight_SOURCE_DIR}/Tracing/MDL_XL )

SET( cbin_SRCS
  GCBinarization/cell_binarization.h
  GCBinarization/cell_binarization.cpp    
  GCBinarization/new_graph.h
  GCBinarization/block_b.h   
)

#CHECK TO MAKE SURE ITK LIBRARIES ARE NEEDED FOR THESE:

SET( xl_srcs
  ${cbin_SRCS}
  mdlTypes.h 
  mdlUtils.h mdlUtils.cpp
  mdlVolumeProcess.h mdlVolumeProcess.cpp
  mdlIntegratedSkeleton.h mdlIntegratedSkeleton.cpp
  mdlMST.h mdlMST.cpp WeightedMahalsnobisDistance.cxx
  mdlBSplineFitting.cpp mdlBSplineFitting.h
  pipeline_main2.cpp
)


ADD_EXECUTABLE(mdl_XL_pipe ${xl_srcs} )
TARGET_LINK_LIBRARIES(mdl_XL_pipe         ${ITK_LIBRARIES} )

ADD_EXECUTABLE(Demo_Bspline Demo_Bspline.cpp  mdlBSplineFitting.cpp mdlBSplineFitting.h   mdlTypes.h 
  mdlUtils.h mdlUtils.cpp)
TARGET_LINK_LIBRARIES(Demo_Bspline mdlUtils          ${ITK_LIBRARIES} )

