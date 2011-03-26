SET( MODELSEGMENTATION_SRCS
  model_nucleus_seg.cpp
)

SET( MODELSEGMENTATION_HDRS
  model_nucleus_seg.h
)


SET( GRAPHS_SRCS
  seg_graphs.cpp
)

SET( GRAPHS_HDRS
  seg_graphs.h
)

SET( GLPK_HDRS
  glpk.h	
)


ADD_LIBRARY( ModelSeg
                ${MODELSEGMENTATION_SRCS}
		${GLPK_HDRS}
		${MODELSEGMENTATION_HDRS} 
	      )



ADD_LIBRARY( GraphsnTrees
                ${GRAPHS_SRCS} 
                ${GRAPHS_HDRS}
	      )

	
TARGET_LINK_LIBRARIES( ModelSeg GraphsnTrees ftkGUI ftkImage ftkCommon ftkFeatures libSVM NuclearSegmentation Nuclear_Association ITKBasicFilters ITKCommon ITKIO libagf GLPKLib CytoplasmSegmentation )


ADD_EXECUTABLE(modelSegExe StatModelSegMain.cpp)
TARGET_LINK_LIBRARIES(modelSegExe ModelSeg GraphsnTrees)  

