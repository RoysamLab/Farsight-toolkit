SET( MODELSEGMENTATION_SRCS
  model_nucleus_seg.cpp
)

SET( MODELSEGMENTATION_HDRS
  model_nucleus_seg.h
)


SET( GLPK_HDRS
  glpk.h	
)


ADD_LIBRARY( ModelSeg
                ${MODELSEGMENTATION_SRCS}
		${GLPK_HDRS}
		${MODELSEGMENTATION_HDRS} 
	      )



TARGET_LINK_LIBRARIES( ModelSeg ftkGUI ftkImage ftkCommon ftkFeatures libSVM NuclearSegmentation Nuclear_Association libagf ${LIB_GLPK} CytoplasmSegmentation ${ITK_LIBRARIES} )

