SET( cb3_SRCS
	alpha_expansion_3d.h
	alpha_expansion_3d.cpp
	GCoptimization.h
	GCoptimization.cpp	
	LinkedBlockList.h
	LinkedBlockList.cpp
	#GraphCut3dConstr.cpp	
	#GraphCutMex.cpp
	graph.cpp
	maxflow.cpp
	block.h
	graph.h
	GraphCut.h
)

ADD_LIBRARY( alpha_expansion_3d ${cb3_SRCS} )
TARGET_LINK_LIBRARIES ( alpha_expansion_3d        ${ITK_LIBRARIES} )
#INSTALL( TARGETS alpha_expansion_3d DESTINATION . )

