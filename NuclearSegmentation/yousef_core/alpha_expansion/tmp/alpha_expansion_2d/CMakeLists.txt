#cmake_minimum_required (VERSION 2.6) 

#Project( alpha_expansion_2d )

SET( cb2_SRCS
	block.h
	graph.h
	GraphCut.h	
	gvc.h
	LinkedBlockList.cpp
	#GraphCutConstr.cpp	
	#GraphCutMex.cpp
	graph.cpp
	maxflow.cpp
	GCoptimization.h
	LinkedBlockList.h
	Multi_Color_Graph_Learning_2D.h
	GCoptimization.cpp	
	alpha_expansion_2d.cpp
	Multi_Color_Graph_Learning_2d.cpp
	gvc.cpp			
)


ADD_LIBRARY( alpha_expansion_2d ${cb2_SRCS} )
TARGET_LINK_LIBRARIES ( alpha_expansion_2d      ${ITK_LIBRARIES} )
#INSTALL( TARGETS alpha_expansion_3d DESTINATION . )
