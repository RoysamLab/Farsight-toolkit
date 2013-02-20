SET( SRCS
	Multi_Color_Graph_Learning_3D.h
	gvc.h
	Multi_Color_Graph_Learning_3D.cpp
	sequential_coloring.cpp
)

ADD_LIBRARY( graphColLearn_3D ${SRCS} )
TARGET_LINK_LIBRARIES ( graphColLearn_3D        ${ITK_LIBRARIES} )
#INSTALL( TARGETS graphColLearn_3D DESTINATION . )
