SET( SRCS
	Local_Max_Clust_2D.h
	Local_Max_Clust_2D.cpp
)

ADD_LIBRARY( clustering_2D ${SRCS} )
TARGET_LINK_LIBRARIES ( clustering_2D        ${ITK_LIBRARIES} )
#INSTALL( TARGETS clustering_2D DESTINATION . )
