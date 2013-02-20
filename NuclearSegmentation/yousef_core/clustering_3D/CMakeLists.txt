SET( SRCS
	local_max_clust_3D.h
	Local_Max_Clust_3D_CPP.cpp
  Local_Max_Clust_3D.cu)


ADD_LIBRARY( clustering_3D ${SRCS} )

TARGET_LINK_LIBRARIES ( clustering_3D        ${ITK_LIBRARIES} )

#INSTALL( TARGETS clustering_3D DESTINATION . )
