SET( SRCS
	local_max_clust_3D.h
	Local_Max_Clust_3D_CPP.cpp
)

ADD_LIBRARY( clustering_3D ${SRCS} )
TARGET_LINK_LIBRARIES ( clustering_3D ITKBasicFilters ITKCommon ITKIO)
#INSTALL( TARGETS clustering_3D DESTINATION . )
