SET( SRCS
	EM_Project_3_cpp_3D_comp.cxx
	EM_Project_3_cpp_3D_comp.h
)

ADD_LIBRARY( EM_GMM ${SRCS} )
TARGET_LINK_LIBRARIES ( EM_GMM ITKBasicFilters ITKCommon ITKIO)
#INSTALL( TARGETS EM_GMM DESTINATION . )

