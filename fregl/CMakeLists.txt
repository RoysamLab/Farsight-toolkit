#IF(vxl_SOURCE_DIR)
#  SET(VXL_MUL_INCLUDE_DIR ${MUL_INCLUDE_DIR})
#ENDIF(vxl_SOURCE_DIR)

#IF( NOT ITK_USE_SYSTEM_VXL AND NOT USE_SUPERBUILD)
    #MESSAGE(FATAL_ERROR "ITK must be built with ITK_USE_SYSTEM_VXL set ON.")
#ENDIF( NOT ITK_USE_SYSTEM_VXL AND NOT USE_SUPERBUILD)

INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )

SET( fregl_sources
	fregl_reg_record.h          	fregl_reg_record.cxx
	fregl_pairwise_register.h   	fregl_pairwise_register.cxx
	fregl_joint_register.h      	fregl_joint_register.cxx
	fregl_space_transformer.h	fregl_space_transformer.cxx
  
	fregl_image_manager.h       	fregl_image_manager.cxx
	#fregl_result_record.h       	fregl_result_record.cxx
 	fregl_util.h                	fregl_util.cxx
 	fregl_roi.h			fregl_roi.cxx		
)

AUX_SOURCE_DIRECTORY( Templates fregl_sources )

ADD_LIBRARY( fregl ${fregl_sources} )

TARGET_LINK_LIBRARIES( fregl ftkImage rgrl vil3d vnl vbl  TinyXML ${ITK_LIBRARIES} ${PTHREAD_LIBRARIES} ${QT_LIBRARIES})

#IF(BUILD_TESTING)
#  ADD_SUBDIRECTORY(tests)
#ENDIF(BUILD_TESTING)

INSTALL( TARGETS fregl DESTINATION ${INSTALL_BIN_DIR} )

#ADD_SUBDIRECTORY( GUI )
ADD_SUBDIRECTORY( exe )
