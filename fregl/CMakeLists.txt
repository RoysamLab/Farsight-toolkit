INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

SET( fregl_sources
	fregl_reg_record.h          	fregl_reg_record.cxx
	fregl_pairwise_register.h   	fregl_pairwise_register.cxx
	fregl_joint_register.h      	fregl_joint_register.cxx
	fregl_space_transformer.h	     fregl_space_transformer.cxx
	fregl_image_manager.h       	fregl_image_manager.cxx
 	fregl_util.h                	fregl_util.cxx
 	fregl_roi.h			     fregl_roi.cxx
)

AUX_SOURCE_DIRECTORY( Templates fregl_sources )

ADD_LIBRARY( fregl ${fregl_sources} )

TARGET_LINK_LIBRARIES( fregl ftkImage rgrl vil3d vnl vbl  TinyXML ${ITK_LIBRARIES} ${QT_LIBRARIES})

INSTALL( TARGETS fregl DESTINATION ${INSTALL_BIN_DIR} )

ADD_SUBDIRECTORY( exe )
