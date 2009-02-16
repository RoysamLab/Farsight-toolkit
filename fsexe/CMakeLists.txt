
#ADD_EXECUTABLE( nucleus_seg nucleus_seg.cxx)
#TARGET_LINK_LIBRARIES( nucleus_seg maciejseg fscommon ITKCommon ITKIO)

ADD_EXECUTABLE( fuse_image_stacks fuse_image_stacks.cxx)
TARGET_LINK_LIBRARIES( fuse_image_stacks fscommon ITKCommon ITKIO)

#ADD_EXECUTABLE( convert_seg_result_from_idl convert_seg_result_from_idl.cxx)
#TARGET_LINK_LIBRARIES( convert_seg_result_from_idl maciejseg vul
#ITKCommon ITKIO)

#ADD_EXECUTABLE( separate_nuclei_by_class separate_nuclei_by_class.cxx)
#TARGET_LINK_LIBRARIES( separate_nuclei_by_class maciejseg vul vnl)

ADD_EXECUTABLE( lsm_to_tiff lsm_to_tiff.cxx)
TARGET_LINK_LIBRARIES( lsm_to_tiff FTKImage ITKCommon ITKIO)

#ADD_EXECUTABLE( vessel_segmentation vessel_segmentation.cpp)

#ADD_EXECUTABLE( npts_to_tiff npts_to_tiff.cxx)
#TARGET_LINK_LIBRARIES( npts_to_tiff vnl ITKCommon ITKIO)

#ADD_EXECUTABLE( transform_segmentation_results transform_segmentation_results.cxx)
#TARGET_LINK_LIBRARIES( transform_segmentation_results maciejseg fregl
#vnl )

ADD_EXECUTABLE( generate_color_stack generate_color_stack.cxx)
TARGET_LINK_LIBRARIES( generate_color_stack ITKCommon ITKIO )

#ADD_EXECUTABLE( assign_signals assign_signals.cxx)
#TARGET_LINK_LIBRARIES( assign_signals maciejseg vnl vul ITKCommon ITKIO )

ADD_EXECUTABLE( max_project_multi_channels max_project_multi_channels.cxx)
TARGET_LINK_LIBRARIES( max_project_multi_channels vnl ITKCommon ITKIO )

ADD_EXECUTABLE( multi_channels_2D multi_channels_2D.cxx)
TARGET_LINK_LIBRARIES( multi_channels_2D vnl ITKCommon ITKIO )

ADD_EXECUTABLE( montage_blending montage_blending.cxx)
TARGET_LINK_LIBRARIES( montage_blending fscommon fregl vnl ITKCommon ITKIO )

ADD_EXECUTABLE( extract_stack_interfaces extract_stack_interfaces.cxx)
TARGET_LINK_LIBRARIES( extract_stack_interfaces fscommon fregl vul ITKCommon ITKIO )

#ADD_EXECUTABLE( update_features update_features.cxx)
#TARGET_LINK_LIBRARIES( update_features maciejseg vul ITKCommon ITKIO )

#IF (BUILD_FREGL)
   ADD_EXECUTABLE( register_pair register_pair.cxx)
   TARGET_LINK_LIBRARIES( register_pair fregl fscommon ITKCommon ITKIO )

   ADD_EXECUTABLE( register_joint register_joint.cxx)
   TARGET_LINK_LIBRARIES( register_joint fregl fscommon ITKCommon
   ITKIO )

   ADD_EXECUTABLE( mosaic_image_pair mosaic_image_pair.cxx)
   TARGET_LINK_LIBRARIES( mosaic_image_pair fregl fscommon ITKCommon
   ITKIO )

   ADD_EXECUTABLE( mosaic_images mosaic_images.cxx)
   TARGET_LINK_LIBRARIES( mosaic_images fregl fscommon ITKCommon
   ITKIO )

   #ADD_EXECUTABLE( update_result_sets update_result_sets.cxx)
   #TARGET_LINK_LIBRARIES( update_result_sets fregl )

   #ADD_EXECUTABLE( nuclear_montage nuclear_montage.cxx)
   #TARGET_LINK_LIBRARIES( nuclear_montage maciejseg fregl vul)

   #ADD_EXECUTABLE( vessel_montage vessel_montage.cxx)
   #TARGET_LINK_LIBRARIES( vessel_montage fscommon maciejseg fregl vul)

   #ADD_EXECUTABLE( update_nearest_neighbor_graph update_nearest_neighbor_graph.cxx)
   #TARGET_LINK_LIBRARIES( update_nearest_neighbor_graph maciejseg
   #fregl vul)

   #ADD_EXECUTABLE( query_distance_to_top query_distance_to_top.cxx)
   #TARGET_LINK_LIBRARIES( query_distance_to_top maciejseg fregl vul)

   #ADD_EXECUTABLE( layer_analysis layer_analysis.cxx)
   #TARGET_LINK_LIBRARIES( layer_analysis maciejseg fregl vul
   #ITKCommon ITKIO )

   #ADD_EXECUTABLE( register_parallel register_parallel.cxx)
   #TARGET_LINK_LIBRARIES( register_parallel fregl)

   ADD_EXECUTABLE( SubsampleVolume SubsampleVolume.cxx)
   TARGET_LINK_LIBRARIES( SubsampleVolume ITKCommon ITKIO )
#ENDIF(BUILD_FREGL)

INSTALL( TARGETS lsm_to_tiff register_pair register_joint mosaic_image_pair mosaic_images
DESTINATION bin )

#copy the gdbicp to the installation directory
#IF (WIN32) 
#   CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/fsexe/gdbicp_win.exe
#                  ${PROJECT_SOURCE_DIR}/fsexe/gdbicp.exe)
#   INSTALL( TARGETS ${PROJECT_SOURCE_DIR}/fsexe/gdbicp.exe DESTINATION . )
#ELSE (WIN32) 
#   CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/fsexe/gdbicp_unix
#                  ${PROJECT_SOURCE_DIR}/fsexe/gdbicp)
#   INSTALL( TARGETS ${PROJECT_SOURCE_DIR}/fsexe/gdbicp DESTINATION . )
#ENDIF (WIN32)

IF (WIN32)
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/fsexe/gdbicp_win.exe
    ${PROJECT_BINARY_DIR}/fsexe/gdbicp.exe
    COPYONLY)
  INSTALL( PROGRAMS ${PROJECT_BINARY_DIR}/fsexe/gdbicp.exe DESTINATION bin )
ELSE (WIN32)
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/fsexe/gdbicp_unix
    ${PROJECT_BINARY_DIR}/fsexe/gdbicp
    COPYONLY)
  INSTALL( PROGRAMS ${PROJECT_BINARY_DIR}/fsexe/gdbicp DESTINATION bin )
ENDIF (WIN32) 

