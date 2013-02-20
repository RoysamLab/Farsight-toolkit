
#ADD_EXECUTABLE( nucleus_seg nucleus_seg.cxx)
#TARGET_LINK_LIBRARIES( nucleus_seg maciejseg fscommon  ${ITK_LIBRARIES} )

ADD_EXECUTABLE( fuse_image_stacks fuse_image_stacks.cxx)
TARGET_LINK_LIBRARIES( fuse_image_stacks  fregl ${ITK_LIBRARIES} )

#ADD_EXECUTABLE( convert_seg_result_from_idl convert_seg_result_from_idl.cxx)
#TARGET_LINK_LIBRARIES( convert_seg_result_from_idl maciejseg vul
#${ITK_LIBRARIES})

#ADD_EXECUTABLE( separate_nuclei_by_class separate_nuclei_by_class.cxx)
#TARGET_LINK_LIBRARIES( separate_nuclei_by_class maciejseg vul vnl)

IF( FARSIGHT_VTK_6 )
  SET( LSM_TO_TIFF_EXEC  )
ELSE( FARSIGHT_VTK_6 )
  ADD_EXECUTABLE( lsm_to_tiff lsm_to_tiff.cxx)
  TARGET_LINK_LIBRARIES( lsm_to_tiff ftkImage  ${ITK_LIBRARIES} )
  SET( LSM_TO_TIFF_EXEC lsm_to_tiff )
ENDIF( FARSIGHT_VTK_6 )

#ADD_EXECUTABLE( vessel_segmentation vessel_segmentation.cpp)

#ADD_EXECUTABLE( npts_to_tiff npts_to_tiff.cxx)
#TARGET_LINK_LIBRARIES( npts_to_tiff vnl  ${ITK_LIBRARIES} )

ADD_EXECUTABLE( transform_segmentation_results transform_segmentation_results.cxx)
TARGET_LINK_LIBRARIES( transform_segmentation_results fregl vnl )

ADD_EXECUTABLE( extract_translations extract_translations.cxx)
TARGET_LINK_LIBRARIES( extract_translations fregl vnl )

ADD_EXECUTABLE( manual_translation manual_translation.cxx)
TARGET_LINK_LIBRARIES( manual_translation fregl vul   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( generate_color_stack generate_color_stack.cxx)
TARGET_LINK_LIBRARIES( generate_color_stack   ${ITK_LIBRARIES} )

#ADD_EXECUTABLE( assign_signals assign_signals.cxx)
#TARGET_LINK_LIBRARIES( assign_signals maciejseg vnl vul   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( max_project_multi_channels max_project_multi_channels.cxx)
TARGET_LINK_LIBRARIES( max_project_multi_channels vnl   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( multi_channels_2D multi_channels_2D.cxx)
TARGET_LINK_LIBRARIES( multi_channels_2D vnl   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( multi_channels_3D multi_channels_3D.cxx)
TARGET_LINK_LIBRARIES( multi_channels_3D vnl   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( montage_blending montage_blending.cxx)
TARGET_LINK_LIBRARIES( montage_blending fregl vil vul vnl   ${ITK_LIBRARIES} )

ADD_EXECUTABLE( extract_stack_interfaces extract_stack_interfaces.cxx)
TARGET_LINK_LIBRARIES( extract_stack_interfaces fregl vul   ${ITK_LIBRARIES} )

#ADD_EXECUTABLE( update_features update_features.cxx)
#TARGET_LINK_LIBRARIES( update_features maciejseg vul   ${ITK_LIBRARIES} )

#IF (BUILD_FREGL)
   ADD_EXECUTABLE( register_pair register_pair.cxx Templates/register_pair_template.cxx Templates/register_pair_template.h)
   TARGET_LINK_LIBRARIES( register_pair fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( register_pair_16 register_pair_16.cxx Templates/register_pair_template.cxx Templates/register_pair_template.h)
   TARGET_LINK_LIBRARIES( register_pair_16 fregl  ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( register_joint register_joint.cxx Templates/register_joint_template.cxx Templates/register_joint_template.h)
   TARGET_LINK_LIBRARIES( register_joint fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( register_joint_16 register_joint_16.cxx Templates/register_joint_template.cxx Templates/register_joint_template.h)
   TARGET_LINK_LIBRARIES( register_joint_16 fregl   ${ITK_LIBRARIES} )
   
   ADD_EXECUTABLE( mosaic_image_pair mosaic_image_pair.cxx)
   TARGET_LINK_LIBRARIES( mosaic_image_pair fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( mosaic_images mosaic_images.cxx Templates/mosaic_images_template.cxx Templates/mosaic_images_template.h)
   TARGET_LINK_LIBRARIES( mosaic_images fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( mosaic_images_16 mosaic_images_16.cxx Templates/mosaic_images_template.cxx Templates/mosaic_images_template.h)
   TARGET_LINK_LIBRARIES( mosaic_images_16 fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( mosaic_roi mosaic_roi.cxx)
   TARGET_LINK_LIBRARIES( mosaic_roi fregl   ${ITK_LIBRARIES} )

   ADD_EXECUTABLE( initialized_register_pair initialized_register_pair.cxx)
   TARGET_LINK_LIBRARIES( initialized_register_pair fregl   ${ITK_LIBRARIES} )

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
   #${ITK_LIBRARIES} )

   #ADD_EXECUTABLE( register_parallel register_parallel.cxx)
   #TARGET_LINK_LIBRARIES( register_parallel fregl)

   ADD_EXECUTABLE( SubsampleVolume SubsampleVolume.cxx)
   TARGET_LINK_LIBRARIES( SubsampleVolume   ${ITK_LIBRARIES} )
#ENDIF(BUILD_FREGL)

INSTALL( TARGETS ${LSM_TO_TIFF_EXEC} register_pair register_pair_16 register_joint register_joint_16 mosaic_image_pair mosaic_images mosaic_images_16 mosaic_roi
DESTINATION ${INSTALL_BIN_DIR} )

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
    ${PROJECT_SOURCE_DIR}/fregl/exe/gdbicp_win.exe
    ${PROJECT_BINARY_DIR}/fregl/exe/gdbicp.exe
    COPYONLY)
  INSTALL( PROGRAMS ${PROJECT_BINARY_DIR}/fregl/exe/gdbicp.exe DESTINATION ${INSTALL_BIN_DIR} )
ELSE (WIN32)
  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/fregl/exe/gdbicp_unix
    ${PROJECT_BINARY_DIR}/fregl/exe/gdbicp
    COPYONLY)
  # This executable was built against shared libraries.  Install fails on Linux
  # because fixup_bundle can't find these libraries.
  IF(NOT UNIX)
    INSTALL( PROGRAMS ${PROJECT_BINARY_DIR}/fregl/exe/gdbicp DESTINATION ${INSTALL_BIN_DIR} )
  ENDIF()
ENDIF (WIN32) 

