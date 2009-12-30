ADD_EXECUTABLE( segment_nuclei yousef_main.cpp)
ADD_EXECUTABLE( compute_nuclei_features compute_nuclei_features_main.cpp)
ADD_EXECUTABLE( compute_associative_measures compute_associative_measures.cpp)
ADD_EXECUTABLE( classify_nuclei split_by_class.cpp )
ADD_EXECUTABLE( graph_coloring_on_output graph_coloring_on_output.cpp )
ADD_EXECUTABLE( convert_dat_to_label convert_dat_to_label.cpp )

TARGET_LINK_LIBRARIES( segment_nuclei  Yousef_Nucleus_Seg)# ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( compute_nuclei_features NuclearSegmentation ftkCommon)# ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( compute_associative_measures ftkCommon NuclearSegmentation)#  ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( classify_nuclei ftkCommon NuclearSegmentation )
TARGET_LINK_LIBRARIES( graph_coloring_on_output ftkCommon NuclearSegmentation )
TARGET_LINK_LIBRARIES( convert_dat_to_label ftkCommon NuclearSegmentation )

INSTALL( TARGETS segment_nuclei DESTINATION ${INSTALL_BIN_DIR} )
INSTALL( TARGETS compute_nuclei_features DESTINATION ${INSTALL_BIN_DIR} )
INSTALL( TARGETS compute_associative_measures DESTINATION ${INSTALL_BIN_DIR} )
INSTALL( TARGETS classify_nuclei DESTINATION ${INSTALL_BIN_DIR} )
INSTALL( TARGETS convert_dat_to_label DESTINATION ${INSTALL_BIN_DIR} )

#ADD_EXECUTABLE( nucseg2 extended_main.cpp )
#TARGET_LINK_LIBRARIES( nucseg2 nuclei )



