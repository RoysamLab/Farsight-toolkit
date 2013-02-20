ADD_EXECUTABLE( segment_nuclei yousef_main.cpp)

ADD_EXECUTABLE( segment_nuclei_16 yousef_main16.cpp)

ADD_EXECUTABLE( compute_nuclei_features compute_nuclei_features_main.cpp)

ADD_EXECUTABLE( compute_associative_measures compute_associative_measures.cpp)

ADD_EXECUTABLE( classify_nuclei split_by_class.cpp )
ADD_EXECUTABLE( graph_coloring_on_output graph_coloring_on_output.cpp )

ADD_EXECUTABLE( convert_dat_to_label convert_dat_to_label.cpp )

IF(BUILD_MODEL_SEG)
ADD_EXECUTABLE(MBSeg StatModelSegMain.cpp)
TARGET_LINK_LIBRARIES(MBSeg ModelSeg) 
ENDIF(BUILD_MODEL_SEG)

ADD_LIBRARY( segment_soma SomaExtraction.cpp SomaExtraction.h)
TARGET_LINK_LIBRARIES( segment_soma ftkCommon ${ITK_LIBRARIES} Yousef_Nucleus_Seg)
ADD_EXECUTABLE(SomaExtraction SOMA_MAIN.cpp)
TARGET_LINK_LIBRARIES( SomaExtraction  segment_soma )
INSTALL( TARGETS SomaExtraction DESTINATION ${INSTALL_BIN_DIR} )

TARGET_LINK_LIBRARIES( segment_nuclei  Yousef_Nucleus_Seg     ${ITK_LIBRARIES} )#      ${ITK_LIBRARIES} )

TARGET_LINK_LIBRARIES( segment_nuclei_16  Yousef_Nucleus_Seg     ${ITK_LIBRARIES} )#      ${ITK_LIBRARIES} )

TARGET_LINK_LIBRARIES( compute_nuclei_features NuclearSegmentation ftkCommon ftkFeatures     ${ITK_LIBRARIES} )#      ${ITK_LIBRARIES} )

TARGET_LINK_LIBRARIES( compute_associative_measures ftkCommon ftkFeatures NuclearSegmentation     ${ITK_LIBRARIES} )#       ${ITK_LIBRARIES} )

TARGET_LINK_LIBRARIES( classify_nuclei ftkCommon ftkFeatures NuclearSegmentation )

TARGET_LINK_LIBRARIES( graph_coloring_on_output ftkCommon ftkFeatures NuclearSegmentation )

TARGET_LINK_LIBRARIES( convert_dat_to_label ftkCommon NuclearSegmentation )



IF (BUILD_OPENCL)
    IF (MSVC)
	INCLUDE_DIRECTORIES( $ENV{CUDA_INC_PATH} )
	LINK_DIRECTORIES( $ENV{CUDA_LIB_PATH} )
	TARGET_LINK_LIBRARIES(segment_nuclei OpenCL.lib)
    ENDIF(MSVC)
ENDIF(BUILD_OPENCL)


INSTALL( TARGETS segment_nuclei DESTINATION ${INSTALL_BIN_DIR} )

INSTALL( TARGETS segment_nuclei_16 DESTINATION ${INSTALL_BIN_DIR} )

INSTALL( TARGETS compute_nuclei_features DESTINATION ${INSTALL_BIN_DIR} )

INSTALL( TARGETS compute_associative_measures DESTINATION ${INSTALL_BIN_DIR} )

INSTALL( TARGETS classify_nuclei DESTINATION ${INSTALL_BIN_DIR} )

INSTALL( TARGETS convert_dat_to_label DESTINATION ${INSTALL_BIN_DIR} )

#INSTALL( TARGETS MBSeg DESTINATION ${INSTALL_BIN_DIR} )


#ADD_EXECUTABLE( nucseg2 extended_main.cpp )

#TARGET_LINK_LIBRARIES( nucseg2 nuclei )



