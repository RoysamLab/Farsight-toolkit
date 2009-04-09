ADD_EXECUTABLE( segment_nuclei yousef_main.cpp)
ADD_EXECUTABLE( compute_nuclei_features compute_nuclei_features_main.cpp)
ADD_EXECUTABLE( compute_associative_measures compute_associative_measures.cpp)
ADD_EXECUTABLE( classify_nuclei split_by_class.cpp )
ADD_EXECUTABLE( graph_coloring_on_output graph_coloring_on_output.cpp )

TARGET_LINK_LIBRARIES( segment_nuclei  Yousef_Nucleus_Seg)# ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( compute_nuclei_features NuclearSegmentation SegmentationCommon)# ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( compute_associative_measures SegmentationCommon NuclearSegmentation)#  ITKBasicFilters ITKCommon ITKIO)
TARGET_LINK_LIBRARIES( classify_nuclei SegmentationCommon NuclearSegmentation )
TARGET_LINK_LIBRARIES( graph_coloring_on_output SegmentationCommon NuclearSegmentation )

INSTALL( TARGETS segment_nuclei DESTINATION bin )
INSTALL( TARGETS compute_nuclei_features DESTINATION bin )
INSTALL( TARGETS compute_associative_measures DESTINATION bin )
INSTALL( TARGETS classify_nuclei DESTINATION bin )

#ADD_EXECUTABLE( nucseg2 extended_main.cpp )
#TARGET_LINK_LIBRARIES( nucseg2 nuclei )



