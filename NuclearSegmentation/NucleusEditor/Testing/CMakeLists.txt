if(QtTestingFound)
  if(APPLE)
    set(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor.app/Contents/MacOS/NucleusEditor)
    set(AppName "Farsight")
  else()
    set(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor)
    set(AppName "NucleusEditor")
  endif()

  set(NucleusEditor_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/segment_nuclei.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/segment_nuclei.xml
                 @ONLY)
  add_test(
    NucleusEditor-SegmentNuclei
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/segment_nuclei.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/segment_nuclei_baseline.png
  )
  
  add_test(
    NucleusEditor-ROIStatistics
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/ROI_statistics.xml
  )
  set_tests_properties(NucleusEditor-ROIStatistics PROPERTIES
    PASS_REGULAR_EXPRESSION "variance: 152.21"
  )
  
  add_test(
    NucleusEditor-SliceViewer
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/slice_viewer.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/slice_viewer_baseline.png
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/merge_segments.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/merge_segments.xml
                 @ONLY)
  add_test(
    NucleusEditor-MergeSegments
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/merge_segments.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/merge_segments_baseline.png
  )
  set_tests_properties(NucleusEditor-MergeSegments PROPERTIES
    DEPENDS NucleusEditor-SegmentNuclei
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/delete_segments.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/delete_segments.xml
                 @ONLY)
  add_test(
    NucleusEditor-DeleteSegments
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/delete_segments.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/delete_segments_baseline.png
  )
  set_tests_properties(NucleusEditor-DeleteSegments PROPERTIES
    DEPENDS NucleusEditor-SegmentNuclei
  )
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/show_centroids.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/show_centroids.xml
                 @ONLY)
  add_test(
    NucleusEditor-ShowCentroids
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${CMAKE_CURRENT_BINARY_DIR}/show_centroids.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/show_centroids_baseline.png
  )
  set_tests_properties(NucleusEditor-ShowCentroids PROPERTIES
    DEPENDS NucleusEditor-SegmentNuclei
  )

endif()

