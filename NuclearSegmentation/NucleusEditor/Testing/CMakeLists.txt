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
endif()

