if(QtTestingFound)
  if(APPLE)
    SET(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor.app/Contents/MacOS/NucleusEditor)
  else()
    set(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor)
  endif()

  set(NucleusEditor_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(
    NucleusEditor-SegmentNuclei
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/segment_nuclei.xml
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
endif()

