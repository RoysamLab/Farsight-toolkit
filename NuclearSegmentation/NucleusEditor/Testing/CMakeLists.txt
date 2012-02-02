IF(QtTestingFound)
  IF(APPLE)
    SET(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor.app/Contents/MacOS/NucleusEditor)
  ELSE()
    SET(NucleusEditor_EXE ${Farsight_BINARY_DIR}/exe/NucleusEditor)
  ENDIF()

  SET(NucleusEditor_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  ADD_TEST(
    NucleusEditor-SegmentNuclei
    ${NucleusEditor_EXE}
    --image
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    --test
    ${NucleusEditor_TEST_DIR}/segment_nuclei.xml
    --baseline
    ${NucleusEditor_TEST_DIR}/segment_nuclei_baseline.png
  )
ENDIF()

