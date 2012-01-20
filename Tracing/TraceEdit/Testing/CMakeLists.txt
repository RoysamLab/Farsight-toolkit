
include_directories(${TraceEdit_SOURCE_DIR})

ADD_EXECUTABLE(TraceEdit-StartStop TraceEditStartStop.cxx)
TARGET_LINK_LIBRARIES(TraceEdit-StartStop Trace)
ADD_TEST(TraceEdit-StartStop ${Farsight_BINARY_DIR}/exe/TraceEdit-StartStop)
SET_TESTS_PROPERTIES(TraceEdit-StartStop PROPERTIES
  FAIL_REGULAR_EXPRESSION "vtkDebugLeaks has detected LEAKS!"
)

IF(QtTestingFound)
  IF(APPLE)
    SET(TraceEdit_EXE ${Farsight_BINARY_DIR}/exe/TraceEdit.app/Contents/MacOS/TraceEdit) 
  ELSE()
    SET(TraceEdit_EXE ${Farsight_BINARY_DIR}/exe/TraceEdit)
  ENDIF()

  SET(TraceEdit_TEST_DIR ${TraceEdit_SOURCE_DIR}/Testing)

  ADD_TEST(
    TraceEdit-AutoMerge
    ${TraceEdit_EXE}
    ${FARSIGHT_DATA_ROOT}/TraceData/C5TracedPoints.swc
    ${TraceEdit_TEST_DIR}/test_automerge.xml
    ${TraceEdit_TEST_DIR}/baseline_automerge.png
  )
ENDIF()
