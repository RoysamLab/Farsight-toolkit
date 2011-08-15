
include_directories(${TraceEdit_SOURCE_DIR})

ADD_EXECUTABLE(TraceEdit-StartStop TraceEditStartStop.cxx)
TARGET_LINK_LIBRARIES(TraceEdit-StartStop Trace)
ADD_TEST(TraceEdit-StartStop ${Farsight_BINARY_DIR}/exe/TraceEdit-StartStop)
SET_TESTS_PROPERTIES(TraceEdit-StartStop PROPERTIES
  FAIL_REGULAR_EXPRESSION "vtkDebugLeaks has detected LEAKS!"
)

