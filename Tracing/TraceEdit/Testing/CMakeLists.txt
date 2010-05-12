
include_directories(${TraceEdit_SOURCE_DIR})

ADD_EXECUTABLE(TraceEdit-StartStop TraceEditStartStop.cxx)
TARGET_LINK_LIBRARIES(TraceEdit-StartStop Trace)
#ADD_TEST(TraceEdit-StartStop TraceEdit-StartStop)

