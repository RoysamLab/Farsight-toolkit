PROJECT(Batch_Distance2Device)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Tracing/TraceEdit)

SET(Batch_Distance2Device_HDRS 
)

SET(Batch_Distance2Device_SRCS 
	main.cpp
)

ADD_EXECUTABLE(Batch_Distance2Device ${Batch_Distance2Device_HDRS} ${Batch_Distance2Device_SRCS})
TARGET_LINK_LIBRARIES(Batch_Distance2Device Trace)
INSTALL( TARGETS Batch_Distance2Device RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
