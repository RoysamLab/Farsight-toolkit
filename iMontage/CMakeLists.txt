PROJECT(iMontage)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Tracing/TraceEdit)

SET(iMontage_HDRS 
)

SET(iMontage_SRCS 
	main.cpp
)

ADD_EXECUTABLE(iMontage main.cpp ${iMontage_HDRS} ${iMontage_SRCS})
TARGET_LINK_LIBRARIES(iMontage Trace)
INSTALL( TARGETS iMontage RUNTIME DESTINATION bin )