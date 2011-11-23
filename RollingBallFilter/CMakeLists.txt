PROJECT(RollingBallFilter)

#INCLUDE_DIRECTORIES()

SET(RollingBallFilter_HDRS 
)

SET(RollingBallFilter_SRCS 
	main.cpp
)

ADD_EXECUTABLE(RollingBallFilter main.cpp ${RollingBallFilter_HDRS} ${RollingBallFilter_SRCS})
TARGET_LINK_LIBRARIES(RollingBallFilter ${ITK_LIBRARIES})
INSTALL( TARGETS RollingBallFilter RUNTIME DESTINATION bin )