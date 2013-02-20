PROJECT(RollingBallFilter)

#INCLUDE_DIRECTORIES()

SET(RollingBallFilter_HDRS 
	RollingBallFilter.h
)

SET(RollingBallFilter_SRCS 
	RollingBallFilter.cpp
)

ADD_LIBRARY(RollingBallFilter ${RollingBallFilter_HDRS} ${RollingBallFilter_SRCS} )
ADD_EXECUTABLE(RollingBallFilterTester main.cpp)

TARGET_LINK_LIBRARIES(RollingBallFilterTester RollingBallFilter ${ITK_LIBRARIES})
INSTALL( TARGETS RollingBallFilterTester RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
