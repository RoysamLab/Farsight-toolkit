PROJECT(RollingBallFilter)

#INCLUDE_DIRECTORIES()

SET(LibRollingBallFilter_HDRS 
	RollingBallFilter.h
)

SET(LibRollingBallFilter_SRCS 
	RollingBallFilter.cpp
)

ADD_LIBRARY(LibRollingBallFilter ${LibRollingBallFilter_HDRS} ${LibRollingBallFilter_SRCS} )
ADD_EXECUTABLE(RollingBallFilter main.cpp)

TARGET_LINK_LIBRARIES(RollingBallFilter LibRollingBallFilter ${ITK_LIBRARIES})
INSTALL( TARGETS RollingBallFilter RUNTIME DESTINATION bin )