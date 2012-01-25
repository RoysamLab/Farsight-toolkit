PROJECT(MicrogliaRegionTracer)

INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )

SET(MicrogliaRegionTracer_HDRS 
	MicrogliaRegionTracer.h
	Seed.h
)

SET(MicrogliaRegionTracer_SRCS 
	MicrogliaRegionTracer.cpp
	Seed.cpp
	main.cpp
)

ADD_EXECUTABLE(MicrogliaRegionTracer main.cpp ${MicrogliaRegionTracer_HDRS} ${MicrogliaRegionTracer_SRCS})
TARGET_LINK_LIBRARIES(MicrogliaRegionTracer fregl)
INSTALL( TARGETS MicrogliaRegionTracer RUNTIME DESTINATION bin )