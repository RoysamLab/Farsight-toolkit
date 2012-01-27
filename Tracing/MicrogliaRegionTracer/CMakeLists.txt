PROJECT(MicrogliaRegionTracer)

if(NOT BUILD_REGISTRATION)
  message(FATAL_ERROR "MicrogliaRegionTracer requires fregl.  Turn on BUILD_REGISTRATION")
endif()

INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )

SET(MicrogliaRegionTracer_HDRS 
	MicrogliaRegionTracer.h
	Seed.h
	ROIGrabber.h
	LoG.h
)

SET(MicrogliaRegionTracer_SRCS 
	MicrogliaRegionTracer.cpp
	Seed.cpp
	main.cpp
	ROIGrabber.cpp
	LoG.cpp
)

ADD_EXECUTABLE(MicrogliaRegionTracer main.cpp ${MicrogliaRegionTracer_HDRS} ${MicrogliaRegionTracer_SRCS})
TARGET_LINK_LIBRARIES(MicrogliaRegionTracer fregl)
INSTALL( TARGETS MicrogliaRegionTracer RUNTIME DESTINATION bin )
