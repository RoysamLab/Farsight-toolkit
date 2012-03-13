PROJECT(MicrogliaRegionTracer)

if(NOT BUILD_REGISTRATION)
  message(FATAL_ERROR "MicrogliaRegionTracer requires fregl.  Turn on BUILD_REGISTRATION")
endif()

INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )

SET(MicrogliaRegionTracer_HDRS
	MicrogliaRegionTracer.h
	Cell.h
	ROIGrabber.h
	LoG.h
	Node.h
	Tree.h
)

SET(MicrogliaRegionTracer_SRCS
	MicrogliaRegionTracer.cpp
	Cell.cpp
	main.cpp
	ROIGrabber.cpp
	LoG.cpp
	Node.cpp
	Tree.cpp
)

ADD_EXECUTABLE(MicrogliaRegionTracer main.cpp ${MicrogliaRegionTracer_HDRS} ${MicrogliaRegionTracer_SRCS})
ADD_EXECUTABLE(itkRolloverOpenMPTest itkRolloverOpenMPTest.cpp)

TARGET_LINK_LIBRARIES(MicrogliaRegionTracer fregl)
TARGET_LINK_LIBRARIES(itkRolloverOpenMPTest ${ITK_LIBRARIES})

ADD_TEST(itkRolloverOpenMPTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkRolloverOpenMPTest ${FARSIGHT_DATA_ROOT}/TraceData/C02.1344x1024x18.tif)
INSTALL( TARGETS MicrogliaRegionTracer RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
