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
    BoundingBox.h
    #BoundingVolumeHierarchyNode.h
    BoundingBoxFromTraceCalculator.h
	EigenAnalysis.h
	AspectRatioResampler.h
)

SET(MicrogliaRegionTracer_SRCS
	MicrogliaRegionTracer.cpp
	Cell.cpp
	main.cpp
	ROIGrabber.cpp
	LoG.cpp
	Node.cpp
	Tree.cpp
    BoundingBox.cpp
    #BoundingVolumeHierarchyNode.cpp
    BoundingBoxFromTraceCalculator.cpp
	EigenAnalysis.cpp
	AspectRatioResampler.txx
)

ADD_EXECUTABLE(MicrogliaRegionTracer main.cpp ${MicrogliaRegionTracer_HDRS} ${MicrogliaRegionTracer_SRCS})

TARGET_LINK_LIBRARIES(MicrogliaRegionTracer fregl)

INSTALL( TARGETS MicrogliaRegionTracer RUNTIME DESTINATION ${INSTALL_BIN_DIR} )

SET(MicrogliaRegionTracer_TEST_DATA_PATH "" CACHE
  FILEPATH "Path to the MicrogliaRegionTracer test data.")
SET(RUN_MicrogliaRegionTracer_TESTS OFF CACHE
  BOOL "Run the MicrogliaRegionTracer test.")
MARK_AS_ADVANCED(MicrogliaRegionTracer_TEST_DATA_PATH
  RUN_MicrogliaRegionTracer_TESTS)
IF(BUILD_TESTING AND RUN_MicrogliaRegionTracer_TESTS)
  IF(NOT EXISTS ${MicrogliaRegionTracer_TEST_DATA_PATH})
    MESSAGE(SEND_ERROR "The MicrogliaRegionTracer tests require the
    MicrogliaRegionTracer_TEST_DATA_PATH variable to be set.")
  ENDIF()
  ADD_TEST(MicrogliaRegionTracerTest ${EXE_DIR}/MicrogliaRegionTracer
    ${MicrogliaRegionTracer_TEST_DATA_PATH}/joint_transforms.xml
    ${MicrogliaRegionTracer_TEST_DATA_PATH}/
    8bitkt06045_w311GFPdsu.TIF
    ${MicrogliaRegionTracer_TEST_DATA_PATH}/seedpoints.txt
    )
ENDIF()
