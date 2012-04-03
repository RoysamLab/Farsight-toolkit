PROJECT(image_dicer)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Tracing/TraceEdit)

SET(image_dicer_HDRS 
)

SET(image_dicer_SRCS 
	main.cpp
)

# ADD_EXECUTABLE(image_dicer main.cpp ${image_dicer_HDRS} ${image_dicer_SRCS})
# TARGET_LINK_LIBRARIES(image_dicer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg )
# INSTALL( TARGETS image_dicer RUNTIME DESTINATION ${INSTALL_BIN_DIR} )

# ADD_EXECUTABLE(darpa_tracer darpa_pipeline.cpp ${image_dicer_HDRS})
# TARGET_LINK_LIBRARIES(darpa_tracer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg MultipleNeuronTracerLib segment_soma)

IF(BUILD_OPENMP AND BUILD_CLUSCLUS)
#   ADD_EXECUTABLE(darpa_tracer_w_seg darpa_pipeline_w_seg.cpp )
	ADD_EXECUTABLE(darpa_tracer_w_seg_v2 darpa_pipeline_w_seg_v2.cpp )
	TARGET_LINK_LIBRARIES(darpa_tracer_w_seg_v2 vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)

	ADD_EXECUTABLE(darpa_tracer_w_seg darpa_pipeline_w_seg.cpp )
	TARGET_LINK_LIBRARIES(darpa_tracer_w_seg vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)
ELSE()
  MESSAGE(WARNING "BUILD_image_dicer is ON, but either BUILD_OPENMP or BUILD_CLUSCLUS are OFF.  image_dicer executables cannot be built without these modules")
ENDIF()
