PROJECT(image_dicer)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Tracing/TraceEdit)

SET(image_dicer_HDRS 
)

SET(image_dicer_SRCS 
	main.cpp
)

ADD_EXECUTABLE(image_dicer main.cpp ${image_dicer_HDRS} ${image_dicer_SRCS})
TARGET_LINK_LIBRARIES(image_dicer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg )
INSTALL( TARGETS image_dicer RUNTIME DESTINATION bin )

ADD_EXECUTABLE(darpa_tracer darpa_pipeline.cpp ${image_dicer_HDRS})
TARGET_LINK_LIBRARIES(darpa_tracer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg Multiple_Neuron_Tracer )
