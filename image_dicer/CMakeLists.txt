project(image_dicer)

include_directories(${Farsight_SOURCE_DIR}/Tracing/TraceEdit)

set(image_dicer_HDRS 
  )

set(image_dicer_SRCS 
  main.cpp
  )

# ADD_EXECUTABLE(image_dicer main.cpp ${image_dicer_HDRS} ${image_dicer_SRCS})
# TARGET_LINK_LIBRARIES(image_dicer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg )
# INSTALL( TARGETS image_dicer RUNTIME DESTINATION ${INSTALL_BIN_DIR} )

# ADD_EXECUTABLE(darpa_tracer darpa_pipeline.cpp ${image_dicer_HDRS})
# TARGET_LINK_LIBRARIES(darpa_tracer vul ${ITK_LIBRARIES} Yousef_Nucleus_Seg MultipleNeuronTracerLib segment_soma)

if(BUILD_OPENMP AND BUILD_CLUSCLUS AND BUILD_NUCLEI AND BUILD_TracingSystem)
  add_executable(darpa_tracer_w_seg_v2 darpa_pipeline_w_seg_v2.cpp )
  target_link_libraries(darpa_tracer_w_seg_v2 vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)

  add_executable(darpa_tracer_w_seg darpa_pipeline_w_seg.cpp )
  target_link_libraries(darpa_tracer_w_seg vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)
else()
  message(WARNING "BUILD_image_dicer is ON, but

  BUILD_OPENMP or BUILD_CLUSCLUS or BUILD_NUCLEI or BUILD_TracingSystem
  
  are OFF.  image_dicer executables cannot be built without these modules")
endif()
