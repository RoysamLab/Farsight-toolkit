project(image_dicer)

include_directories(${Farsight_SOURCE_DIR}/Tracing/TraceEdit)
include_directories(${Farsight_SOURCE_DIR}/NuclearSegmentation)

# set(image_dicer_HDRS 
#   )
# 
# set(image_dicer_SRCS 
# 	main.cpp
# 	)

SET(MAINDARPA_HDRS
	ftkMainDarpa.h
	ftkMainDarpa.hxx
	ftkMainDarpaTemplates.h
	ftkMainDarpaDeclaration.h

	ftkMainDarpaSegment.h
	ftkMainDarpaAstroTrace.h
	ftkMainDarpaAstroTrace.hxx
	ftkMainDarpaTrace.h
	ftkMainDarpaTrace.hxx
)

SET(MAINDARPA_SRCS
	main.cpp
	ftkMainDarpa.cpp

	ftkMainDarpaSegment.cpp
	ftkMainDarpaAstroTrace.cpp
	ftkMainDarpaTrace.cpp
)


#if(BUILD_OPENMP AND BUILD_CLUSCLUS AND BUILD_NUCLEI AND BUILD_TracingSystem)
if(BUILD_OPENMP AND BUILD_NUCLEI)

# 	add_executable(darpa_tracer_w_seg_v2 darpa_pipeline_w_seg_v2.cpp )
# 	target_link_libraries(darpa_tracer_w_seg_v2 vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)
# 
# 	add_executable(darpa_tracer_w_seg_v5 darpa_pipeline_w_seg_v5.cpp )
# 	target_link_libraries(darpa_tracer_w_seg_v5 vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)

	ADD_LIBRARY(LIBftkMainDarpa ${MAINDARPA_HDRS} ${MAINDARPA_SRCS})

	add_executable(	ftkMainDarpa main.cpp )
	target_link_libraries( ftkMainDarpa LIBftkMainDarpa ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib )

# 	add_executable(darpa_tracer_w_seg darpa_pipeline_w_seg.cpp )
# 	target_link_libraries(darpa_tracer_w_seg vul ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib segment_soma)
else()
  message(WARNING "BUILD_image_dicer is ON, but

  #BUILD_OPENMP or BUILD_CLUSCLUS or BUILD_NUCLEI or BUILD_TracingSystem
  BUILD_OPENMP or BUILD_NUCLEI
  
  are OFF.  image_dicer executables cannot be built without these modules")
endif()
