if (NOT BUILD_AstroTracer)
	message (FATAL_ERROR "BUILD_PIPELINE requires BUILD_AstroTracer")
endif()

if (NOT BUILD_MultipleNeuronTracer)
	message (FATAL_ERROR "BUILD_PIPELINE requires BUILD_MultipleNeuronTracer")
endif()

project(ftkPipeline)

include_directories(${Farsight_SOURCE_DIR}/Tracing/TraceEdit)
include_directories(${Farsight_SOURCE_DIR}/NuclearSegmentation)

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


if(BUILD_NUCLEI)

	ADD_LIBRARY(LIBftkMainDarpa ${MAINDARPA_HDRS} ${MAINDARPA_SRCS})

	add_executable(	ftkMainDarpa main.cpp )
	target_link_libraries( ftkMainDarpa LIBftkMainDarpa ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib AstroTracerLib)

else()
  message(FATAL_ERROR "BUILD_image_dicer is ON, but BUILD_NUCLEI is OFF.  image_dicer executables cannot be built without these modules")
endif()

add_executable( ftkSpatialStats ftkSpatialStats.cpp)
target_link_libraries( ftkSpatialStats ftkGraphs )



