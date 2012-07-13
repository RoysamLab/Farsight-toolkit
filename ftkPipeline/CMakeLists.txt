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

if(BUILD_OPENMP AND BUILD_NUCLEI)

	ADD_LIBRARY(LIBftkMainDarpa ${MAINDARPA_HDRS} ${MAINDARPA_SRCS})

	add_executable(	ftkMainDarpa main.cpp )
	target_link_libraries( ftkMainDarpa LIBftkMainDarpa ${ITK_LIBRARIES} Project_Processor MultipleNeuronTracerLib AstroTracerLib)

else()
  message(WARNING "BUILD_image_dicer is ON, but

  BUILD_OPENMP or BUILD_NUCLEI
  
  are OFF.  image_dicer executables cannot be built without these modules")
endif()
