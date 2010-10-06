cmake_minimum_required(VERSION 2.6)
project( MultipleNeuronTracer )

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR
  "Cannot build without ITK.  Please set ITK_DIR.")
endif(ITK_FOUND)

add_executable( MultipleNeuronTracer mainTracer3a.cxx )
target_link_libraries(MultipleNeuronTracer
  ITKNumerics ITKStatistics ITKIO ITKBasicFilters )

if(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION bin )
endif(APPLE)
