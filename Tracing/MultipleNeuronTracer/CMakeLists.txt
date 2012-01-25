cmake_minimum_required(VERSION 2.6)
project( MultipleNeuronTracer )

ADD_EXECUTABLE( MultipleNeuronTracer MNT_MAIN.cpp MultipleNeuronTracer.cpp MultipleNeuronTracer.h)
TARGET_LINK_LIBRARIES(MultipleNeuronTracer
  ${ITK_LIBRARIES} )
  
ADD_LIBRARY( MultipleNeuronTracerLib MultipleNeuronTracer.cpp MultipleNeuronTracer.h)
TARGET_LINK_LIBRARIES(MultipleNeuronTracerLib
  ${ITK_LIBRARIES} )

#MultipleNeuronTracer's test depends on ExtractSomas' test.
if(BUILD_ExtractSomas AND BUILD_TESTING AND FARSIGHT_DATA_ROOT AND PYTHON_EXECUTABLE)
  add_test(TraceSingleMicroglia
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/MultipleNeuronTracer
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/soma-centroids.txt
    1000
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/C1-1unmixed-crop-somas.tif
  )

  set_tests_properties(TraceSingleMicroglia PROPERTIES
    PASS_REGULAR_EXPRESSION "Total time to tracing"
    DEPENDS ExtractSomas
    TIMEOUT 600) #10 minutes
endif()

if(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION bin )
endif(APPLE)
