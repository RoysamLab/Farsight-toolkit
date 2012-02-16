cmake_minimum_required(VERSION 2.6)
project( MultipleNeuronTracer )

ADD_EXECUTABLE( MultipleNeuronTracer MNT_MAIN.cpp MultipleNeuronTracer.cpp MultipleNeuronTracer.h)
TARGET_LINK_LIBRARIES(MultipleNeuronTracer
  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )
  
ADD_LIBRARY( MultipleNeuronTracerLib MultipleNeuronTracer.cpp MultipleNeuronTracer.h)
TARGET_LINK_LIBRARIES(MultipleNeuronTracerLib
  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

if(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS MultipleNeuronTracer DESTINATION bin )
endif(APPLE)
