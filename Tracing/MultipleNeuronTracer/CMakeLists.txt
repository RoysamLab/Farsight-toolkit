cmake_minimum_required(VERSION 2.6)
project( MultipleNeuronTracer )

ADD_EXECUTABLE( MultipleNeuronTracer MNT_MAIN.cpp MultipleNeuronTracer.cpp PointOperation.cpp PointOperation.h	
	MultipleNeuronTracer.h )
TARGET_LINK_LIBRARIES(MultipleNeuronTracer
  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )
  
ADD_LIBRARY( MultipleNeuronTracerLib MultipleNeuronTracer.cpp MultipleNeuronTracer.h PointOperation.h PointOperation.cpp)
TARGET_LINK_LIBRARIES(MultipleNeuronTracerLib ftkCommon
  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

install(TARGETS MultipleNeuronTracer RUNTIME DESTINATION ${INSTALL_BIN_DIR})

