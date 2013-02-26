if(NOT BUILD_VESSEL_TRACING)
  message(FATAL_ERROR "AstroTracer requires ftkVesselTracer.  Turn on BUILD_VESSEL_TRACING")
endif()

SET(AstroTracer_sources
	AstroTracer.cxx
)

SET(AstroTracer_headers
	AstroTracer.h
)

ADD_LIBRARY( AstroTracerLib ${AstroTracer_sources} ${AstroTracer_headers} AstroTracer_main.cxx)
TARGET_LINK_LIBRARIES(AstroTracerLib ${ITK_LIBRARIES} ${VTK_LIBRARIES} vnl vnl_algo vnl_io MCLR ftkCommon VesselTracer_lib)

ADD_EXECUTABLE( AstroTracer ${AstroTracer_sources} ${AstroTracer_headers} AstroTracer_main.cxx)
TARGET_LINK_LIBRARIES(AstroTracer ${ITK_LIBRARIES} ${VTK_LIBRARIES}	vnl	vnl_algo vnl_io	MCLR ftkCommon VesselTracer_lib)

if(APPLE)
   install(TARGETS AstroTracer DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS AstroTracer DESTINATION bin )
endif(APPLE)
