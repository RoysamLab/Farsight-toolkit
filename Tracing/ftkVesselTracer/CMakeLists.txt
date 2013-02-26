SET(ftkVesselTracer_sources
	Common.cxx
	ftkVesselTracer.cxx
	)
	
SET(ftkVesselTracer_headers
	Common.h
	ftkVesselTracer.h
	)	

ADD_EXECUTABLE ( ftkVesselTracer ${ftkVesselTracer_sources} ${ftkVesselTracer_headers} Tracer_main.cxx)

TARGET_LINK_LIBRARIES( ftkVesselTracer ${VTK_LIBRARIES} ${ITK_LIBRARIES} vcl vnl vnl_algo vnl_io mbl)

ADD_LIBRARY(VesselTracer_lib ${ftkVesselTracer_sources} ${ftkVesselTracer_headers} Tracer_main.cxx)

TARGET_LINK_LIBRARIES(VesselTracer_lib ${ITK_LIBRARIES} ${VTK_LIBRARIES} vcl vnl vnl_algo vnl_io mbl)

