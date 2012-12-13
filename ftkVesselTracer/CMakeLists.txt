PROJECT (ftkVesselTracer)

SET(ftkVesselTracer_sources
	Common.cxx
	ftkVesselTracer.cxx
	itkHessianToObjectnessMeasureImageFilter.txx	
	itkMultiScaleHessianBasedMeasureImageFilter.txx
	)
	
SET(ftkVesselTracer_headers
	Common.h
	ftkVesselTracer.h
	itkHessianToObjectnessMeasureImageFilter.h
	itkMultiScaleHessianBasedMeasureImageFilter.h
	)	

ADD_EXECUTABLE ( ftkVesselTracer ${ftkVesselTracer_sources} ${ftkVesselTracer_headers} Tracer_main.cxx)

ADD_LIBRARY(VesselTracer_lib ${ftkVesselTracer_sources} ${ftkVesselTracer_headers} Tracer_main.cxx)

TARGET_LINK_LIBRARIES( ftkVesselTracer ${VTK_LIBRARIES} ${ITK_LIBRARIES} vcl vnl vnl_algo vnl_io rtvl)
