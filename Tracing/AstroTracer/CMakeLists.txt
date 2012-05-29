cmake_minimum_required(VERSION 2.6)
PROJECT( AstroTracer )


FIND_PACKAGE(VTK)
FIND_PACKAGE(ITK)
FIND_PACKAGE(VXL)
IF(VTK_FOUND AND VXL_FOUND AND ITK_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
    INCLUDE( ${ITK_USE_FILE} )
    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ELSEIF(VTK_FOUND AND VXL_FOUND AND ITK_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build without VTK & ITK & VXL.")
ENDIF(VTK_FOUND AND VXL_FOUND AND ITK_FOUND)


SET(AstroTracer_sources
	AstroTracer.cxx
    itkHessianToObjectnessMeasureImageFilter.txx	
	itkMultiScaleHessianBasedMeasureImageFilter.txx
)

SET(AstroTracer_headers
	AstroTracer.h
	itkHessianToObjectnessMeasureImageFilter.h
	itkMultiScaleHessianBasedMeasureImageFilter.h
)


ADD_EXECUTABLE( AstroTracer ${AstroTracer_sources} ${AstroTracer_headers} AstroTracer_main.cxx)
TARGET_LINK_LIBRARIES(AstroTracer
  ${ITK_LIBRARIES} vnl vnl_algo vnl_io rtvl)



# MNT testing code
#MultipleNeuronTracer's test depends on ExtractSomas' test.
#if(BUILD_ExtractSomas AND BUILD_TESTING AND FARSIGHT_DATA_ROOT AND PYTHON_EXECUTABLE)
#  add_test(TraceSingleMicroglia
#    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/MultipleNeuronTracer
#    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
#    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/soma-centroids.txt
#    1000
#    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/C1-1unmixed-crop-somas.tif
#  )

#  set_tests_properties(TraceSingleMicroglia PROPERTIES
#    PASS_REGULAR_EXPRESSION "Total time to tracing"
#    DEPENDS ExtractSomas
#    TIMEOUT 600) #10 minutes
#endif()


if(APPLE)
   install(TARGETS AstroTracer DESTINATION ${FARSIGHT_BUNDLE_LOCATION})
else(APPLE)
   install(TARGETS AstroTracer DESTINATION bin )
endif(APPLE)
