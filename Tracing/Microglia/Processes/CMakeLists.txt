add_executable( MicrogliaProcessTracer
  MPT_MAIN.cpp MicrogliaProcessTracer.cpp MicrogliaProcessTracer.h)
target_link_libraries(MicrogliaProcessTracer ${ITK_LIBRARIES})
install(TARGETS MicrogliaProcessTracer RUNTIME DESTINATION ${INSTALL_BIN_DIR})

#MicrogliaProcessTracer's test depends on ExtractSomas' test.
if(BUILD_ExtractSomas AND BUILD_TESTING AND FARSIGHT_DATA_ROOT
                      AND PYTHON_EXECUTABLE)
  add_test(MicrogliaProcessTracer-SingleMicroglia
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/MicrogliaProcessTracer
    ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
    ${Farsight_BINARY_DIR}/Testing/Temporary/soma-centroids.txt
    ${Farsight_BINARY_DIR}/Testing/Temporary/C1-1unmixed-crop-somas.mhd
    3
    ${Farsight_BINARY_DIR}/Testing/Temporary/C1-1unmixed-crop-trace.swc
    0.3
    0.5
  )

  set_tests_properties(MicrogliaProcessTracer-SingleMicroglia PROPERTIES
    PASS_REGULAR_EXPRESSION "Total time to segmentation"
    DEPENDS ExtractSomas
    TIMEOUT 600) #10 minutes
endif()
