add_executable(MedianImageFilter MedianImageFilter.cxx)
target_link_libraries(MedianImageFilter ${ITK_LIBRARIES})
install(TARGETS MedianImageFilter RUNTIME DESTINATION ${INSTALL_BIN_DIR})

add_executable(RobustAutomaticThreshold RobustAutomaticThreshold.cxx)
target_link_libraries(RobustAutomaticThreshold ${ITK_LIBRARIES})
install(TARGETS RobustAutomaticThreshold RUNTIME DESTINATION ${INSTALL_BIN_DIR})

add_executable(Open Open.cxx)
target_link_libraries(Open ${ITK_LIBRARIES})
install(TARGETS Open RUNTIME DESTINATION ${INSTALL_BIN_DIR})

add_executable(IsolateCentroids IsolateCentroids.cxx)
target_link_libraries(IsolateCentroids ${ITK_LIBRARIES})
install(TARGETS IsolateCentroids RUNTIME DESTINATION ${INSTALL_BIN_DIR})

add_custom_command(
  TARGET IsolateCentroids
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
       -D PROJECT_SOURCE_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${EXE_DIR}
       -D CMAKE_RUNTIME_OUTPUT_DIRECTORY:STRING=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
       -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
       -P ${CMAKE_CURRENT_SOURCE_DIR}/ConfigureExtractSomas.cmake)
install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ExtractSomas.py
        DESTINATION ${INSTALL_BIN_DIR}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

if(BUILD_TESTING AND FARSIGHT_DATA_ROOT)
  find_package(PythonInterp)
  if(PYTHON_EXECUTABLE)
    add_test(ExtractSomas
      ${PYTHON_EXECUTABLE}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ExtractSomas.py
      ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
      ${Farsight_BINARY_DIR}/Testing/Temporary/C1-1unmixed-crop-somas.mhd
      0.3
      0.5)
    set_tests_properties(ExtractSomas PROPERTIES
      PASS_REGULAR_EXPRESSION "I count 1 distinct somas."
      TIMEOUT 900) #15 minutes
  endif()
endif()

