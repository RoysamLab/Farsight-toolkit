add_executable(MedianImageFilter1 MedianImageFilter.cxx)
target_link_libraries(MedianImageFilter1 ${ITK_LIBRARIES})

add_executable(RobustAutomaticThreshold RobustAutomaticThreshold.cxx)
target_link_libraries(RobustAutomaticThreshold ${ITK_LIBRARIES})

add_executable(HoleFilling HoleFilling.cxx)
target_link_libraries(HoleFilling ${ITK_LIBRARIES})

add_executable(Open Open.cxx)
target_link_libraries(Open ${ITK_LIBRARIES})

add_executable(IsolateCentroids1 IsolateCentroids.cxx)
target_link_libraries(IsolateCentroids1 ${ITK_LIBRARIES})

add_custom_command(
  TARGET IsolateCentroids1
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
       -D PROJECT_SOURCE_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${EXE_DIR}
       -D CMAKE_RUNTIME_OUTPUT_DIRECTORY:STRING=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
       -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
       -P ${CMAKE_CURRENT_SOURCE_DIR}/ConfigureExtractSomas.cmake)
