cmake_minimum_required(VERSION 2.6)

project(ExtractSomas)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if(NOT ITK_USE_REVIEW)
  message(FATAL_ERROR "ExtractSomas requires ITK to be built with ITK_USE_REVIEW turned on.  Please reconfigure ITK if you'd like to use ExtractSoma")
endif()

if(APPLE)
  set(INSTALL_DEST ${FARSIGHT_BUNDLE_LOCATION}/bin)
else(APPLE)
  set(INSTALL_DEST bin)
endif(APPLE)


add_executable(MedianImageFilter MedianImageFilter.cxx)
target_link_libraries(MedianImageFilter ITKIO ITKStatistics)
install(TARGETS MedianImageFilter DESTINATION ${INSTALL_DEST} )

add_executable(RobustAutomaticThresholdImageFilterTest RobustAutomaticThresholdImageFilterTest.cxx)
target_link_libraries(RobustAutomaticThresholdImageFilterTest ITKIO ITKStatistics)
install(TARGETS RobustAutomaticThresholdImageFilterTest DESTINATION ${INSTALL_DEST} )

add_executable(PadImage PadImage.cxx)
target_link_libraries(PadImage ITKIO ITKStatistics)
install(TARGETS PadImage DESTINATION ${INSTALL_DEST} )

add_executable(HoleFillingFilterTest HoleFillingFilterTest.cxx)
target_link_libraries(HoleFillingFilterTest ITKIO ITKStatistics)
install(TARGETS HoleFillingFilterTest DESTINATION ${INSTALL_DEST} )

add_executable(ExtractImageFromPadding ExtractImageFromPadding.cxx)
target_link_libraries(ExtractImageFromPadding ITKIO ITKStatistics)
install(TARGETS ExtractImageFromPadding DESTINATION ${INSTALL_DEST} )

add_executable(DistanceMapFilter DistanceMapFilter.cxx)
target_link_libraries(DistanceMapFilter ITKIO ITKStatistics)
install(TARGETS DistanceMapFilter DESTINATION ${INSTALL_DEST} )

add_executable(BinaryThresholdImageFilter BinaryThresholdImageFilter.cxx)
target_link_libraries(BinaryThresholdImageFilter ITKIO ITKStatistics)
install(TARGETS BinaryThresholdImageFilter DESTINATION ${INSTALL_DEST} )

add_executable(PadDistanceMap PadDistanceMap.cxx)
target_link_libraries(PadDistanceMap ITKIO ITKStatistics)
install(TARGETS PadDistanceMap DESTINATION ${INSTALL_DEST} )

add_executable(DownwardFrontPropagation DownwardFrontPropagation.cxx)
target_link_libraries(DownwardFrontPropagation ITKIO ITKStatistics)
install(TARGETS DownwardFrontPropagation DESTINATION ${INSTALL_DEST} )

add_executable(IsolateCentroids IsolateCentroids.cxx)
target_link_libraries(IsolateCentroids ITKIO ITKStatistics)
install(TARGETS IsolateCentroids DESTINATION ${INSTALL_DEST} )

add_executable(VoronoiFilter VoronoiFilter.cxx)
target_link_libraries(VoronoiFilter ITKIO ITKStatistics)
install(TARGETS VoronoiFilter DESTINATION ${INSTALL_DEST} )

include(CTest)

add_custom_command(
  TARGET VoronoiFilter
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
       -D PROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${PROJECT_BINARY_DIR}
       -D CMAKE_RUNTIME_OUTPUT_DIRECTORY:STRING=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
       -P ${PROJECT_SOURCE_DIR}/ConfigureExtractSomas.cmake)
install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ExtractSomas.py DESTINATION ${INSTALL_DEST})

if(BUILD_TESTING AND FARSIGHT_DATA_ROOT)
  find_package(PythonInterp)
  if( PYTHON_EXECUTABLE)
    add_test(ExtractSomas
      ${PYTHON_EXECUTABLE}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ExtractSomas.py
      ${FARSIGHT_DATA_ROOT}/TraceData/C1-1unmixed-crop.tif
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/C1-1unmixed-crop-somas.tif)
    set_tests_properties(ExtractSomas PROPERTIES
      PASS_REGULAR_EXPRESSION "I count 1 distinct somas."
      TIMEOUT 900) #15 minutes
  endif()
endif()
