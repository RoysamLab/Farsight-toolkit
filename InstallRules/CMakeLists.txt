if(APPLE)

#list of directories where CMake can search for dependent libs
set(BUNDLE_LIB_DIRS
  ${LIB_DIR}
  ${VXL_DIR}/lib
)

#plugins loaded by TraceEdit at runtime
get_filename_component(qgif_plugin ${QT_QGIF_PLUGIN_RELEASE} NAME)
get_filename_component(qico_plugin ${QT_QICO_PLUGIN_RELEASE} NAME)
get_filename_component(qjpeg_plugin ${QT_QJPEG_PLUGIN_RELEASE} NAME)
get_filename_component(qmng_plugin ${QT_QMNG_PLUGIN_RELEASE} NAME)
get_filename_component(qsvg_plugin ${QT_QSVG_PLUGIN_RELEASE} NAME)
get_filename_component(qtiff_plugin ${QT_QTIFF_PLUGIN_RELEASE} NAME)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CompleteBundles.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake"
  @ONLY
)

install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake")

endif(APPLE)
