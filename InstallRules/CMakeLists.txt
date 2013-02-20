#list of directories where CMake can search for dependent libs
set(BUNDLE_LIB_DIRS
  ${QT_LIBRARY_DIR}
  ${LIB_DIR}
  ${VXL_DIR}/lib
  ${ITK_DIR}/lib
  ${VTK_DIR}/bin
  ${QtTesting_LIBRARY_DIR}
  ${OPENMP_LIB_DIR}
)

#plugins loaded by TraceEdit at runtime
get_filename_component(qgif_plugin ${QT_QGIF_PLUGIN_RELEASE} NAME)
get_filename_component(qico_plugin ${QT_QICO_PLUGIN_RELEASE} NAME)
get_filename_component(qjpeg_plugin ${QT_QJPEG_PLUGIN_RELEASE} NAME)
get_filename_component(qmng_plugin ${QT_QMNG_PLUGIN_RELEASE} NAME)
get_filename_component(qsvg_plugin ${QT_QSVG_PLUGIN_RELEASE} NAME)
get_filename_component(qtiff_plugin ${QT_QTIFF_PLUGIN_RELEASE} NAME)

  
#specify where TraceEdit's required Qt plugins should be installed
if(APPLE)
  set(TraceEdit_plugins_dir "@FARSIGHT_BUNDLE_LOCATION@/TraceEdit.app/Contents/plugins/imageformats")
else()
  set(TraceEdit_plugins_dir "bin/plugins/imageformats")
endif()

set(suffix "${CMAKE_EXECUTABLE_SUFFIX}")
if(UNIX)
  set(suffix "-real")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CompleteBundles.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake"
  @ONLY
)

install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake")
