# This function creates forwarding executables on Linux.  This allows us to
# create redistrutable executables that depend on shared libraries.  At the
# time of this writing, the primary source of shared libraries that Farsight
# depends on is Qt.  Therefore, you should use this function to create any
# executables that depend on Qt.

function(add_farsight_executable exe_name srcs libs)

if(APPLE)
  set(EXE_TYPE MACOSX_BUNDLE)
  set(install_dir ${FARSIGHT_BUNDLE_LOCATION})
else()
  set(install_dir bin)
endif()

#make forwarding executable on Linux
if(NOT WIN32 AND NOT APPLE)
  set(exe_output_path ${Farsight_BINARY_DIR}/exe)

  set(FARSIGHT_EXE_SUFFIX -real)
  set(FARSIGHT_FORWARD_DIR_BUILD "${exe_output_path}")
  set(FARSIGHT_FORWARD_DIR_INSTALL ".")
  set(FARSIGHT_FORWARD_PATH_BUILD "\"${FARSIGHT_FORWARD_DIR_BUILD}\"")
  set(FARSIGHT_FORWARD_PATH_INSTALL "\"${FARSIGHT_FORWARD_DIR_INSTALL}\"")
  set(FARSIGHT_FORWARD_EXE ${exe_name}-real)

  configure_file(
    ${Farsight_SOURCE_DIR}/CMake/farsight-forward.c.in
    ${CMAKE_CURRENT_BINARY_DIR}/${exe_name}-forward.c
    @ONLY IMMEDIATE)

  add_executable(${exe_name} ${CMAKE_CURRENT_BINARY_DIR}/${exe_name}-forward.c)

  install(TARGETS ${exe_name}
    RUNTIME DESTINATION ${install_dir}
    BUNDLE DESTINATION ${install_dir})
  add_dependencies(${exe_name} ${exe_name}${FARSIGHT_EXE_SUFFIX})
endif(NOT WIN32 AND NOT APPLE)

add_executable(${exe_name}${FARSIGHT_EXE_SUFFIX} ${EXE_TYPE} ${srcs})

target_link_libraries(${exe_name}${FARSIGHT_EXE_SUFFIX} ${libs})

install(TARGETS ${exe_name}${FARSIGHT_EXE_SUFFIX}
  RUNTIME DESTINATION ${install_dir}
  BUNDLE DESTINATION ${install_dir})

endfunction()
