function(add_farsight_executable exe_name srcs hdrs libs install_dir)

if(APPLE)
  set(EXE_TYPE MACOSX_BUNDLE)
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

add_executable(${exe_name}${FARSIGHT_EXE_SUFFIX} ${EXE_TYPE} ${hdrs} ${srcs})

target_link_libraries(${exe_name}${FARSIGHT_EXE_SUFFIX} ${libs})

install(TARGETS ${exe_name}${FARSIGHT_EXE_SUFFIX}
  RUNTIME DESTINATION ${install_dir}
  BUNDLE DESTINATION ${install_dir})

endfunction()
