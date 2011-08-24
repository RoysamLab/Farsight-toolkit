if(APPLE)
  #
  # Run the CompleteBundles.cmake script during "make install"
  #
  # This will pull all the required libraries into each bundle listed in
  # CompleteBundles.cmake.
  #

  #
  # WARNING:
  #   The file CompleteBundles.cmake.in is based on
  #   ParaView3/Applications/Client/CompleteBundle.cmake.in.
  #   If that file changes, this one will likely need a similar change.
  #

set(BUNDLE_LIB_DIRS
  ${LIB_DIR}
  ${VXL_DIR}/lib
)

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/CompleteBundles.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake"
    @ONLY
  )

  install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/CompleteBundles.cmake")

endif(APPLE)
