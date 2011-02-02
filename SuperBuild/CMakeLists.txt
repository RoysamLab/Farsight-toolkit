cmake_minimum_required(VERSION 2.8.2)

project(FarsightSuperBuild)

include(CTest)
include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")

mark_as_advanced(
  BUILD_TESTING
  CMAKE_INSTALL_PREFIX
  DART_TESTING_TIMEOUT
)
