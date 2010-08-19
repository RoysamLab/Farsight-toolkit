#
# ITK.cmake
#

# Initial settings for building ITK for use with FARSIGHT.

set(BUILD_DOXYGEN OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(BUILD_EXAMPLES OFF CACHE BOOL "Set to OFF for FARSIGHT")
#set(BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Set to ON for FARSIGHT")
set(BUILD_TESTING OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(ITK_USE_SYSTEM_VXL ON CACHE BOOL "Set to ON for FARSIGHT")
set(VXL_DIR $ENV{BUILD_DIR}/vxl CACHE STRING "Set to VXL build directory")

#set(ITK_CSWIG_PYTHON ON CACHE BOOL "Set to ON for FARSIGHT")

# force generation of 32-bit code on Mac OS X
#set(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "Do 32-bit on Mac OS X")
