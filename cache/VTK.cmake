#
# VTK.cmake
#

# Initial settings for building VTK for use with FARSIGHT.

#set(VTK_USE_TK OFF CACHE BOOL "Set to OFF for FARSIGHT")
#set(VTK_WRAP_PYTHON ON CACHE BOOL "Set to ON for FARSIGHT")

set(BUILD_DOCUMENTATION OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(BUILD_EXAMPLES OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(BUILD_TESTING OFF CACHE BOOL "Set to OFF for FARSIGHT")
set(VTK_USE_BOOST ON CACHE BOOL "Set to ON for FARSIGHT")
set(VTK_USE_GUISUPPORT ON CACHE BOOL "Set to ON for FARSIGHT")
set(VTK_USE_QT ON CACHE BOOL "Set to ON for FARSIGHT")
set(VTK_USE_SYSTEM_JPEG ON CACHE BOOL "Set to ON for FARSIGHT")
set(VTK_USE_SYSTEM_TIFF ON CACHE BOOL "Set to ON for FARSIGHT")
set(VTK_USE_SYSTEM_ZLIB ON CACHE BOOL "Set to ON for FARSIGHT")

#set(JPEG_INCLUDE_DIR $ENV{VXL_DIR}/v3p/jpeg CACHE STRING
#  "Set to VXL jpeg source")
#set(JPEG_LIBRARY $ENV{BUILD_DIR}/vxl/lib/libjpeg.* CACHE STRING
#  "Set to VXL jpeg library")
#set(TIFF_INCLUDE_DIR $ENV{VXL_DIR}/v3p/tiff CACHE STRING
#  "Set to VXL tiff source")
#set(TIFF_LIBRARY $ENV{BUILD_DIR}/vxl/lib/libtiff.* CACHE STRING
#  "Set to VXL tiff library")
#set(ZLIB_INCLUDE_DIR $ENV{VXL_DIR}/v3p/zlib CACHE STRING
#  "Set to VXL zlib source")
#set(ZLIB_LIBRARY $ENV{BUILD_DIR}/vxl/lib/libz.* CACHE STRING
#  "Set to VXL zlib library")

set(Boost_INCLUDE_DIR $ENV{BOOST_DIR} CACHE STRING "Set to Boost include path")
set(DESIRED_QT_VERSION 4 CACHE STRING "Set to 4 for FARSIGHT")
set(QT_QMAKE_EXECUTABLE $ENV{QMAKE_DIR} CACHE STRING "Set to qmake directory")
set(VTK_USE_CARBON OFF CACHE BOOL "Set to OFF for FARSIGHT on Mac OS X")
set(VTK_USE_COCOA ON CACHE BOOL "Set to ON for FARSIGHT on Mac OS X")
set(VTK_USE_QVTK ON CACHE BOOL "Set to ON for FARSIGHT")

# force generation of 32-bit code on Mac OS X
#set(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "Do 32-bit on Mac OS X")
