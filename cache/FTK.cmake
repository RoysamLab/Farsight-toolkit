#
# FTK.cmake
#

# Initial settings for building FARSIGHT.

set(Boost_INCLUDE_DIR /usr/local/include/boost-1_38 CACHE STRING "Set to default Mac installation directory")
set(BUILD_NUCLEI ON CACHE BOOL "Set to ON for FARSIGHT")
set(ITK_DIR /usr/local/lib/InsightToolkit CACHE STRING "Set to default Mac installation directory")
set(QT_QMAKE_EXECUTABLE /usr/local/Trolltech/Qt-4.5.1/bin/qmake CACHE STRING "Set to default Mac installation directory")
set(VTK_DIR /usr/local/lib/vtk-5.4 CACHE STRING "Set to default Mac installation directory")
