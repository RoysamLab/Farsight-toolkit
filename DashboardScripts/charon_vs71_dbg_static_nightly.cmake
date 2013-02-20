SET(CTEST_SOURCE_NAME farsight)
SET(CTEST_BINARY_NAME farsight-dbg-static-nightly)
SET(CTEST_DASHBOARD_ROOT "c:/projects")
SET(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

SET(CTEST_COMMAND
  "\"c:/Program Files/CMake 2.6/bin/ctest.exe\" -V -VV -D Nightly -A \"${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}\""
  )

SET(CTEST_CMAKE_COMMAND
  "\"c:/Program Files/CMake 2.6/bin/cmake.exe\""
  )

SET(CTEST_INITIAL_CACHE "
SITE:STRING=charon.kitware
BUILDNAME:STRING=vs71-32-dbg-static
CMAKE_GENERATOR:INTERNAL=Visual Studio 7 .NET 2003
MAKECOMMAND:STRING=\"c:/Program Files/Microsoft Visual Studio .NET 2003/Common7/IDE/devenv.com\" Farsight.sln /build Debug /project ALL_BUILD
BUILD_SHARED_LIBS:BOOL=OFF
ITK_DIR:PATH=c:/projects/ITK-3.12.0-shared
VTK_DIR:PATH=c:/projects/VTK-CVS-static
VXL_DIR:PATH=c:/projects/vxl-1.12.0-static
GLUT_INCLUDE_DIR:PATH=c:/projects/glut-3.7/include
GLUT_glut_LIBRARY:FILEPATH=c:/projects/glut-libs/glut32.lib
Boost_INCLUDE_DIR:PATH=c:/projects/boost_1_39_0
QT_QMAKE_EXECUTABLE:FILEPATH=c:/Qt/2009.02/qt/bin/qmake.exe
")

SET(CTEST_CVS_COMMAND "c:/Program Files/TortoiseSVN/bin/TortoiseProc.exe")

