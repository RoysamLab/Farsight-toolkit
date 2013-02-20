SET(CTEST_SOURCE_NAME farsight-trunk)
SET(CTEST_BINARY_NAME farsight-dbg-static-nightly)
SET(CTEST_DASHBOARD_ROOT "c:/Dashboards")
SET(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

SET(CTEST_COMMAND
  "\"c:/Program Files (x86)/CMake 2.8/bin/ctest.exe\" -V -VV -D Nightly -A \"${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}\""
  )

SET(CTEST_CMAKE_COMMAND
  "\"c:/Program Files (x86)/CMake 2.8/bin/cmake.exe\""
  )

SET(CTEST_INITIAL_CACHE "
SITE:STRING=farsight-vista64
BUILDNAME:STRING=vs8-64-dbg-static
CMAKE_GENERATOR:INTERNAL=Visual Studio 8 2005 Win64
MAKECOMMAND:STRING=C:/PROGRA~2/MICROS~4/Common7/IDE/devenv.com Farsight.sln /build Debug /project ALL_BUILD
BUILD_SHARED_LIBS:BOOL=OFF
ITK_DIR:PATH=C:/Dashboards/Dependencies/ITK-3.16.0-static
VTK_DIR:PATH=C:/Dashboards/Dependencies/VTK-5.4.2-static
VXL_DIR:PATH=C:/Dashboards/Dependencies/VXL-SVN-static
Boost_INCLUDE_DIR:PATH=C:/Dashboards/Dependencies/boost_1_41_0
QT_QMAKE_EXECUTABLE:FILEPATH=C:/Dashboards/Dependencies/Qt-4.6.1/bin/qmake.exe
BUILD_VESSEL:BOOL=OFF
")

SET(CTEST_CVS_COMMAND "C:/Program Files/TortoiseSVN/bin/TortoiseProc.exe")

