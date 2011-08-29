SET(CTEST_SOURCE_NAME "src/farsight")
SET(CTEST_BINARY_NAME "bin/farsight-nightly")
SET(CTEST_DASHBOARD_ROOT "C:/dashboard")
SET(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

SET(CTEST_COMMAND
  "\"c:/Program Files (x86)/CMake 2.8/bin/ctest.exe\" -V -VV -D Nightly -A \"${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}\""
  )

SET(CTEST_CMAKE_COMMAND
  "\"c:/Program Files (x86)/CMake 2.8/bin/cmake.exe\""
  )

SET( CMPLR_PATH "C:/Program Files (x86)/Microsoft Visual Studio 9.0/Common7/IDE/devenv.com" )
SET( FPROJ_PATH "C:/dashboard/bin/farsight-nightly/Farsight.sln" )
SET( PRLL_STR_CXX "/DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /MP16" )
SET( PRLL_STR_C "/DWIN32 /D_WINDOWS /W3 /Zm1000 /MP16" )

SET(CTEST_INITIAL_CACHE "
SITE:STRING=farsight-win_7_64
BUILDNAME:STRING=vs9-64-dbg-nightly
CMAKE_GENERATOR:INTERNAL=Visual Studio 9 2008 Win64
MAKECOMMAND:STRING=${CMPLR_PATH} ${FPROJ_PATH} /build Debug /project ALL_BUILD
CMAKE_CXX_FLAGS:STRING=${PRLL_STR_CXX}
CMAKE_C_FLAGS:STRING=${PRLL_STR_C}
BUILD_SHARED_LIBS:BOOL=OFF
ITK_DIR:PATH=C:/dashboard/bin/itk-nightly
VTK_DIR:PATH=C:/dashboard/bin/vtk-nightly
VXL_DIR:PATH=C:/dashboard/bin/vxl-nightly
Boost_INCLUDE_DIR:PATH=C:/dashboard/src/boost
QT_QMAKE_EXECUTABLE:FILEPATH=C:/Qt/4.7.2/bin/qmake.exe
")

SET(CTEST_CVS_COMMAND "C:/Program Files/TortoiseSVN/bin/TortoiseProc.exe")
