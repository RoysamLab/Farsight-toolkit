# Enable C and C++ languages to use things like MSVC_VERSION,
# CMAKE_COMPILER_IS_GNUCC and CMAKE_SIZEOF_VOID_P...
#
enable_language(C)
enable_language(CXX)

set(qt_version "4.7.2")
set(qt_url_md5 "66b992f5c21145df08c99d21847f4fdb")

set(qt_server "http://download.qt.nokia.com")
set(qt_baseurl "${qt_server}/qt/source")
set(qt_tgzname "qt-everywhere-opensource-src-${qt_version}.tar.gz")
set(qt_url "${qt_baseurl}/${qt_tgzname}")

set(suffix "unsupported")

if(MSVC10)
  set(suffix "vs2010")
elseif(MSVC90)
  set(suffix "vs2008")
elseif(MSVC80)
  set(suffix "vs2005")
elseif(MSVC71)
  set(suffix "vs2003")
elseif(CMAKE_COMPILER_IS_GNUCC)
  set(suffix "gcc")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(suffix "${suffix}-x64")
endif()

set(depends "")

set(use_jom 0)


# Base directory : determine default, allow user to choose with cache var:
#
if(WIN32)
  set(QT_BASE_DIR_DEFAULT "C:/Qt")
else()
  set(QT_BASE_DIR_DEFAULT "$ENV{HOME}/Qt")
endif()

set(QT_BASE_DIR "${QT_BASE_DIR_DEFAULT}"
  CACHE PATH "Base of all SuperBuild-built Qt source/build trees.")

set(qt_dir "${QT_BASE_DIR}/${qt_version}-${suffix}")


# Make command : use something that handles parallel builds for best
# Qt build time
#
if(WIN32)
  set(use_jom 1)
  set(JOM_EXECUTABLE "${QT_BASE_DIR}/SuperBuild/tools/jom.exe")
  set(make_cmd "${JOM_EXECUTABLE}")
else()
  set(QT_BASE_DIR_DEFAULT "$ENV{HOME}/Qt")
  set(make_cmd "make;-j9")
endif()


# jom, if needed:
#
if(use_jom AND NOT EXISTS "${JOM_EXECUTABLE}")
  ExternalProject_Add(DownloadJom
    URL "ftp://ftp.qt.nokia.com/jom/jom.zip"
    TIMEOUT 300
    PREFIX ${QT_BASE_DIR}/SuperBuild/${qt_version}-${suffix}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/jom.exe "${JOM_EXECUTABLE}"
    INSTALL_COMMAND ""
  )

  set(depends DEPENDS DownloadJom)
endif()

#a patch to fix an issue with Qt's timers on Linux
if(UNIX)
  set(patch_args PATCH_COMMAND patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/FixQtTimer.patch)
endif()

# Finally, Qt itself:
#
ExternalProject_Add(Qt
  URL ${qt_url}
  URL_MD5 ${qt_url_md5}
  TIMEOUT 2700 # 45 minutes
  PREFIX ${QT_BASE_DIR}/SuperBuild/${qt_version}-${suffix}
  DOWNLOAD_DIR ${QT_BASE_DIR}/Downloads
  SOURCE_DIR ${qt_dir}
  BUILD_IN_SOURCE 1
  ${patch_args}
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
    -D qt_dir=${qt_dir}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/ConfigureQt.cmake
  BUILD_COMMAND ${make_cmd}
  INSTALL_COMMAND ""
  ${depends}
)

message(STATUS "qt_dir='${qt_dir}'")
message(STATUS "qt_url='${qt_url}'")
message(STATUS "patch_args='${patch_args}'")

add_test(SmokeTest-Qt-${qt_version}-moc "${qt_dir}/bin/moc" "-v")
add_test(SmokeTest-Qt-${qt_version}-qmake "${qt_dir}/bin/qmake" "-v")
add_test(SmokeTest-Qt-${qt_version}-uic "${qt_dir}/bin/uic" "-v")

set_property(TEST SmokeTest-Qt-${qt_version}-moc PROPERTY
  PASS_REGULAR_EXPRESSION "Meta Object Compiler")
