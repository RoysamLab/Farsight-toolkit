cmake_minimum_required(VERSION 2.8.4)
include(ExternalProject)

set(base "${CMAKE_BINARY_DIR}/ExternalProjects")
set_property(DIRECTORY PROPERTY EP_BASE ${base})

set(install_dir "${base}/Install")

option(BUILD_SHARED_LIBS "Should Farsight be built with shared libraries? (Not possible on Windows)" OFF)
if(WIN32)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Farsight cannot built with shared libraries on Windows" FORCE)
  mark_as_advanced(BUILD_SHARED_LIBS)
endif()
set(shared ${BUILD_SHARED_LIBS}) # setting to use for BUILD_SHARED_LIBS on all subsequent projects
set(testing OFF) # setting to use for BUILD_TESTING on all subsequent projects

option(BUILD_FARSIGHT_TESTING "Download and enable the Qt GUI testing framework for use in FARSIGHT?" ON)

############################################################################
# Boost
#
ExternalProject_Add(Boost
  URL http://farsight-toolkit.org/support/boost_1_47_0.tar.gz
  URL_MD5 ff180a5276bec773a7625cac7e2288e8
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)  
set(boost "${base}/Source/Boost")

# Compute -G arg for configuring external projects with the same CMake generator:
#
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()

# Set the default build type---this will affect all libraries and
# applications
#
set(build_type "")
if(CMAKE_BUILD_TYPE)
  set(build_type "${CMAKE_BUILD_TYPE}")
endif()

option(DEBUG_LEAKS "SHOW VTK DEBUG LEAKS" ON)
if (DEBUG_LEAKS)
  set(vtk_debug_leaks ON)
endif()
set(mac_args)
if(APPLE)
  set(mac_args
    -DVTK_USE_CARBON:BOOL=ON
    -DVTK_USE_COCOA:BOOL=OFF
    -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=""
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
  #Mac OS 10.7 (Lion) apparently ships with a version of PNG that VXL doesn't like.
  #If we're building on this OS, force VXL to build its own version of PNG
  if(${CMAKE_SYSTEM} MATCHES "Darwin-11")
    set(png_arg "-DVXL_FORCE_V3P_PNG:BOOL=ON")
  endif()
endif()

############################################################################
# VXL
#
ExternalProject_Add(VXL
  SVN_REPOSITORY "http://svn.code.sf.net/p/vxl/svn/trunk"
  SVN_REVISION -r "36611"
  SVN_TRUST_CERT 1
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/VXL
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_BRL:BOOL=OFF
    -DBUILD_CONVERSIONS:BOOL=OFF
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_GEL:BOOL=OFF
    -DBUILD_OUL:BOOL=OFF
    -DBUILD_OXL:BOOL=OFF
    -DBUILD_PRIP:BOOL=OFF
    -DBUILD_TBL:BOOL=OFF
    -DBUILD_RPL:BOOL=ON
    -DBUILD_RPL_RGTL:BOOL=ON
    -DBUILD_RPL_RTVL:BOOL=ON
    -DBUILD_TESTING:BOOL=${testing}
    -DVNL_CONFIG_LEGACY_METHODS:BOOL=ON
    ${mac_args}
    ${png_arg}
    ${VXL_EXTRA_SUPERBUILD_CMAKE_ARGS}
  INSTALL_COMMAND ""
)
set(VXL_DIR ${base}/Build/VXL)

############################################################################
# Qt
#
find_package(Qt4)

if(QT_QMAKE_EXECUTABLE)
  set(USE_SYSTEM_QT_DEFAULT ON)
else()
  set(USE_SYSTEM_QT_DEFAULT OFF)
endif()

option(USE_SYSTEM_QT "Use the system Qt4" ${USE_SYSTEM_QT_DEFAULT})

set(Qt_Target "")

if(NOT USE_SYSTEM_QT)
  unset(QT_QMAKE_EXECUTABLE CACHE)
  include(${CMAKE_CURRENT_SOURCE_DIR}/BuildQt.cmake)
  ExternalProject_Get_Property(Qt binary_dir)
  set(QT_QMAKE_EXECUTABLE "${binary_dir}/bin/qmake${CMAKE_EXECUTABLE_SUFFIX}")
  set(Qt_Target "Qt")
endif()

############################################################################
# VTK
#

set(VTK_PATCH_COMMAND "")
if(APPLE)
  set(VTK_PATCH_COMMAND
    "patch"
    "${base}/Source/VTK/GUISupport/Qt/Chart/vtkQtBarChart.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/vtkQtBarChart.apple.patch"
  )
endif()
message(STATUS "check: ${VTK_PATCH_COMMAND}")
ExternalProject_Add(VTK
  URL http://www.vtk.org/files/release/5.10/vtk-5.10.1.tar.gz
  URL_MD5 264b0052e65bd6571a84727113508789
  PATCH_COMMAND ${VTK_PATCH_COMMAND}
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/VTK
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DVTK_DEBUG_LEAKS:BOOL=${vtk_debug_leaks}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${testing}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DVTK_USE_GUISUPPORT:BOOL=ON
    -DVTK_USE_QT:BOOL=ON
    -DVTK_USE_QTCHARTS:BOOL=ON
    -DCMAKE_SKIP_RPATH:BOOL=OFF
    -DVTK_QT_USE_WEBKIT:BOOL=OFF
    -DVTK_USE_BOOST:BOOL=ON
    -DBoost_INCLUDE_DIR:FILEPATH=${boost}
    ${mac_args}
    ${VTK_EXTRA_SUPERBUILD_CMAKE_ARGS}
  INSTALL_COMMAND ""
  DEPENDS
    ${Qt_Target}
)
set(VTK_DIR ${base}/Build/VTK)

# ITK v4 needs a short path on Windows
# We setup different (shorter) directories for ITK here if necessary
if(WIN32)
  set(ITK_BASE_DIR_DEFAULT "C:/ITKRoot")
  set(ITK_BASE_DIR "${ITK_BASE_DIR_DEFAULT}"
    CACHE PATH "Base of all SuperBuild-built ITK source/build trees.  If this path is too long, ITK will fail to build.")
  set(ITK_DOWNLOAD_DIR "${ITK_BASE_DIR}/Download")
  set(ITK_SOURCE_DIR "${ITK_BASE_DIR}/src")
  set(ITK_BINARY_DIR "${ITK_BASE_DIR}/bin")
else()
  set(ITK_DOWNLOAD_DIR "${base}/Download/ITK")
  set(ITK_SOURCE_DIR "${base}/Source/ITK")
  set(ITK_BINARY_DIR "${base}/Build/ITK")
endif()

############################################################################
# ITK
#
ExternalProject_Add(ITK
  URL http://dl.dropbox.com/u/26629462/InsightToolkit-4.3.1.tar.gz
  URL_MD5 de443086f1f5dd27c3639644cebe488c
  DOWNLOAD_DIR ${ITK_DOWNLOAD_DIR}
  SOURCE_DIR ${ITK_SOURCE_DIR}
  BINARY_DIR ${ITK_BINARY_DIR}
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/ITK
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${testing}
    -DITK_USE_REVIEW:BOOL=ON
    -DITK_USE_SYSTEM_VXL:BOOL=ON
    -DVXL_DIR:FILEPATH=${base}/Build/VXL
    -DModule_ITKVtkGlue:BOOL=ON
    -DVTK_DIR:FILEPATH=${base}/Build/VTK
    ${mac_args}
    ${ITK_EXTRA_SUPERBUILD_CMAKE_ARGS}
  INSTALL_COMMAND ""
  DEPENDS
    "VXL"
    "VTK"
)

#check if we can build vessel
FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  SET(BUILD_VESSEL ON CACHE BOOL "Build Vessel Surface Segmentation")
ELSE()
  SET(BUILD_VESSEL OFF CACHE BOOL "Build Vessel Surface Segmentation")
ENDIF()

option(USE_KPLS "Use KPLS module for classification" ON)
option(BUILD_CURVELETS "Build the Curvelets preprocessing module" ON)
option(BUILD_image_dicer "Build the image dice and trace module" ON)

############################################################################
# fftw
#
ExternalProject_Add(fftw-2.1.5
  URL http://dl.dropbox.com/u/26629462/fftw-2.1.5.tar.gz
  URL_MD5 1a89d7034071875ff703144fda121e9a
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/fftw-2.1.5
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
)

#Used in the fftw search path
SET( FFTW_INCLUDE_DIRSS "${install_dir}/fftw-2.1.5/include" )
SET( FFTW_LIB_SEARCH_PATHS "${install_dir}/fftw-2.1.5/lib" )

############################################################################
# Testing support for Farsight
#
set(testing_args)
set(testing_deps)
if(BUILD_FARSIGHT_TESTING)
  #search for Farsight data directory
  FIND_PATH(FARSIGHT_DATA_ROOT FarsightData.readme 
    ${FarsightSuperBuild_SOURCE_DIR}/../../data 
    $ENV{FARSIGHT_DATA_ROOT} )

  #check it out if not found
  if(NOT FARSIGHT_DATA_ROOT)
    ExternalProject_Add(FarsightData
      SVN_REPOSITORY "https://farsight-svn.ee.uh.edu/repos/farsight/data"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )  
    set(FARSIGHT_DATA_ROOT "${base}/Source/FarsightData" CACHE FILEPATH "Location of farsight data directory" FORCE)
    set(testing_deps FarsightData)
  endif()

  #download and build QtTesting framework
  ExternalProject_Add(QtTesting
    URL http://farsight-toolkit.org/support/QtTesting.tar.gz
    URL_MD5 610e89c2c77e76c41cd5607dfe742638 
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/QtTesting
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      ${mac_args}
    INSTALL_COMMAND ""
    DEPENDS
      ${Qt_Target}
  )
  set(testing_args
    -DBUILD_TESTING:BOOL=ON
    -DQtTestingConfig:FILEPATH=${base}/Build/QtTesting/QtTestingConfig.cmake
    -DFARSIGHT_DATA_ROOT:FILEPATH=${FARSIGHT_DATA_ROOT}
    )
  list(APPEND testing_deps QtTesting)
endif()

############################################################################
# Farsight 
#
option(BUILD_REGISTRATION "Build the registration/mosaicing utilities" ON)

ExternalProject_Add(Farsight
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/Farsight"
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/Farsight
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBoost_INCLUDE_DIR:FILEPATH=${boost}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DVTK_DIR:FILEPATH=${VTK_DIR}
    -DITK_DIR:FILEPATH=${ITK_BINARY_DIR}
    -DVXL_DIR:FILEPATH=${VXL_DIR}
    -DBUILD_VESSEL:BOOL=${BUILD_VESSEL}
    -DBUILD_CURVELETS:BOOL=${BUILD_CURVELETS}
    -DBUILD_REGISTRATION:BOOL=${BUILD_REGISTRATION}
    -DBUILD_image_dicer=${BUILD_image_dicer}
    -DUSE_KPLS:BOOL=${USE_KPLS}
    -DFFTW_INCLUDE_DIRS=${FFTW_INCLUDE_DIRSS}
    -DFFTW_LIB_SEARCH_PATH=${FFTW_LIB_SEARCH_PATHS}
    ${testing_args}
    ${mac_args}
  DEPENDS
    VXL
    ITK
    VTK
    fftw-2.1.5
    ${testing_deps}
)
