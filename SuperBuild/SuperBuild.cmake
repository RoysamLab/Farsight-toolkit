include(ExternalProject)

set(base "${CMAKE_BINARY_DIR}/ExternalProjects")
set_property(DIRECTORY PROPERTY EP_BASE ${base})

set(install_dir "${base}/Install")

set(shared ON) # setting to use for BUILD_SHARED_LIBS on all subsequent projects
set(testing OFF) # setting to use for BUILD_TESTING on all subsequent projects

find_package(Qt4)
if(NOT QT_QMAKE_EXECUTABLE)
  message(FATAL_ERROR "error: install Qt4 or set QT_QMAKE_EXECUTABLE")
endif()
set(qmake "${QT_QMAKE_EXECUTABLE}")

#TODO: add external project stuff to download Boost if its not found.
find_package(Boost)
if(NOT Boost_INCLUDE_DIR)
  message(FATAL_ERROR "error: install Boost or set Boost_INCLUDE_DIR")
endif()
set(boost "${Boost_INCLUDE_DIR}")

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

# Maybe make an option for this, but force developers to see their leaks for now...
set(vtk_debug_leaks ON)

set(mac_args)
if(APPLE)
  set(mac_args
    -DVTK_USE_CARBON:BOOL=ON
    -DVTK_USE_COCOA:BOOL=OFF
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
endif()

ExternalProject_Add(VXL
  SVN_REPOSITORY "https://vxl.svn.sourceforge.net/svnroot/vxl/trunk"
  SVN_REVISION -r "29473"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}/vxl
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
    -DBUILD_TESTING:BOOL=${testing}
    ${mac_args}
)
set(VXL_DIR ${base}/Build/VXL)

ExternalProject_Add(VTK
  GIT_REPOSITORY git://vtk.org/VTK.git
  GIT_TAG v5.6.0
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DVTK_DEBUG_LEAKS:BOOL=${vtk_debug_leaks}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${testing}
    -DDESIRED_QT_VERSION:STRING=4
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${qmake}
    -DVTK_USE_GUISUPPORT:BOOL=ON
    -DVTK_USE_QT:BOOL=ON
    -DVTK_USE_RPATH:BOOL=ON
    -DVTK_QT_USE_WEBKIT:BOOL=OFF
    -DBoost_INCLUDE_DIR:FILEPATH=${boost}
    ${mac_args}
)
set(VTK_DIR ${base}/Build/VTK)

ExternalProject_Add(ITK
  GIT_REPOSITORY git://itk.org/ITK.git
  GIT_TAG v3.20.0
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${testing}
    -DITK_USE_REVIEW:BOOL=ON
    -DITK_USE_SYSTEM_VXL:BOOL=ON
    -DVXL_DIR:FILEPATH=${base}/Build/VXL
    ${mac_args}
  DEPENDS
    "VXL"
)
set(ITK_DIR ${base}/Build/ITK)

set(EXE_DIR ${CMAKE_CURRENT_BINARY_DIR}/exe)
set(LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/libs)

#check if we can build vessel
FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  SET(BUILD_VESSEL ON CACHE BOOL "Build Vessel Surface Segmentation")
ELSE()
  SET(BUILD_VESSEL OFF CACHE BOOL "Build Vessel Surface Segmentation")
ENDIF()

ExternalProject_Add(Farsight
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/Farsight"
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBoost_INCLUDE_DIR:FILEPATH=${boost}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${qmake}
    -DVTK_DIR:FILEPATH=${VTK_DIR}
    -DITK_DIR:FILEPATH=${ITK_DIR}
    -DVXL_DIR:FILEPATH=${VXL_DIR}
    -DEXE_DIR:FILEPATH=${EXE_DIR}
    -DLIB_DIR:FILEPATH=${LIB_DIR}
    -DBUILD_VESSEL:BOOL=${BUILD_VESSEL}
    ${mac_args}
  DEPENDS
    VXL
    ITK
    VTK
)
