cmake_minimum_required(VERSION 2.8.10)

if(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
endif(COMMAND CMAKE_POLICY)

#########################################################################################
# Set compiler when in OSX to use GCC instead of Clang which is the default in OSX10.8+. 
# THIS BLOCK OF CODE MUST HAPPEN BEFORE THE project() DEFINITION
#########################################################################################
if (APPLE)
  if (NOT APPLE_USE_CLANG)
    MESSAGE (AUTHOR_WARNING "Mac OS X detected, setting GCC to default compilers. Define APPLE_USE_CLANG:BOOL=ON BEFORE configuring for the first time to use Clang")
    SET(CMAKE_C_COMPILER "/usr/bin/gcc")
    SET(CMAKE_CXX_COMPILER "/usr/bin/g++")
  else (NOT APPLE_USE_CLANG)
    MESSAGE (AUTHOR_WARNING "Using clang on Apple platform")
    SET(CMAKE_C_COMPILER "/usr/bin/clang")
    SET(CMAKE_CXX_COMPILER "/usr/bin/clang++")
  endif (NOT APPLE_USE_CLANG)
else (APPLE)
  if (DEFINED APPLE_USE_CLANG)  
    MESSAGE (FATAL_ERROR "You appear to have defined APPLE_USE_CLANG but you are not on an Apple platform")
  endif (DEFINED APPLE_USE_CLANG)
endif (APPLE)

###################################
# set project's name
###################################
project( Farsight )

####################################################################
# Make sure that binaries aren't being built in the source directory
####################################################################
if("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")
  message(FATAL_ERROR "Farsight must be built with separate source and binary directories.")
endif("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")


#########################################
# ADD NEW MODULES HERE
#########################################
option( BUILD_TRACING "Build Tracing Tools" ON )
option( BUILD_NUCLEI "Build Nuclear Segmentation" ON )
option( BUILD_VESSEL "Build Vessel Surface Segmentation" OFF )
option( BUILD_REGISTRATION "Build Image Registration (fregl)" OFF )
option( BUILD_SAMPLE_EDITOR "Build Sample Editor" OFF )
#option( BUILD_TISSUENETS "Build TissueNets" OFF )
option( BUILD_TRACKER "Build Cell Tracker" OFF )
option( USE_KPLS "Use KPLS module for classification" ON )
option( WRAP_PYTHON "Wrap Farsight in Python" OFF )
option( BUILD_UNMIXING "Build Unmixing module" OFF )
option( BUILD_OPENMP "Build with OpenMP support (Note that only VS2008 Pro and VS2008 Team System supports OpenMP)" OFF)
option( BUILD_POWEROF2TILING "Build PowerOf2Tiling.exe which will take an image and split it into tiles of 1024 pixels, zero-padding where needed" OFF)
option( BUILD_ROLLINGBALLFILTER "Build Rolling Ball Filter (aka Subtract Background)" OFF)
option( BUILD_IMAGE_MANAGER  "Manage Image File Conversions and Directories" OFF)
option( USE_TRACKING "Use Tracking module with nucleus editor" OFF)
option( BUILD_VESSEL_TRACING "Build Ball Tracing for Vessels" OFF )
option( BUILD_MODEL_SEG "Build Statistical_Model_Segmentation (Note that it has dependency on glpk library)" OFF)
option( BUILD_ZERNIKE  " Build zernike_moments_calculation " OFF)
option( BUILD_HARALICK "Biild Haralick_features_calculation" OFF)
option( BUILD_DENDROGRAM "Dendrogram code" OFF)
option( BUILD_XML_GENERATOR "Build XML_Generator" OFF)
option( BUILD_CLUSCLUS "Build CLUSCLUS" ON)
option( BUILD_SPD "Build SPD" ON)
option( BUILD_ACTIVEVALIDATION "Build ActiveValidation" OFF)
option( BUILD_IMONTAGE "Build iMontage" OFF)
option( BUILD_FTKVOTING "Build Radial Voting" OFF)
option( BUILD_FTKNUCLEARSEGMENTATIONNIC "Build Nicolas Segmentation" OFF)
option( BUILD_Batch_Distance2Device "Build Distance to Device Batch Processing Code" OFF)
option( BUILD_PIPELINE "Build program to dice an image based on a text file centered around a coordinate" OFF)
option( BUILD_CURVELETS "Build the curvelts based preprocessing executables" OFF)

# Make sure OPENCL AND CUDA ARE NOT BOTH ENABLED
option( BUILD_OPENCL "Build with OpenCL support (VS 2008 only at the moment)" OFF) 
option(	BUILD_CUDA "Build with CUDA support" OFF)
if( BUILD_CUDA AND BUILD_OPENCL )
  message( SEND_ERROR "Both BUILD_CUDA and BUILD_OPENCL cannot be enabled at the
  same time." )
endif()

mark_as_advanced(WRAP_PYTHON)
mark_as_advanced(BUILD_ZERNIKE)
mark_as_advanced(BUILD_CURVELETS)
include_directories( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )


#########################################################
# Throw errors about dependencies of BUILD options
#########################################################

if( BUILD_SAMPLE_EDITOR )
  if( NOT BUILD_SPD )
    message(FATAL_ERROR "Sample Editor requires Sample Progression Discovery. Turn on BUILD_SPD")
  endif()

  if( NOT BUILD_CLUSCLUS )
    message(FATAL_ERROR "SampleEditor requires Clusclus. Turn on BUILD_Clusclus")
  endif()
endif()

if( BUILD_NUCLEI )
  if( NOT BUILD_CLUSCLUS )
    message(FATAL_ERROR "NucleusEditor requires Clusclus. Turn on BUILD_Clusclus")
  endif()
  set( BUILD_ACTIVEVALIDATION ON CACHE BOOL "Build ActiveValidation" FORCE )
endif()


#########################################################
# Define the compiler flags for the USE variables here
#########################################################

if( BUILD_CLUSCLUS )
  add_definitions( -DUSE_Clusclus )
  include_directories( ${PROJECT_SOURCE_DIR}/ClusClus )
endif()

option( BUILD_SPD "Build SPD" OFF)
if( BUILD_SPD )
  add_definitions( -DUSE_SPD )
  include_directories( ${PROJECT_SOURCE_DIR}/SPD )
endif()

if( BUILD_VESSEL_TRACING )
  add_definitions( -DUSE_BALL_TRACER )
  include_directories( ${PROJECT_SOURCE_DIR}/ftkVesselTracer )
endif()


#####################################################################################
# CPACK SETUP FOR BUILDING INSTALLERS
######################################################################################
set(CPACK_PACKAGE_VERSION_MAJOR 0)
set(CPACK_PACKAGE_VERSION_MINOR 4)
set(CPACK_PACKAGE_VERSION_PATCH 5)
set(CPACK_PACKAGE_VENDOR "Bio-Image Analytics Laboratory")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${Farsight_SOURCE_DIR}/farsight-description.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolkit of image analysis methods")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Farsight ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

add_definitions(-DCPACK_PACKAGE_VERSION_MAJOR=${CPACK_PACKAGE_VERSION_MAJOR})
add_definitions(-DCPACK_PACKAGE_VERSION_MINOR=${CPACK_PACKAGE_VERSION_MINOR})
add_definitions(-DCPACK_PACKAGE_VERSION_PATCH=${CPACK_PACKAGE_VERSION_PATCH})

#first screen.  make longer
set(CPACK_RESOURCE_FILE_WELCOME "${Farsight_SOURCE_DIR}/welcome.txt")
#second screen?
set(CPACK_RESOURCE_FILE_README "${Farsight_SOURCE_DIR}/readme.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${Farsight_SOURCE_DIR}/license.txt") 

set(CPACK_NSIS_MENU_LINKS "http://www.farsight-toolkit.org" "Farsight Online Help")

set(CPACK_PACKAGING_INSTALL_PREFIX "/")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/Applications")
include(CPack)

#bundle location is used to make drag-n-droppable Mac OS X applications
set(FARSIGHT_BUNDLE_LOCATION
  "FARSIGHT ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

if(APPLE)
  set(EXE_TYPE MACOSX_BUNDLE)
  set(MAKE_BUNDLE MACOSX_BUNDLE)
  set(MACOSX_BUNDLE_BUNDLE_NAME "${PROJECT_NAME}")

  set(FARSIGHT_INSTALL_DESTINATION_ARGS
    RUNTIME DESTINATION "${FARSIGHT_BUNDLE_LOCATION}")

  set(FARSIGHT_INSTALL_DESTINATION_ARGS
    BUNDLE DESTINATION "${FARSIGHT_BUNDLE_LOCATION}")

  set(INSTALL_BIN_DIR ${FARSIGHT_BUNDLE_LOCATION}/bin)
  set(INSTALL_LIB_DIR ${FARSIGHT_BUNDLE_LOCATION}/lib)

else(APPLE)
  set(INSTALL_BIN_DIR bin)
  set(INSTALL_LIB_DIR lib)
endif(APPLE)

configure_file(
  ${PROJECT_SOURCE_DIR}/FarsightConfig.h.in
  ${PROJECT_BINARY_DIR}/FarsightConfig.h
  )

set( EXE_DIR ${Farsight_BINARY_DIR}/exe CACHE PATH
  "Single Directory for all executables")
set( LIB_DIR ${Farsight_BINARY_DIR}/libs CACHE PATH
  "Single Directory for all static Libraries")

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXE_DIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR} )

# the add_farsight_executable function allows us to make redistributable
# executables on Linux
include(CMake/AddFarsightExecutable.cmake)

################################################################################
# FARSIGHT DEPENDS UPON A NUMBER OF OUTSIDE LIBRARIES
# THESE ARE ITK, VTK, VXL, QT, and CUDA
# WE MAKE SURE THEY ARE PRESENT HERE
################################################################################

#Find ITK
find_package(ITK REQUIRED)
include( ${ITK_USE_FILE} )

#Find VTK
find_package(VTK REQUIRED)
include( ${VTK_USE_FILE} )
if("${VTK_MAJOR_VERSION}" EQUAL 6)
  message("Debug: VTK MAJOR VER=${VTK_MAJOR_VERSION}")
  set(FARSIGHT_VTK_6 1)
else()
  message("Debug: VTK MAJOR VER=${VTK_MAJOR_VERSION}")
  set(FARSIGHT_VTK_6 0)
endif()

#Find VXL
find_package(VXL REQUIRED)
set(VXL_PROVIDE_STANDARD_optionS 1)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)
# UseITK.cmake utilizes include_directories(BEFORE..., so we do the same here
# because we want the VXL vxl_config.h to be found first.
include_directories(BEFORE ${VXL_VCL_INCLUDE_DIR} ${VXL_CORE_INCLUDE_DIR} ${VXL_RPL_INCLUDE_DIR})

# Find Qt
find_package(Qt4)
if(QT4_FOUND)
    include( ${QT_USE_FILE} )
else(QT4_FOUND)
    message(FATAL_ERROR "Cannot build without Qt4")
endif(QT4_FOUND)

# Find Boost
find_package(Boost)
include_directories ( ${Boost_INCLUDE_DIRS})

#####################################################################################
# Optionally find Farsight's data directory, used for testing and installers
# You can check out the data directory from subversion like so:
# svn co https://farsight-svn.ee.uh.edu/repos/farsight/data
#####################################################################################

find_path(FARSIGHT_DATA_ROOT FarsightData.readme
  ${Farsight_SOURCE_DIR}/data
  ${Farsight_SOURCE_DIR}/../data
  ${Farsight_SOURCE_DIR}/../../data
  $ENV{FARSIGHT_DATA_ROOT}
  DOC "The repository for data used for testing.  To obtain from SVN: \"svn co https://farsight-svn.ee.uh.edu/repos/farsight/data\""
  )


##########################################################################
# Testing stuff
###########################################################################
enable_testing()
include(CTest)
if(BUILD_TESTING)
  configure_file(${Farsight_SOURCE_DIR}/CMake/CTestCustom.ctest.in
    ${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest @ONLY)
  file(WRITE ${Farsight_BINARY_DIR}/CTestCustom.cmake
    "include(\"${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest\")\n")
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)

  #Find resources needed to test Qt GUIs
  find_file(QtTestingConfig QtTestingConfig.cmake DOC "Path to QtTestingConfig.cmake.  This file can be found in the QtTesting binary dir")
  if(QtTestingConfig)
    include(${QtTestingConfig})
    include_directories( ${QtTesting_INCLUDE_DIRS} )
    link_directories( ${QtTesting_LIBRARY_DIR} )
    set(QtTestingFound TRUE)
  endif()

  # CSV Testing Support.
  add_subdirectory( ftkTesting )
endif(BUILD_TESTING)

################################################################################
################################################################################
# COMPILER FLAGS AND OTHERS FOR PARALLELIZATION
#
################################################################################

if(BUILD_OPENMP)
  if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif (MSVC)
endif(BUILD_OPENMP)

if(BUILD_OPENCL)
  if(MSVC)
    add_definitions( -DOPENCL )
  endif (MSVC)
endif(BUILD_OPENCL)

if(BUILD_CUDA)
  find_package(CUDA)
  if(CUDA_FOUND)
    if(MSVC)
      add_definitions( -DCUDA )
      set( CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20)
      include_directories( $ENV{CUDA_INC_PATH} )
    endif(MSVC)
  endif(CUDA_FOUND)
endif(BUILD_CUDA)


################################################################################
################################################################################
# FARSIGHT ALSO HAS A NUMBER OF COMMON LIBRARIES THAT ARE INTERNAL TO
# THE project, WE INCLUDE THEM HERE:
################################################################################

if(USE_KPLS)
  add_definitions(-DUSE_KPLS)
  if(WIN32)
    if(CMAKE_CL_64)
      set(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/win64")
    else(CMAKE_CL_64)
      set(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/win32")
    endif(CMAKE_CL_64)
    if(MSVC90)
      find_library(LIB_KPLS_RELEASE NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Release/VS9)
      find_library(LIB_KPLS_DEBUG NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Debug/VS9)
    elseif(MSVC10)
      find_library(LIB_KPLS_RELEASE NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Release/VS10)
      find_library(LIB_KPLS_DEBUG NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Debug/VS10)
    elseif(MSVC11)
      find_library(LIB_KPLS_RELEASE NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Release/VS11)
      find_library(LIB_KPLS_DEBUG NAMES libKPLS PATHS ${KPLS_LIB_DIR}/Debug/VS11)
    else(MSVC90)
      message(FATAL_ERROR "Unsupported compiler for USE_KPLS ON.")
    endif(MSVC90)
    set(LIB_KPLS
      debug ${LIB_KPLS_DEBUG}
      optimized ${LIB_KPLS_RELEASE}
      )
  else(WIN32)
    if(UNIX)
      if(APPLE)
        set(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/osx")
        
      else(APPLE)
        set(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/unix")
      endif()
    endif()
      find_library(LIB_KPLS NAMES libKPLS HINTS "${KPLS_LIB_DIR}")
  endif(WIN32)
  if(MSVC10)
    configure_file(
      ${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/kpls_vs10.h
      ${PROJECT_BINARY_DIR}/PatternAnalysis/embrex/kpls.h
      COPYONLY)
  else(MSVC10)
    configure_file(
      ${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/kpls_default.h
      ${PROJECT_BINARY_DIR}/PatternAnalysis/embrex/kpls.h
      COPYONLY)
  endif(MSVC10)
endif(USE_KPLS)

if(BUILD_ZERNIKE)
  add_subdirectory( Zernike )
endif(BUILD_ZERNIKE)

set( common_libraries
  ftkImage
  ColorSegmentation
  ftkCommon
  ftkPreprocess
  ftkFeatures
  ftkSpectralUnmixing
  ftkGUI
  ftkGraphs
  tinyxml
  PatternAnalysis
  )
foreach( _lib ${common_libraries} )
  include_directories( ${PROJECT_SOURCE_DIR}/${_lib} )
endforeach()
foreach( _lib ${common_libraries} )
  add_subdirectory( ${_lib} )
endforeach()

if(BUILD_VESSEL)
  if(WIN32)
    set( GLUT_INCLUDE_DIR "freeglut/include")
    set( GLUT_glut_LIBRARY ${CMAKE_SOURCE_DIR}/freeglut/lib/freeglut.lib)
    find_package(GLUT)
  endif( WIN32 )
endif(BUILD_VESSEL)

if(BUILD_MODEL_SEG)
  if(WIN32)
    if(CMAKE_CL_64)
      link_directories(GLPK_LIB_DIR ${PROJECT_SOURCE_DIR}/PatternAnalysis/glpk )
    endif(CMAKE_CL_64)
  endif(WIN32)
  find_library(LIB_GLPK NAMES glpk HINTS "${GLPK_LIB_DIR}" )
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(FINDGLPK DEFAULT_MSG LIB_GLPK)
endif(BUILD_MODEL_SEG)

if(USE_TRACKING)
  add_subdirectory( CellTrackerLib )
endif(USE_TRACKING)

####################################################################################
####################################################################################
# NOW WE CAN INCLUDE THE MODULES ACCORDING TO THE OPTIONS SELECTED ABOVE:
####################################################################################
#Note IF A USER WANTS TO ADD AN EXECUTABLE TO AN INSTALLER, ALL OPTIONS
#USED TO BUILD SAID EXECUTABLE HAVE TO BE TURNED ON
if(BUILD_TRACING)
  add_subdirectory(Tracing)
  if( BUILD_TraceEdit )
    set(CPACK_PACKAGE_EXECUTABLES "TraceEdit" "Trace Editor" ${CPACK_PACKAGE_EXECUTABLES})
  endif()
endif(BUILD_TRACING)

if(BUILD_TRACKER)
  add_subdirectory(Tracking)
endif(BUILD_TRACKER)

if(BUILD_NUCLEI)
  add_subdirectory( SQLite )
  add_subdirectory( SQLite/NESqlite )
  add_subdirectory( NuclearSegmentation )
  set(CPACK_PACKAGE_EXECUTABLES "NucleusEditor" "Nucleus Editor" ${CPACK_PACKAGE_EXECUTABLES})
endif(BUILD_NUCLEI)

if(BUILD_VESSEL_TRACING)
  set(CPACK_PACKAGE_EXECUTABLES "ftkVesselTracer" "Vessel Tracer" ${CPACK_PACKAGE_EXECUTABLES})
  add_subdirectory( ftkVesselTracer )
endif(BUILD_VESSEL_TRACING)

#if(BUILD_SuperEllipsoidTracer)
 # 	set(CPACK_PACKAGE_EXECUTABLES "Trace3D" "Super" ${CPACK_PACKAGE_EXECUTABLES})
#endif(BUILD_SuperEllipsoidTracer)

if(BUILD_VESSEL)
 add_subdirectory( Vessel )
endif(BUILD_VESSEL)

if(BUILD_REGISTRATION)
  add_subdirectory( fregl )
endif(BUILD_REGISTRATION)

if(BUILD_SAMPLE_EDITOR)
  add_subdirectory( SampleEditor )
endif(BUILD_SAMPLE_EDITOR)

#if(BUILD_TISSUENETS)
#  include_directories( ${PROJECT_SOURCE_DIR}/tissuenets )
#  add_subdirectory(tissuenets)
#endif(BUILD_TISSUENETS)

if(BUILD_UNMIXING)
  add_subdirectory(Unmixing)
endif(BUILD_UNMIXING)

if(BUILD_IMAGE_MANAGER)
  add_subdirectory( ftkImageManager )
endif(BUILD_IMAGE_MANAGER)

if(BUILD_HARALICK)
  add_subdirectory(Haralick)
endif(BUILD_HARALICK)

if(BUILD_DENDROGRAM)
  add_subdirectory(Dendrogram)
endif(BUILD_DENDROGRAM)

if(BUILD_SPD)
  add_subdirectory(SPD)
endif(BUILD_SPD)

if(BUILD_XML_GENERATOR)
  add_subdirectory( Xml_Generator )
endif(BUILD_XML_GENERATOR)

if(BUILD_CLUSCLUS)
  add_subdirectory(ClusClus)
endif(BUILD_CLUSCLUS)

if(BUILD_ACTIVEVALIDATION)
  add_subdirectory(Active_Validation)
endif(BUILD_ACTIVEVALIDATION)

if(BUILD_IMONTAGE)
  add_subdirectory(iMontage)
endif(BUILD_IMONTAGE)

if(BUILD_FTKVOTING)
  add_subdirectory( ftkAdapSeedDetection )
endif(BUILD_FTKVOTING)

if(BUILD_FTKNUCLEARSEGMENTATIONNIC)
  add_subdirectory( ftkNuclearSegmentationNic )
endif(BUILD_FTKNUCLEARSEGMENTATIONNIC)

if(BUILD_POWEROF2TILING)
  add_subdirectory( PowerOf2Tiling )
endif(BUILD_POWEROF2TILING)

if(BUILD_ROLLINGBALLFILTER)
  add_subdirectory( RollingBallFilter )
endif(BUILD_ROLLINGBALLFILTER)

if(BUILD_Batch_Distance2Device)
  add_subdirectory( Batch_Distance2Device )
endif(BUILD_Batch_Distance2Device)

if(BUILD_PIPELINE)
  add_subdirectory( ftkPipeline)
#	set( BUILD_CURVELETS ON CACHE BOOL "Build the curvelts based preprocessing executables" FORCE )
endif(BUILD_PIPELINE)

if(BUILD_CURVELETS)
  add_subdirectory( Tracing/Preprocessing/CurveletFiltering/fdct_wrapping_cpp/src )
endif(BUILD_CURVELETS)

#####################################################################################
# SOME install RULES
######################################################################################

#python scripts
if(APPLE)
  install(DIRECTORY python/
          DESTINATION ${FARSIGHT_BUNDLE_LOCATION}/python
          PATTERN "*.svn" EXCLUDE)
else(APPLE)
  install(DIRECTORY python/ DESTINATION python PATTERN "*.svn" EXCLUDE)
endif(APPLE)

#install rules for sample data and results are kept in a separate
#CMakeLists.txt file in the data directory, so we don't have to keep changing
#the top-level CMakeLists.txt file whenever you want to include a new sample
#data set in the installers
if(FARSIGHT_DATA_ROOT)
	add_subdirectory( ${FARSIGHT_DATA_ROOT} ${Farsight_BINARY_DIR}/data )
endif(FARSIGHT_DATA_ROOT)

include(InstallRequiredSystemLibraries)

#find the OPENMP library directory for MSVC
if(MSVC60)
  set(MSVC_ver "60")
elseif(MSVC70)
  set(MSVC_ver "70")
elseif(MSVC71)
  set(MSVC_ver "71")
elseif(MSVC80)
  set(MSVC_ver "80")
elseif(MSVC90)
  set(MSVC_ver "90")
elseif(MSVC10)
  set(MSVC_ver "10")
endif()

if(MSVC_ver)
  set(MSVC_REDIST_DIR "${MSVC${MSVC_ver}_REDIST_DIR}")
  set(OPENMP_LIB_DIR "${MSVC_REDIST_DIR}/${CMAKE_MSVC_ARCH}/Microsoft.VC${MSVC_ver}.OPENMP")
  set(OPENMP_MANIFEST "${OPENMP_LIB_DIR}/Microsoft.VC${MSVC_ver}.OpenMP.manifest")
  install(FILES ${OPENMP_MANIFEST} DESTINATION ${INSTALL_BIN_DIR})
endif()

# IMPORTANT:
# Keep this chunk of code here at the very bottom of the main
# CMakeLists.txt file so that the install rules that install the
# bundles and any files inside of them are already run before these
# script steps are run.
add_subdirectory(InstallRules)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D GLUT_BUILDING_LIB")
endif(MSVC)
