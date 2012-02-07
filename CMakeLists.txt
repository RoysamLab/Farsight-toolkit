CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( Farsight )

IF ("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")
  MESSAGE(FATAL_ERROR "Farsight must be built with separate source and binary directories.")
ENDIF ("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")
#######################################################################################
#########################################
# ADD NEW MODULES HERE
#########################################
OPTION( BUILD_TRACING "Build Tracing Tools" ON )
OPTION( BUILD_NUCLEI "Build Nuclear Segmentation" ON )
OPTION( BUILD_VESSEL "Build Vessel Surface Segmentation" OFF )
OPTION( BUILD_REGISTRATION "Build Image Registration (fregl)" OFF )
OPTION( BUILD_SAMPLE_EDITOR "Build Sample Editor" OFF )
IF(BUILD_SAMPLE_EDITOR)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SampleEditor )
  SET( BUILD_SPD ON CACHE BOOL "Build SPD" FORCE )
  SET( BUILD_CLUSCLUS ON CACHE BOOL "Build CLUSCLUS" FORCE )
ENDIF()
#OPTION( BUILD_TISSUENETS "Build TissueNets" OFF )
OPTION( BUILD_TRACKER "Build Cell Tracker" OFF )
OPTION( USE_KPLS "Use KPLS module for classification" ON )
OPTION( USE_OPENSLIDE "Use openslide libraries for reading images" OFF )
OPTION( WRAP_PYTHON "Wrap Farsight in Python" OFF )
OPTION( BUILD_EMBREX "Build Mark Embrecht's KPLS module" OFF )
OPTION( BUILD_UNMIXING "Build Unmixing module" OFF )
OPTION( BUILD_OPENMP "Build with OpenMP support (Note that only VS2008 Pro and VS2008 Team System supports OpenMP)" OFF)
OPTION( BUILD_POWEROF2TILING "Build PowerOf2Tiling.exe which will take an image and split it into tiles of 1024 pixels, zero-padding where needed" OFF)
OPTION( BUILD_ROLLINGBALLFILTER "Build Rolling Ball Filter (aka Subtract Background)" OFF)
OPTION( BUILD_IMAGE_MANAGER  "Manage Image File Conversions and Directories" OFF)
OPTION( USE_TRACKING "Use Tracking module with nucleus editor" OFF)
OPTION( BUILD_VESSEL_TRACING "Build Ball Tracing for Vessels" ON )
OPTION( BUILD_MODEL_SEG "Build Statistical_Model_Segmentation (Note that it has dependency on glpk library)" OFF)
OPTION( BUILD_ZERNIKE  " Build zernike_moments_calculation " OFF)
OPTION( BUILD_HARALICK "Biild Haralick_features_calculation" OFF)
OPTION( BUILD_DENDROGRAM "Dendrogram code" OFF)
OPTION( BUILD_SPD "Build SPD" OFF)
IF( BUILD_SPD )
	ADD_DEFINITIONS( -DUSE_SPD )
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SPD )
ENDIF()
OPTION( BUILD_XML_GENERATOR "Build XML_Generator" OFF)
OPTION( BUILD_CLUSCLUS "Build CLUSCLUS" OFF)
IF( BUILD_CLUSCLUS )
  ADD_DEFINITIONS( -DUSE_Clusclus )
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ClusClus )
ENDIF()
OPTION( BUILD_IMONTAGE "Build iMontage" OFF)
OPTION( BUILD_FTKVOTING "Build Radial Voting" OFF)
OPTION( BUILD_FTKNUCLEARSEGMENTATIONNIC "Build Nicolas Segmentation" OFF)
OPTION( BUILD_Batch_Distance2Device "Build Distance to Device Batch Processing Code" OFF)
OPTION( BUILD_image_dicer "Build program to dice an image based on a text file centered around a coordinate" OFF)
OPTION( BUILD_CURVELETS "Build the curvelts based preprocessing executables" OFF)

#Make sure OPENCL AND CUDA ARE NOT BOTH ENABLED
OPTION( BUILD_OPENCL "Build with OpenCL support (VS 2008 only at the moment)" OFF) 
OPTION(	BUILD_CUDA "Build with CUDA support" OFF)

MARK_AS_ADVANCED(BUILD_EMBREX)
MARK_AS_ADVANCED(WRAP_PYTHON)
MARK_AS_ADVANCED(BUILD_ZERNIKE)
MARK_AS_ADVANCED(BUILD_CURVELETS)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

if(APPLE)
  set(EXE_TYPE MACOSX_BUNDLE)
  set(MAKE_BUNDLE MACOSX_BUNDLE)
  set(MACOSX_BUNDLE_BUNDLE_NAME "${PROJECT_NAME}")

  set(FARSIGHT_INSTALL_DESTINATION_ARGS
    RUNTIME DESTINATION "${FARSIGHT_BUNDLE_LOCATION}")

  set(FARSIGHT_INSTALL_DESTINATION_ARGS
    BUNDLE DESTINATION "${FARSIGHT_BUNDLE_LOCATION}")

  set(INSTALL_BIN_DIR ${FARSIGHT_BUNDLE_LOCATION}/bin)

else(APPLE)
  set(INSTALL_BIN_DIR bin)
endif(APPLE)

#####################################################################################
# CPACK SETUP FOR BUILDING INSTALLERS
######################################################################################
SET(CPACK_PACKAGE_VERSION_MAJOR 0)
SET(CPACK_PACKAGE_VERSION_MINOR 4)
SET(CPACK_PACKAGE_VERSION_PATCH 4)
SET(CPACK_PACKAGE_VENDOR "Bio-Image Analytics Laboratory")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${Farsight_SOURCE_DIR}/farsight-description.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolkit of image analysis methods")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Farsight ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_MAJOR=${CPACK_PACKAGE_VERSION_MAJOR})
ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_MINOR=${CPACK_PACKAGE_VERSION_MINOR})
ADD_DEFINITIONS(-DCPACK_PACKAGE_VERSION_PATCH=${CPACK_PACKAGE_VERSION_PATCH})

#first screen.  make longer
SET(CPACK_RESOURCE_FILE_WELCOME "${Farsight_SOURCE_DIR}/welcome.txt")
#second screen?
SET(CPACK_RESOURCE_FILE_README "${Farsight_SOURCE_DIR}/readme.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${Farsight_SOURCE_DIR}/license.txt") 

SET(CPACK_NSIS_MENU_LINKS "http://www.farsight-toolkit.org" "Farsight Online Help")

SET(CPACK_PACKAGING_INSTALL_PREFIX "/")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/Applications")
INCLUDE(CPack)

CONFIGURE_FILE (
  ${PROJECT_SOURCE_DIR}/FarsightConfig.h.in
  ${PROJECT_BINARY_DIR}/FarsightConfig.h
)

SET (EXE_DIR ${Farsight_BINARY_DIR}/exe CACHE PATH
  "Single Directory for all executables")

SET (LIB_DIR ${Farsight_BINARY_DIR}/libs CACHE PATH
  "Single Directory for all static Libraries")
  
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXE_DIR})
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR})
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})

#bundle location is used to make drag-n-droppable Mac OS X applications
set(FARSIGHT_BUNDLE_LOCATION 
  "FARSIGHT ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")


################################################################################
# FARSIGHT DEPENDS UPON A NUMBER OF OUTSIDE LIBRARIES
# THESE ARE ITK, VTK, VXL, QT, and CUDA
# WE MAKE SURE THEY ARE PRESENT HERE
################################################################################

#Find ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${ITK_USE_FILE} )

#Find VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${VTK_USE_FILE} )

#Find VXL
FIND_PACKAGE(VXL REQUIRED)
SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

FIND_PACKAGE(Qt4)
IF(QT4_FOUND)
    INCLUDE( ${QT_USE_FILE} )
ELSE(QT4_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without Qt4")
ENDIF(QT4_FOUND)

FIND_PACKAGE(Boost)

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS})
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})

IF(USE_OPENSLIDE)
        FIND_PACKAGE(OpenSlide)
        IF(OPENSLIDE_FOUND)
		ADD_DEFINITIONS(-DUSE_OPENSLIDE)
                LINK_DIRECTORIES( ${OPENSLIDE_BINARY_DIR} )
                INCLUDE_DIRECTORIES( ${OPENSLIDE_INCLUDE_DIR} )
        ELSE(OPENSLIDE_FOUND)
                MESSAGE(FATAL_ERROR "OpenSlide not found. Please set the 
                directories to OPENSLIDE_LIBRARY and OPENSLIDE_INCLUDE_DIR or
                turn off USE_OPENSLIDE")
        ENDIF(OPENSLIDE_FOUND)
ELSE(USE_OPENSLIDE)
        SET(OPENSLIDE_LIBRARY )
        SET(OPENSLIDE_INCLUDE_DIR )
ENDIF(USE_OPENSLIDE)

#####################################################################################
# Optionally find Farsight's data directory, used for testing and installers
# You can check out the data directory from subversion like so:
# svn co https://farsight-svn.ee.uh.edu/repos/farsight/data 
#####################################################################################

FIND_PATH(FARSIGHT_DATA_ROOT FarsightData.readme 
  ${Farsight_SOURCE_DIR}/data 
  ${Farsight_SOURCE_DIR}/../data 
  ${Farsight_SOURCE_DIR}/../../data 
  $ENV{FARSIGHT_DATA_ROOT}
  DOC "The repository for data used for testing.  To obtain from SVN: \"svn co https://farsight-svn.ee.uh.edu/repos/farsight/data\"")
  

##########################################################################
# Testing stuff
###########################################################################
ENABLE_TESTING()
IF(BUILD_TESTING)
  CONFIGURE_FILE(${Farsight_SOURCE_DIR}/CMake/CTestCustom.ctest.in
    ${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest @ONLY)
  FILE(WRITE ${Farsight_BINARY_DIR}/CTestCustom.cmake
    "INCLUDE(\"${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest\")\n")
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)

  #Find resources needed to test Qt GUIs
  FIND_FILE(QtTestingConfig QtTestingConfig.cmake DOC "Path to QtTestingConfig.cmake.  This file can be found in the QtTesting binary dir")
  IF(QtTestingConfig)
    INCLUDE(${QtTestingConfig})
    INCLUDE_DIRECTORIES( ${QtTesting_INCLUDE_DIRS} )
    LINK_DIRECTORIES( ${QtTesting_LIBRARY_DIR} )
    SET(QtTestingFound TRUE)
  ENDIF()

  # CSV Testing Support.
  ADD_SUBDIRECTORY( ftkTesting )
ENDIF(BUILD_TESTING)

################################################################################
################################################################################
# COMPILER FLAGS AND OTHERS FOR PARALLELIZATION
#
################################################################################

IF (BUILD_OPENMP)
    IF (MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    
endif (MSVC)
ENDIF(BUILD_OPENMP)


IF (BUILD_OPENCL)
    IF (MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D OPENCL")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D OPENCL")
    
endif (MSVC)
ENDIF(BUILD_OPENCL)

IF (BUILD_CUDA)
    FIND_PACKAGE(CUDA)
        IF (CUDA_FOUND)
    	    IF (MSVC)
       	        #SET (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
	        SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D CUDA")
 	        SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D CUDA")
	        SET (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20)
	    INCLUDE_DIRECTORIES( $ENV{CUDA_INC_PATH} )
        ENDIF(MSVC)
    ENDIF (CUDA_FOUND)
ENDIF(BUILD_CUDA)


IF(USE_KPLS)
ADD_DEFINITIONS(-DUSE_KPLS)
ENDIF(USE_KPLS)
################################################################################
################################################################################
# FARSIGHT ALSO HAS A NUMBER OF COMMON LIBRARIES THAT ARE INTERNAL TO
# THE PROJECT:
################################################################################
set( common_libraries
  VesselTracer
  ItkVtkGlue
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
FOREACH( _lib ${common_libraries} )
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/${_lib} )
ENDFOREACH()
FOREACH( _lib ${common_libraries} )
  ADD_SUBDIRECTORY( ${_lib} )
ENDFOREACH()

IF(USE_KPLS)
IF( NOT BUILD_EMBREX )
  IF(WIN32)
    IF(CMAKE_CL_64)
      SET(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/win64")
    ELSE(CMAKE_CL_64)
      SET(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/win32")
    ENDIF(CMAKE_CL_64)
  ELSEIF(APPLE)
    SET(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/osx")
  ELSEIF(UNIX)
    SET(KPLS_LIB_DIR "${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex/unix")
  ENDIF(WIN32)
  FIND_LIBRARY(LIB_KPLS NAMES libKPLS HINTS "${KPLS_LIB_DIR}")
ENDIF( NOT BUILD_EMBREX )
ENDIF(USE_KPLS)

IF(BUILD_VESSEL)
  IF(WIN32)
	SET (GLUT_INCLUDE_DIR "freeglut/include")
	SET (GLUT_glut_LIBRARY ${CMAKE_SOURCE_DIR}/freeglut/lib/freeglut.lib)
	FIND_PACKAGE(GLUT)
  ENDIF( WIN32 )
ENDIF(BUILD_VESSEL)



IF(BUILD_MODEL_SEG)
  IF(WIN32)
    IF(CMAKE_CL_64)
      LINK_DIRECTORIES(GLPK_LIB_DIR ${PROJECT_SOURCE_DIR}/PatternAnalysis/glpk )
    ENDIF(CMAKE_CL_64)
  ENDIF(WIN32)
  FIND_LIBRARY(LIB_GLPK NAMES glpk HINTS "${GLPK_LIB_DIR}" )
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(FINDGLPK DEFAULT_MSG LIB_GLPK)
ENDIF(BUILD_MODEL_SEG)

IF (USE_TRACKING)
  ADD_DEFINITIONS(-DUSE_TRACKING)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/CellTrackerLib )
  ADD_SUBDIRECTORY( CellTrackerLib )
ENDIF (USE_TRACKING)

IF (WIN32)
    IF ( ${CMAKE_SIZEOF_VOID_P} EQUAL "4")
        INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/pthreads/win32/ )
        LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/pthreads/win32/ )
        SET( PTHREAD_LIBRARIES optimized pthreadVC2 debug pthreadVC2d)
    ELSE ( ${CMAKE_SIZEOF_VOID_P} EQUAL "4")
        INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/pthreads/win64/ )
        LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/pthreads/win64/ )
        SET( PTHREAD_LIBRARIES optimized pthreadVC2 debug pthreadVC2d)
    ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL "4")
ENDIF(WIN32)

####################################################################################
####################################################################################
# NOW WE CAN INCLUDE THE MODULES ACCORDING TO THE OPTIONS SELECTED ABOVE:
####################################################################################
#Note IF A USER WANTS TO ADD AN EXECUTABLE TO AN INSTALLER, ALL OPTIONS
#USED TO BUILD SAID EXECUTABLE HAVE TO BE TURNED ON 
IF(BUILD_TRACING)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Tracing )
  ADD_SUBDIRECTORY(Tracing)
  if( BUILD_TraceEdit )
    SET(CPACK_PACKAGE_EXECUTABLES "TraceEdit" "Trace Editor" ${CPACK_PACKAGE_EXECUTABLES})
  endif()
ENDIF(BUILD_TRACING)

IF(BUILD_TRACKER)
  ADD_SUBDIRECTORY(Tracking)
ENDIF(BUILD_TRACKER)

IF(BUILD_NUCLEI)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/NuclearSegmentation )
  ADD_SUBDIRECTORY( NuclearSegmentation )
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SQLite )
  ADD_SUBDIRECTORY( SQLite )
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SQLite/NESqlite )
  ADD_SUBDIRECTORY( SQLite/NESqlite )
  SET(CPACK_PACKAGE_EXECUTABLES "NucleusEditor" "Nucleus Editor" ${CPACK_PACKAGE_EXECUTABLES})
ENDIF(BUILD_NUCLEI)

IF(BUILD_VESSEL_TRACING)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/VesselTracer )
  SET(CPACK_PACKAGE_EXECUTABLES "ftkVesselTracer" "Vessel Tracer" ${CPACK_PACKAGE_EXECUTABLES})
ENDIF(BUILD_VESSEL_TRACING)

IF(BUILD_VESSEL)
 INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Vessel )
  ADD_SUBDIRECTORY( Vessel )
ENDIF(BUILD_VESSEL)

IF(BUILD_REGISTRATION)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/fregl )
  ADD_SUBDIRECTORY( fregl )
ENDIF(BUILD_REGISTRATION)

IF(BUILD_SAMPLE_EDITOR)
  ADD_SUBDIRECTORY( SampleEditor )
ENDIF(BUILD_SAMPLE_EDITOR)

#IF(BUILD_TISSUENETS)
#  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tissuenets )
#  ADD_SUBDIRECTORY(tissuenets)
#ENDIF(BUILD_TISSUENETS)

IF(BUILD_UNMIXING)
	ADD_SUBDIRECTORY(Unmixing)
ENDIF(BUILD_UNMIXING)

IF(BUILD_IMAGE_MANAGER)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkImageManager )
	ADD_SUBDIRECTORY( ftkImageManager )
ENDIF(BUILD_IMAGE_MANAGER)

IF (BUILD_ZERNIKE)
 	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Zernike)
	ADD_SUBDIRECTORY(Zernike)
	IF(MSVC)
        SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D _ZERNIKE")
 	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _ZERNIKE")
	ENDIF(MSVC)
ENDIF(BUILD_ZERNIKE)

IF (BUILD_HARALICK)
 	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Haralick)
	ADD_SUBDIRECTORY(Haralick)
ENDIF(BUILD_HARALICK)

IF (BUILD_DENDROGRAM)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Dendrogram )
	ADD_SUBDIRECTORY(Dendrogram)
ENDIF (BUILD_DENDROGRAM)


IF (BUILD_SPD)
	ADD_SUBDIRECTORY(SPD)
ENDIF (BUILD_SPD)

IF (BUILD_XML_GENERATOR)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Xml_Generator )
	ADD_SUBDIRECTORY( Xml_Generator )
ENDIF (BUILD_XML_GENERATOR)

IF (BUILD_CLUSCLUS)
    ADD_SUBDIRECTORY(ClusClus)
ENDIF (BUILD_CLUSCLUS)

IF (BUILD_IMONTAGE)
    INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/iMontage )
    ADD_SUBDIRECTORY(iMontage)
ENDIF (BUILD_IMONTAGE)

IF (BUILD_FTKVOTING)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkAdapSeedDetection )
	ADD_SUBDIRECTORY( ftkAdapSeedDetection )
ENDIF (BUILD_FTKVOTING)

IF (BUILD_FTKNUCLEARSEGMENTATIONNIC)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkNuclearSegmentationNic )
	ADD_SUBDIRECTORY( ftkNuclearSegmentationNic )
ENDIF (BUILD_FTKNUCLEARSEGMENTATIONNIC)

IF (BUILD_POWEROF2TILING)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/PowerOf2Tiling )
	ADD_SUBDIRECTORY( PowerOf2Tiling )
ENDIF (BUILD_POWEROF2TILING)

IF (BUILD_ROLLINGBALLFILTER)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/RollingBallFilter )
	ADD_SUBDIRECTORY( RollingBallFilter )
ENDIF (BUILD_ROLLINGBALLFILTER)

IF (BUILD_Batch_Distance2Device)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Batch_Distance2Device )
	ADD_SUBDIRECTORY( Batch_Distance2Device )
ENDIF (BUILD_Batch_Distance2Device)

IF (BUILD_image_dicer)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/image_dicer)
	ADD_SUBDIRECTORY( image_dicer)
#	SET( BUILD_CURVELETS ON CACHE BOOL "Build the curvelts based preprocessing executables" FORCE )
ENDIF(BUILD_image_dicer)

IF (BUILD_CURVELETS)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Tracing/Preprocessing/CurveletFiltering/fdct_wrapping_cpp/src )
  ADD_SUBDIRECTORY( Tracing/Preprocessing/CurveletFiltering/fdct_wrapping_cpp/src )
ENDIF (BUILD_CURVELETS)

#####################################################################################
# SOME INSTALL RULES
######################################################################################

#Qt libraries 
IF(WIN32)
  INSTALL(FILES ${QT_LIBRARY_DIR}/QtCore4.dll ${QT_LIBRARY_DIR}/QtGui4.dll DESTINATION bin)
ENDIF(WIN32)
#IF(APPLE)
#  INSTALL(DIRECTORY ${QT_QTCORE_LIBRARY_RELEASE} ${QT_QTGUI_LIBRARY_RELEASE} ${QT_QTSQL_LIBRARY_RELEASE} ${QT_QTNETWORK_LIBRARY_RELEASE} DESTINATION bin)
#ENDIF(APPLE)

#python scripts
if(APPLE)
  INSTALL(DIRECTORY python/
          DESTINATION ${FARSIGHT_BUNDLE_LOCATION}/python
          PATTERN "*.svn" EXCLUDE)
else(APPLE)
  INSTALL(DIRECTORY python/ DESTINATION python PATTERN "*.svn" EXCLUDE)
endif(APPLE)

#install rules for sample data and results are kept in a separate
#CMakeLists.txt file in the data directory, so we don't have to keep changing
#the top-level CMakeLists.txt file whenever you want to include a new sample
#data set in the installers
IF(FARSIGHT_DATA_ROOT)
	#  ADD_SUBDIRECTORY( ${FARSIGHT_DATA_ROOT} )
ENDIF(FARSIGHT_DATA_ROOT)

INCLUDE(InstallRequiredSystemLibraries)
# IMPORTANT:
# Keep this chunk of code here at the very bottom of the main
# CMakeLists.txt file so that the install rules that install the
# bundles and any files inside of them are already run before these
# script steps are run.
add_subdirectory(InstallRules)


if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D GLUT_BUILDING_LIB")
endif(MSVC)

#message(STATUS "module path: ${CMAKE_MODULE_PATH} ${CMAKE_ROOT}")
