CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
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
OPTION( BUILD_VESSEL "Build Vessel Surface Segmentation" ON )
OPTION( BUILD_REGISTRATION "Build Image Registration (fregl)" ON )
OPTION( BUILD_BIONET "Build Bionet" OFF )
OPTION( BUILD_SAMPLE_EDITOR "Build Sample Editor" OFF )
OPTION( BUILD_TISSUENETS "Build TissueNets" OFF )
OPTION( BUILD_TRACKER "Build Cell Tracker" OFF)
OPTION( WRAP_PYTHON "Wrap Farsight in Python" OFF )
MARK_AS_ADVANCED(WRAP_PYTHON)


INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/ParseParameters.py
#               ${PROJECT_BINARY_DIR}/Common/ParseParameters.py)

CONFIGURE_FILE (
  ${PROJECT_SOURCE_DIR}/FarsightConfig.h.in
  ${PROJECT_BINARY_DIR}/FarsightConfig.h
)
################################################################################
################################################################################
# FARSIGHT DEPENDS UPON A NUMBER OF OUTSIDE LIBRARIES
# THESE ARE ITK, VTK, VXL, and QT :::: WE MAKE SURE THEY ARE PRESENT HERE:
################################################################################
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE( ${ITK_USE_FILE} )
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
FIND_PACKAGE(VXL)
IF(VTK_FOUND AND VXL_FOUND)
    INCLUDE( ${VTK_USE_FILE} )
    SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
    INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ELSEIF(VTK_FOUND AND VXL_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
ENDIF(VTK_FOUND AND VXL_FOUND)

FIND_PACKAGE(Qt4)
IF(QT4_FOUND)
    INCLUDE( ${QT_USE_FILE} )
ELSE(QT4_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without Qt4")
ENDIF(QT4_FOUND)
################################################################################
################################################################################
# FARSIGHT ALSO HAS A NUMBER OF COMMON LIBRARIES THAT ARE INTERNAL TO
# THE PROJECT: 
################################################################################
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkImage )
SUBDIRS( ftkImage )

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SegmentationCommon )
SUBDIRS( SegmentationCommon )

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tinyxml )
SUBDIRS( tinyxml )

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/PatternAnalysis )
LINK_DIRECTORIES( ${PROJECT_SOURCE_DIR}/PatternAnalysis/embrex )
SUBDIRS( PatternAnalysis )

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkGUI )
SUBDIRS( ftkGUI )

####################################################################################
####################################################################################
# NOW WE CAN INCLUDE THE MODULES ACCORDING TO THE OPTIONS SELECTED ABOVE:
####################################################################################
SET(CPACK_NSIS_MENU_LINKS "http://www.farsight-toolkit.org" "Farsight Online Help")

IF(BUILD_TRACING)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Tracing )
  SUBDIRS(Tracing)
ENDIF(BUILD_TRACING)

IF(BUILD_TRACKER)
  SUBDIRS(Tracking)
ENDIF(BUILD_TRACKER)

IF(BUILD_NUCLEI)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/NuclearSegmentation )
  SUBDIRS( NuclearSegmentation )
  SET(CPACK_PACKAGE_EXECUTABLES "NucleusEditor;NucleusEditor")  
ENDIF(BUILD_NUCLEI)

IF(BUILD_VESSEL)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Vessel )
  SUBDIRS( Vessel )
ENDIF(BUILD_VESSEL)

IF(BUILD_REGISTRATION)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/fregl )
  SUBDIRS( fregl )
ENDIF(BUILD_REGISTRATION)

IF(BUILD_SAMPLE_EDITOR)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SampleEditor )
  SUBDIRS( SampleEditor )
ENDIF(BUILD_SAMPLE_EDITOR)

IF(BUILD_BIONET)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tissuenets )
  SUBDIRS(tissuenets)
ENDIF(BUILD_BIONET)

IF(BUILD_TISSUENETS)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tissuenets )
  SUBDIRS(tissuenets)
ENDIF(BUILD_TISSUENETS)

#####################################################################################
IF(WRAP_PYTHON)
  IF(NOT BUILD_SHARED_LIBS)
    MESSAGE(FATAL_ERROR "Wrapping requires a shared build, change BUILD_SHARED_LIBS to ON")
  ENDIF(NOT BUILD_SHARED_LIBS)
  FIND_PACKAGE(PythonLibs)
  IF(PYTHON_INCLUDE_PATH)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  ELSE(PYTHON_INCLUDE_PATH)
    MESSAGE(FATAL_ERROR "Python libraries must be installed for python wrapping.")
  ENDIF(PYTHON_INCLUDE_PATH)
  # CableSwig is included in the source distribution of farsight
  SET(CMAKE_MODULE_PATH ${Farsight_SOURCE_DIR}/PythonWrapping/CableSwig/SWIG/CMake)
  SET(CableSwig_DIR ${Farsight_BINARY_DIR}/PythonWrapping/CableSwig CACHE PATH "CableSwig_DIR: The directory containing CableSwigConfig.cmake.")
  SET(CableSwig_INSTALL_ROOT ${Farsight_INSTALL_LIB_DIR}/CSwig)
  INCLUDE(${CableSwig_DIR}/CableSwigConfig.cmake OPTIONAL)
  SUBDIRS(PythonWrapping/CableSwig)
ENDIF(WRAP_PYTHON)

#####################################################################################
#####################################################################################
# SOME INSTALL RULES
######################################################################################
IF(WIN32)
  INSTALL(FILES ${QT_LIBRARY_DIR}/QtCore4.dll ${QT_LIBRARY_DIR}/QtGui4.dll DESTINATION bin)
ENDIF(WIN32)
IF(APPLE)
  INSTALL(DIRECTORY ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTSQL_LIBRARY} ${QT_QTNETWORK_LIBRARY} DESTINATION bin)
ENDIF(APPLE)

INSTALL(DIRECTORY python/ DESTINATION python PATTERN "*.svn" EXCLUDE)

INCLUDE(InstallRequiredSystemLibraries)


#####################################################################################
#####################################################################################
# CPACK SETUP FOR BUILDING INSTALLERS
######################################################################################
SET(CPACK_PACKAGE_VERSION_MAJOR 0)
SET(CPACK_PACKAGE_VERSION_MINOR 1)
SET(CPACK_PACKAGE_VERSION_PATCH 0)
INCLUDE(CPack)

###########################################################################
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
ENDIF(BUILD_TESTING)


