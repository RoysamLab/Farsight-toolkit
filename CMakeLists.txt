CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# set project's name
PROJECT( Farsight )

IF ("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")
  MESSAGE(FATAL_ERROR "Farsight must be built with separate source and binary directories.")
ENDIF ("${Farsight_BINARY_DIR}" MATCHES "^${Farsight_SOURCE_DIR}$")

#########################################
# ADD NEW MODULES HERE
# FOLLOW BUILD_YOUSEF EXAMPLE
#########################################

OPTION( BUILD_GUI "Build Farsight GUI" ON )
OPTION( BUILD_NUCLEI "Build Nuclear Segmentation" ON )
OPTION( BUILD_FUZZY_CLUSTERING "Build Nuclear Classification" ON )
OPTION( BUILD_REGISTRATION "Build Charlene's Executables and FREGL" ON )
OPTION( BUILD_MDL "Build Xiasong's MDL Module" ON )
OPTION( BUILD_TRACING "Build Yousef's Tracing Module" ON )
OPTION( WRAP_PYTHON "Wrap Farsight in Python" OFF)
OPTION( BUILD_VESSEL_SEGMENTATION "Build vessel surface segmentation" ON)
OPTION( BUILD_RENDERING "Build rendering module" ON)
OPTION( BUILD_TRACE_EDITOR "Build trace editor" ON)
OPTION( BUILD_TRACK_EDITOR "Build tracking editor" OFF)
OPTION( BUILD_BIONET "Build bionet" ON)
OPTION( BUILD_TRACE3D "Build Trace3D module" ON)
MARK_AS_ADVANCED(WRAP_PYTHON)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/Common/ParseParameters.py
#               ${PROJECT_BINARY_DIR}/Common/ParseParameters.py)

################################################################################
################################################################################
################################################################################
################################################################################
# MACROS (functions)
################################################################################
MACRO(MACRO_INCLUDE_COMMON_LIBS)

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
    #INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR} )
    #INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )
  ELSEIF(VTK_FOUND AND VXL_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without VTK & VXL.")
  ENDIF(VTK_FOUND AND VXL_FOUND)
  
ENDMACRO(MACRO_INCLUDE_COMMON_LIBS)
################################################################################
MACRO(MACRO_INCLUDE_QT4)
  FIND_PACKAGE(Qt4)
  IF(QT4_FOUND)
    INCLUDE( ${QT_USE_FILE} )
    IF (QT_QTXML_FOUND)
      INCLUDE_DIRECTORIES(${QT_QTXML_INCLUDE_DIR})
    ELSE(QT_QTXML_FOUND)
      MESSAGE(FATAL ERROR "Qt XML libraries missing.")
    ENDIF(QT_QTXML_FOUND)
    #IF (QT_QTHELP_FOUND)
    #  INCLUDE_DIRECTORIES(${QT_QTHELP_INCLUDE_DIR})
    #ELSE(QT_QTHELP_FOUND)
    #  MESSAGE(FATAL ERROR "Qt HELP libraries missing.")
    #ENDIF(QT_QTHELP_FOUND)
  ELSE(QT4_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without Qt4")
  ENDIF(QT4_FOUND)
ENDMACRO(MACRO_INCLUDE_QT4)
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
# BUILD THE PACKAGES THAT ARE SET TO BE BUILD
################################################################################
# These libraries are used in nearly all modules
################################################################################
MACRO_INCLUDE_COMMON_LIBS()
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/FTKImage )
SUBDIRS( FTKImage )
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/SegmentationCommon )
SUBDIRS( SegmentationCommon )
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/tinyxml )
SUBDIRS( tinyxml )
#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/libsvm )
SUBDIRS( libsvm )

#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/utils )
#SUBDIRS( utils )

#INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/ftkTest )
#SUBDIRS( ftkTest )

IF(BUILD_NUCLEI)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/NuclearSegmentation )
  SUBDIRS( NuclearSegmentation )  
ENDIF(BUILD_NUCLEI)
IF(BUILD_FUZZY_CLUSTERING)
  IF(BUILD_NUCLEI)
    INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/fuzzyClustering )
    SUBDIRS( fuzzyClustering)
  ELSE(BUILD_NUCLEI)
    MESSAGE(FATAL_ERROR "Cannot build fuzzy clustering without BUILD_NUCLEI set to ON ")
  ENDIF(BUILD_NUCLEI)
ENDIF(BUILD_FUZZY_CLUSTERING)


IF(BUILD_VESSEL_SEGMENTATION)

  FIND_PACKAGE(OpenGL)
  IF(OPENGL_FOUND)
	 INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIR} )
  ELSE(OPENGL_FOUND)
	MESSAGE(FATAL_ERROR "Cannot build without OpenGL")
  ENDIF(OPENGL_FOUND)

  FIND_PACKAGE(GLUT)
  IF(GLUT_FOUND)
	 INCLUDE_DIRECTORIES( ${GLUT_INCLUDE_DIR} )
  ELSE(GLUT_FOUND)
	MESSAGE(FATAL_ERROR "Cannot build without GLUT")
  ENDIF(GLUT_FOUND)

  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Vessel )
  SUBDIRS( Vessel )
  
ENDIF(BUILD_VESSEL_SEGMENTATION)
IF(BUILD_RENDERING)
SUBDIRS(Rendering)
ENDIF(BUILD_RENDERING)

IF(BUILD_MDL)
  FIND_PACKAGE( Boost REQUIRED )
#  INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/Modules/MDL )
#  SUBDIRS ( Modules/MDL )
ENDIF(BUILD_MDL) 

############################################################################
# ADD NEW MODULES ABOVE
############################################################################
#IF (BUILD_MODULE1)
#  INCLUDES: (Be smart so that previously included modules aren't included again)
#  SUBDIRS( Module1_Dir )
#ENDIF
#############################################################################
# ALSO NEED TO CHANGE FarsightConfig.h.in accordingly.
# ADD ifdef BUILD_MODULE1 statements to GUI source for automatically including
# your new module.
# AND will need to add TARGET_LINK_LIBRARIES( fsgui YOUR_MODULE ) in gui cmakelists.txt
#############################################################################


IF(BUILD_GUI)
  MACRO_INCLUDE_QT4()
  SUBDIRS( ftkGUI ) 
  #IF (BUILD_REGISTRATION)
  #   SUBDIRS( fregl_GUI )
  #ENDIF(BUILD_REGISTRATION)
  SET(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES}
    "Farsight" "Farsight")
#  SUBDIRS( Qt/Widgets )
#  SUBDIRS( Qt/Python )
ENDIF(BUILD_GUI)


IF(BUILD_REGISTRATION)
  IF( NOT ITK_USE_SYSTEM_VXL)
    MESSAGE(FATAL_ERROR "ITK must be built with ITK_USE_SYSTEM_VXL set ON.")
  ENDIF(NOT ITK_USE_SYSTEM_VXL)
  
  INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR} )
  INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )
  #SUBDIRS( Wrappers )
  #SUBDIRS( maciej_seg )
  #SUBDIRS( Common )
  #SUBDIRS( fsexe ) 
  SUBDIRS( fregl )
ENDIF(BUILD_REGISTRATION)

IF(BUILD_TRACING)
  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/Tracing )
  SUBDIRS(Tracing)
ENDIF(BUILD_TRACING)

#IF(BUILD_TRACE_EDITOR)
# SUBDIRS(TraceEdit)
#ENDIF(BUILD_TRACE_EDITOR)

#IF(BUILD_TRACK_EDITOR)
#  MACRO_INCLUDE_QT4()
#  SUBDIRS(TrackEdit)
#ENDIF(BUILD_TRACK_EDITOR)
IF(BUILD_BIONET)
  SUBDIRS(tissuenets)
ENDIF(BUILD_BIONET)

#IF(BUILD_TRACE3D)
#  SUBDIRS(ftl3d)
#ENDIF(BUILD_TRACE3D)

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

#These install rules are for use when linking against ITK/VTK build trees
#They won't work without modification if linked against ITK/VTK installs
SET(ITKBASE_SRC "")
SET(ITKBASE_SRC_NEW "")

#INSTALL(DIRECTORY ${VTK_DIR}/bin/\${BUILD_TYPE}/ DESTINATION .)
#INSTALL(DIRECTORY ${ITK_DIR}/bin/\${BUILD_TYPE}/ DESTINATION .)

#INSTALL(DIRECTORY ${ITK_DIR}/Wrapping/CSwig/Python/\${BUILD_TYPE}/ DESTINATION .)
#INSTALL(CODE "
#  file(READ \${CMAKE_INSTALL_PREFIX}/itkbase.py ITKBASE_SRC)
#  string(REGEX REPLACE \"pkgdir = '[^']*'\" \"pkgdir = '.'\" ITKBASE_SRC_NEW \"\${ITKBASE_SRC}\")
#  file(WRITE \${CMAKE_INSTALL_PREFIX}/itkbase.py \"\${ITKBASE_SRC_NEW}\")
#")

#INSTALL(DIRECTORY ${QT_LIBRARY_DIR}/ DESTINATION .)
IF(WIN32)
  INSTALL(FILES ${QT_LIBRARY_DIR}/QtCore4.dll ${QT_LIBRARY_DIR}/QtGui4.dll DESTINATION bin)
ENDIF(WIN32)
INSTALL(FILES ${QT_QTCORE_LIBRARY_RELEASE} ${QT_QTGUI_LIBRARY_RELEASE} DESTINATION bin)
INSTALL(DIRECTORY python/ DESTINATION python PATTERN "*.svn" EXCLUDE)
SET(CPACK_NSIS_MENU_LINKS "http://www.farsight-toolkit.org" "Farsight Online Help")
INCLUDE(CPack)

################################################################################
# FarsightConfig.h.in is modified and saved as FarsightConfig.h
# based on the configuration settings.
# This allows us to control what code is compiled based on the settings from
# the CMake GUI
################################################################################
CONFIGURE_FILE (
  ${PROJECT_SOURCE_DIR}/FarsightConfig.h.in
  ${PROJECT_BINARY_DIR}/FarsightConfig.h
)

##########################################################################
# Testing stuff
ENABLE_TESTING()
IF(BUILD_TESTING)
  CONFIGURE_FILE(${Farsight_SOURCE_DIR}/CMake/CTestCustom.ctest.in
    ${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest @ONLY)
  FILE(WRITE ${Farsight_BINARY_DIR}/CTestCustom.cmake
    "INCLUDE(\"${Farsight_BINARY_DIR}/CMake/CTestCustom.ctest\")\n")
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)


