PROJECT(QtPython)

# Tell vtkWrapPython.cmake to set VTK_PYTHON_LIBRARIES for us.
SET(VTK_WRAP_PYTHON_FIND_LIBS 1)
INCLUDE("${VTK_CMAKE_DIR}/vtkWrapPython.cmake")

INCLUDE_DIRECTORIES(
  ${MPI_INCLUDE_PATH}
  ${PYTHON_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
#  "${ParaView_BINARY_DIR}/VTK"
#  "${ParaView_SOURCE_DIR}/Servers/Common"
#  "${ParaView_SOURCE_DIR}/Servers/Filters"
#  "${ParaView_SOURCE_DIR}/Servers/ServerManager"
#  "${ParaView_BINARY_DIR}/Servers/ServerManager"
#  "${ParaView_BINARY_DIR}/Servers/Filters"
#  "${ParaView_BINARY_DIR}/Servers/Common"
#  "${ParaView_SOURCE_DIR}/Utilities/VTKClientServer"
#  "${ParaView_BINARY_DIR}/Utilities/VTKClientServer"
  ${QtWidgets_SOURCE_DIR}
  ${QtWidgets_BINARY_DIR}
  ${QtPython_BINARY_DIR}
  ${VTK_INCLUDE_DIR}
  ${PYTHON_INCLUDE_PATH}
  ${PARAVIEW_INCLUDE_DIRS}
  ${VTKPythonWrapping_SOURCE_DIR}/Executable
  )

# create the VTK/Python  executable
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/vtkPythonAppInitConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkPythonAppInitConfigure.h)

ADD_LIBRARY(vtkPVPythonInterpretor
  vtkPVPythonInterpretor.cxx
  vtkPVPythonInteractiveInterpretor.cxx)

# When farsight is linked statically (for example during
# cross-compiling) the python wrapper module "init" functions must be
# referenced and called at startup.  The necessary code is written
# into pvpythonmodules.h and included in vtkPVPythonInterpretor.cxx
# for compilation.  Then the vtkPVPythonInterpretor static library
# "links" to all the module archives.  This creates a circular
# dependency chain between the interpreter and the wrapper modules.
#
# In the future we should change this to have a library whose purpose
# is to reference all the module init functions.  That library may
# then depend on the modules and vtkPVPythonInterpretor without
# creating a cycle.
IF(UNIX)
  SET(PYTHON_MODULE_PREFIX "lib")
ELSE(UNIX)
  SET(PYTHON_MODULE_PREFIX "")
ENDIF(UNIX)
PYTHON_WRITE_MODULES_HEADER(${CMAKE_CURRENT_BINARY_DIR}/qtpythonmodules.h)
IF(PY_STATIC_MODULES_LIST)
  TARGET_LINK_LIBRARIES(vtkPVPythonInterpretor ${PY_STATIC_MODULES_LIST})
ENDIF(PY_STATIC_MODULES_LIST)

TARGET_LINK_LIBRARIES(vtkPVPythonInterpretor
  vtkCommon
  ${VTK_PYTHON_LIBRARIES}
  ${PYTHON_UTIL_LIBRARY_LIB})

# Link against all the kit wrappers.
SET(_QTPYTHON_LINK_LIBS 
#  vtkPVServerManagerPythonD
  vtkCommon 
  vtkFiltering
  vtkIO
  vtkGraphics
  vtkImaging
  vtkPVPythonInterpretor
#  vtkPVServerManager   # dont add this one here, because then cmake doesn't add the automatic dependencies of this lib (if linked static)
  ${VTK_PYTHON_LIBRARIES}
  ${PYTHON_UTIL_LIBRARY_LIB})

IF(VTK_WRAP_TCL)
  SET(_QTPYTHON_LINK_LIBS ${_PVPYTHON_LINK_LIBS}  ${VTK_TK_LIBRARIES})
ENDIF(VTK_WRAP_TCL)

#TARGET_LINK_LIBRARIES(pvpython${PV_EXE_SUFFIX} ${_PVPYTHON_LINK_LIBS})
#TARGET_LINK_LIBRARIES(pvbatch${PV_EXE_SUFFIX} ${_PVPYTHON_LINK_LIBS})

IF(BORLAND)
  WRITE_FILE(${LIBRARY_OUTPUT_PATH}/vtkCommonython.def
    "EXPORTS\ninitvtkCommonPython=_initvtkCommonPython\n")
ENDIF(BORLAND)

# Install vtkPVPythonInterpretor library.
INSTALL(TARGETS vtkPVPythonInterpretor
  RUNTIME DESTINATION bin COMPONENT Runtime
  LIBRARY DESTINATION lib COMPONENT Runtime
  ARCHIVE DESTINATION lib COMPONENT Development)

# Install the pvpython executable.
#SET(PV_EXE_LIST pvpython pvbatch)
#FOREACH(pvexe ${PV_EXE_LIST})
#  IF(NOT PV_INSTALL_NO_RUNTIME)
#      INSTALL(TARGETS ${pvexe}${PV_EXE_SUFFIX} 
#        DESTINATION bin
#        COMPONENT Runtime)
#  ENDIF(NOT PV_INSTALL_NO_RUNTIME)
#ENDFOREACH(pvexe)

# Allow the user to customize their build with some local options
#
INCLUDE (${VTK_BINARY_DIR}/Wrapping/Tcl/LocalUserOptions.cmake OPTIONAL)
INCLUDE (${VTK_SOURCE_DIR}/Wrapping/Tcl/LocalUserOptions.cmake OPTIONAL)

# Add shared link forwarding executables if necessary.
#IF(PV_NEED_SHARED_FORWARD)
#  FOREACH(pvexe ${PV_EXE_LIST})
#    SET(PV_FORWARD_EXE ${pvexe}${PV_EXE_SUFFIX})
#    CONFIGURE_FILE(
#      ${ParaView_SOURCE_DIR}/Servers/Executables/pv-forward.c.in
#      ${CMAKE_CURRENT_BINARY_DIR}/${pvexe}-forward.c
#      @ONLY IMMEDIATE)
#    ADD_EXECUTABLE(${pvexe} ${CMAKE_CURRENT_BINARY_DIR}/${pvexe}-forward.c)
#    ADD_DEPENDENCIES(${pvexe} ${pvexe}${PV_EXE_SUFFIX})
#    IF(NOT PV_INSTALL_NO_RUNTIME)
#      INSTALL(TARGETS ${pvexe}
#        DESTINATION ${PV_INSTALL_BIN_DIR}
#        COMPONENT Runtime)
#    ENDIF(NOT PV_INSTALL_NO_RUNTIME)
#  ENDFOREACH(pvexe)
#ENDIF(PV_NEED_SHARED_FORWARD)


# Any build configured headers should be added to this 
# variable so that they also get installed.
#SET (dynamicHeaders)

# Reduce the number of dirs that get included on moc command line
# since it causes issues on Windows 2000.
GET_DIRECTORY_PROPERTY(include_dirs_tmp INCLUDE_DIRECTORIES)
SET_DIRECTORY_PROPERTIES(PROPERTIES INCLUDE_DIRECTORIES "${MOC_INCLUDE_DIRS}")

QT4_WRAP_CPP(QtPython_MOC_BUILT_SOURCES
  QtPythonDialog.h
  QtPythonShell.h
)

QT4_WRAP_UI(QtPython_UI_BUILT_SOURCES
  QtPythonDialog.ui
)

SET_DIRECTORY_PROPERTIES(PROPERTIES INCLUDE_DIRECTORIES "${include_dirs_tmp}")

ADD_LIBRARY(QtPython
  QtPythonDialog.cxx
  QtPythonDialog.h
  QtPythonShell.cxx
  QtPythonShell.h
  QtPythonExport.h
  ${QtPython_MOC_BUILT_SOURCES}
  ${QtPython_UI_BUILT_SOURCES}
)

TARGET_LINK_LIBRARIES(QtPython
  pqWidgets
  ${QT_LIBRARIES}
  ${PYTHON_LIBRARIES}
  vtkPVPythonInterpretor
  QVTK
)

CONFIGURE_FILE(
  ${QtPython_SOURCE_DIR}/QtPythonConfig.h.in
  ${QtPython_BINARY_DIR}/QtPythonConfig.h)
SET(dynamicHeaders 
    "${dynamicHeaders};${QtPython_BINARY_DIR}/QtPythonConfig.h")

SOURCE_GROUP("Generated Files" FILES
  ${QtPython_MOC_BUILT_SOURCES}
  ${QtPython_UI_BUILT_SOURCES})

# Install
#--------------------------------------------------------------------------------------------------
# Libraries
#IF(NOT PV_INSTALL_NO_LIBRARIES)
#  SET (targets_to_install QtPython)
#  INSTALL(TARGETS ${targets_to_install}
#    RUNTIME DESTINATION ${PV_INSTALL_BIN_DIR} COMPONENT Runtime
#    LIBRARY DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Runtime
#    ARCHIVE DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Development)
#ENDIF(NOT PV_INSTALL_NO_LIBRARIES)
# Development
#IF (NOT PV_INSTALL_NO_DEVELOPMENT)
#  # Headers
#  GLOB_INSTALL_DEVELOPMENT(
#      ${CMAKE_CURRENT_SOURCE_DIR}
#      ${PV_INSTALL_INCLUDE_DIR} 
#      "*.h;*.hxx;*.txx")
#  INSTALL(
#      FILES ${dynamicHeaders}
#      DESTINATION ${PV_INSTALL_INCLUDE_DIR}
#      COMPONENT Development)
  # QT resources
#  GLOB_INSTALL_DEVELOPMENT( 
#    ${CMAKE_CURRENT_SOURCE_DIR}
#    ${PV_INSTALL_INCLUDE_DIR}/Qt/Python/Resources
#    "*.ui;*.ico;*.png;*.svg;*.xpm;*.rc;*.qrc")
#ENDIF (NOT PV_INSTALL_NO_DEVELOPMENT)
