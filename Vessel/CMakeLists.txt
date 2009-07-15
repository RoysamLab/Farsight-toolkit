PROJECT(Vessel)

FIND_PACKAGE(OpenGL)
IF(OPENGL_FOUND)
  INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "Cannot build the Vessel Segmentation module without OpenGL.")
ENDIF()

FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
ELSE()
  MESSAGE(FATAL_ERROR "Cannot build the Vessel Segmentation module without GLUT.")
ENDIF()

INCLUDE_DIRECTORIES(Common)

SUBDIRS(Segmentation)
SUBDIRS(Visualization)
SUBDIRS(Centerline)
