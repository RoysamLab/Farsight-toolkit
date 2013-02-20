PROJECT (RPITrace3D)

SET( RPITrace3D_sources
C3dimage.cpp
Cimage.cpp
Config.cpp
Cpoints.cpp
Ctree.cpp
Cvessel.cpp
draw.cpp
Init.cpp
Main.cpp
Response.cpp
Seeds.cpp
Soma.cpp
StrEle.cpp
Template.cpp
Tools.cpp
Track.cpp
Vector.cpp
Dl.cpp
Branches.cpp
 C3dimage.h
 Cimage.h
 Config.h
 CONSTANTS.h
 Cpoints.h
 Ctree.h
 Cvector.h
 Cvessel.h
 Dllist.h
 Extern.h
 Global.h
 Letters.h
 Main.h
 Mytypes.h
 Numbers.h
 Profile.h
 Queue.h
 Soma.h
 StrEle.h
 Template.h
 timer.h
 Branches.h
)

ADD_EXECUTABLE ( RPITrace3D ${RPITrace3D_sources} )
TARGET_LINK_LIBRARIES( RPITrace3D    TinyXML      ${ITK_LIBRARIES} )
INSTALL(TARGETS RPITrace3D DESTINATION ${INSTALL_BIN_DIR} )

#IF (WIN32)
#	TARGET_LINK_LIBRARIES( RPITrace3D libtiff    TinyXML      ${ITK_LIBRARIES} )
#ELSE (WIN32)
#	TARGET_LINK_LIBRARIES( RPITrace3D tiff    TinyXML      ${ITK_LIBRARIES} )
#ENDIF (WIN32)

ADD_CUSTOM_COMMAND(
  TARGET RPITrace3D
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -D CMAKE_CFG_INTDIR:STRING=${CMAKE_CFG_INTDIR}
  ARGS -D CMAKE_EXECUTABLE_SUFFIX:STRING=${CMAKE_EXECUTABLE_SUFFIX}
  ARGS -D PROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
       -D PROJECT_BINARY_DIR:STRING=${PROJECT_BINARY_DIR}
       -P ${PROJECT_SOURCE_DIR}/ConfigureRPITrace3D.cmake)
