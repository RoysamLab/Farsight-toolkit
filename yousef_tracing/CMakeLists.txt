PROJECT ( RPITrace3D )

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Set ITK_DIR")
ENDIF(ITK_FOUND)

SET( RPITrace3D_sources
 
Branches.cpp
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
 Branches.h
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
)

ADD_EXECUTABLE ( RPITrace3D ${RPITrace3D_sources} )

INSTALL( TARGETS RPITrace3D DESTINATION bin )

IF (WIN32)
	TARGET_LINK_LIBRARIES( RPITrace3D ITKBasicFilters ITKCommon ITKIO )
ELSE (WIN32)
	TARGET_LINK_LIBRARIES( RPITrace3D ITKBasicFilters ITKCommon ITKIO )
ENDIF (WIN32)
