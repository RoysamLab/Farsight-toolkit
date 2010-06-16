OPTION( BUILD_MDL "Build MDL Tracing" ON )

SUBDIRS( Focus )

IF( BUILD_MDL )
    #SUBDIRS( MDL )
    SUBDIRS( MDL2 )
ENDIF()

OPTION( BUILD_MDL_XL "Build MDL_XL" OFF )
IF( BUILD_MDL_XL )
    SUBDIRS( MDL_XL )
ENDIF()

OPTION( BUILD_ftl2d "Build FTL2D Tracing" OFF )
IF( BUILD_ftl2d )
  SUBDIRS( ftl2d )
ENDIF()

OPTION( BUILD_TraceEdit "Build Trace Editor" ON )
IF( BUILD_TraceEdit )
  SUBDIRS( TraceEdit )
ENDIF()

OPTION( BUILD_ftl3d "Build FTL3D Tracing" ON )
IF( BUILD_ftl3d )
    SUBDIRS( ftl3d )
ENDIF()

OPTION( BUILD_RPITrace3D "Build RPI Trace3D" ON )
IF( BUILD_RPITrace3D )
  SUBDIRS( RPITrace3D )
ENDIF()

OPTION( BUILD_SpineRing "Build Spine Ring" ON )
IF( BUILD_SpineRing )
  SUBDIRS( SpineRing )
ENDIF()

#NewTrackEdit is not ready for the use yet.
#OPTION( BUILD_TrackEdit "Build Track Editor" OFF )
#IF( BUILD_TrackEdit )
#  SUBDIRS( NewTrackEdit )
#ENDIF()

OPTION( BUILD_Rendering "Build 5-Label Result Renderer" ON )
IF( BUILD_Rendering )
  SUBDIRS( Rendering )
ENDIF()
