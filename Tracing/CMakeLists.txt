OPTION( BUILD_MDL "Build MDL Tracing" OFF )
IF( BUILD_MDL )
    SUBDIRS( MDL )
ENDIF()

OPTION( BUILD_ftl2d "Build FTL2D Tracing" OFF )
IF( BUILD_ftl2d )
  SUBDIRS( ftl2d )
ENDIF()

OPTION( BUILD_ftl3d "Build FTL3D Tracing" ON )
IF( BUILD_ftl3d )
    SUBDIRS( ftl3d )
ENDIF()

OPTION( BUILD_RPITrace3D "Build RPI Trace3D" OFF )
IF( BUILD_RPITrace3D )
  SUBDIRS( RPITrace3D )
ENDIF()

OPTION( BUILD_TraceEdit "Build Trace Editor" OFF )
IF( BUILD_TraceEdit )
  SUBDIRS( TraceEdit )
ENDIF()

OPTION( BUILD_SpineRing "Build Spine Ring" OFF )
IF( BUILD_SpineRing )
  SUBDIRS( SpineRing )
ENDIF()

OPTION( BUILD_TrackEdit "Build Track Editor" OFF )
IF( BUILD_TrackEdit )
  SUBDIRS( NewTrackEdit )
ENDIF()

OPTION( BUILD_Rendering "Build 5-Label Result Renderer" ON )
IF( BUILD_Rendering )
  SUBDIRS( Rendering )
ENDIF()
