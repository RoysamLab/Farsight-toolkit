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

OPTION( BUILD_Preprocessing "Build image preprocessing applications(only curvelets works as of now)" OFF )
IF( BUILD_Preprocessing )
  SUBDIRS( Preprocessing/CurveletFiltering/fdct_wrapping_cpp/src )
  SUBDIRS( Preprocessing/SHtensor_cpp/2D )
ENDIF()

OPTION( BUILD_Rendering "Build 5-Label Result Renderer" ON )
IF( BUILD_Rendering )
  SUBDIRS( Rendering )
ENDIF()

OPTION( BUILD_ExtractSomas "Build Soma Extraction executables" ${ITK_USE_REVIEW} )
IF( BUILD_ExtractSomas )
  SUBDIRS( ExtractSomas )
ENDIF()

OPTION( BUILD_MultipleNeuronTracer "Build Multiple Neuron Tracer" ON )
IF( BUILD_MultipleNeuronTracer )
  SUBDIRS( MultipleNeuronTracer )
ENDIF()
